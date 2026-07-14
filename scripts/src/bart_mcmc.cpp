// bart_mcmc.cpp — Kim-Rao weighted BART MCMC (one-step random-weight)
// Draws fresh Kim, Rao & Wang (2024) bootstrap weights each MCMC iteration
// and uses them in leaf likelihood / Gibbs updates (raw or K-scaled).
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

using namespace Rcpp;

struct Node {
  bool alive;
  bool is_leaf;
  int depth;
  int parent;     // index into Tree::nodes, -1 if root
  int left;       // index into Tree::nodes, -1 if none
  int right;
  int split_var;  // 0-based column, -1 if leaf
  int split_cut;  // cut index into cuts[split_var], -1 if leaf
  double mu;
  std::vector<int> obs; // 0-based row indices
};

struct Tree {
  std::vector<Node> nodes;
  std::vector<int> free_list;

  int alloc_node() {
    if (!free_list.empty()) {
      int id = free_list.back();
      free_list.pop_back();
      nodes[id].alive = true;
      nodes[id].parent = -1;
      nodes[id].obs.clear();
      return id;
    }
    Node nd;
    nd.alive = true;
    nd.is_leaf = true;
    nd.depth = 0;
    nd.parent = -1;
    nd.left = -1;
    nd.right = -1;
    nd.split_var = -1;
    nd.split_cut = -1;
    nd.mu = 0.0;
    nodes.push_back(nd);
    return static_cast<int>(nodes.size()) - 1;
  }

  void free_node(int id) {
    nodes[id].alive = false;
    nodes[id].is_leaf = true;
    nodes[id].parent = -1;
    nodes[id].left = -1;
    nodes[id].right = -1;
    nodes[id].split_var = -1;
    nodes[id].split_cut = -1;
    nodes[id].obs.clear();
    free_list.push_back(id);
  }
};

// Weighted leaf sufficient statistics: W = Σ w*_i, sum = Σ w*_i R_i
struct Stats {
  double W;
  double sum;
};

struct PsuGroup {
  std::vector<int> idx;
};

struct StratumGroup {
  std::vector<PsuGroup> psus;
};

static inline double runif01() { return R::runif(0.0, 1.0); }
static inline int sample_int(int n) {
  return static_cast<int>(std::floor(runif01() * n));
}

// ---------------------------------------------------------------------------
// Kim-24 design structure + bootstrap weight draw
// ---------------------------------------------------------------------------

static void build_design(
    const IntegerVector& strata,
    const IntegerVector& psu,
    std::vector<StratumGroup>& design
) {
  int n = strata.size();
  // stratum_code -> (psu_code -> member indices)
  std::map<int, std::map<int, std::vector<int> > > tmp;
  for (int i = 0; i < n; ++i)
    tmp[strata[i]][psu[i]].push_back(i);

  design.clear();
  design.reserve(tmp.size());
  for (std::map<int, std::map<int, std::vector<int> > >::iterator it = tmp.begin();
       it != tmp.end(); ++it) {
    StratumGroup sg;
    sg.psus.reserve(it->second.size());
    for (std::map<int, std::vector<int> >::iterator jt = it->second.begin();
         jt != it->second.end(); ++jt) {
      PsuGroup pg;
      pg.idx.swap(jt->second);
      sg.psus.push_back(pg);
    }
    design.push_back(sg);
  }
}

// Kim, Rao & Wang (2024) Example 1:
//   For stratum h with n_h PSUs: Multinomial(n_h-1, uniform), k_h = n_h/(n_h-1),
//   w*_i = K * k_h * m*_hi * w_i  (units in PSU i).
// Optional normalize: rescale so Σ w* = n (Condition A).
static void draw_kim24_weights(
    const NumericVector& w,
    const std::vector<StratumGroup>& design,
    double K,
    bool normalize,
    NumericVector& w_star
) {
  int n = w.size();
  std::fill(w_star.begin(), w_star.end(), 0.0);

  for (size_t h = 0; h < design.size(); ++h) {
    const StratumGroup& sg = design[h];
    int n_h = static_cast<int>(sg.psus.size());
    if (n_h <= 0) continue;

    if (n_h == 1) {
      // Cannot resample; pass through design weights
      const std::vector<int>& idx = sg.psus[0].idx;
      for (size_t j = 0; j < idx.size(); ++j)
        w_star[idx[j]] = K * w[idx[j]];
      continue;
    }

    double k_h = static_cast<double>(n_h) / static_cast<double>(n_h - 1);
    std::vector<double> prob(n_h, 1.0 / n_h);
    std::vector<int> counts(n_h, 0);
    rmultinom(n_h - 1, &prob[0], n_h, &counts[0]);

    for (int i = 0; i < n_h; ++i) {
      double r_star = k_h * static_cast<double>(counts[i]);
      const std::vector<int>& idx = sg.psus[i].idx;
      for (size_t j = 0; j < idx.size(); ++j)
        w_star[idx[j]] = K * r_star * w[idx[j]];
    }
  }

  if (normalize) {
    double Wtot = 0.0;
    for (int i = 0; i < n; ++i) Wtot += w_star[i];
    if (Wtot > 1e-12) {
      double scale = static_cast<double>(n) / Wtot;
      for (int i = 0; i < n; ++i) w_star[i] *= scale;
    } else {
      // Degenerate draw: fall back to design weights normalized to n
      double Wdes = 0.0;
      for (int i = 0; i < n; ++i) Wdes += w[i];
      double scale = (Wdes > 1e-12) ? (static_cast<double>(n) / Wdes) : 1.0;
      for (int i = 0; i < n; ++i) w_star[i] = scale * w[i];
    }
  }
}

static Stats stats_from_obs(
    const NumericVector& r,
    const NumericVector& w_star,
    const std::vector<int>& obs
) {
  Stats s;
  s.W = 0.0;
  s.sum = 0.0;
  for (size_t i = 0; i < obs.size(); ++i) {
    int row = obs[i];
    double wi = w_star[row];
    s.W += wi;
    s.sum += wi * r[row];
  }
  return s;
}

// Chipman integrated leaf log-likelihood with weight sum W in place of n
static double lh(double W, double sy, double sigma, double tau) {
  if (W <= 1e-12) return 0.0;
  double s2 = sigma * sigma;
  double t2 = tau * tau;
  double k = W * t2 + s2;
  return -0.5 * std::log(k) + (t2 * sy * sy) / (2.0 * s2 * k);
}

static double drawnodemu(double W, double sy, double tau, double sigma) {
  if (W <= 1e-12) return R::rnorm(0.0, tau);
  double s2 = sigma * sigma;
  double b = W / s2;
  double a = 1.0 / (tau * tau);
  return (sy / s2) / (a + b) + R::rnorm(0.0, 1.0) / std::sqrt(a + b);
}

// Chipman makexinfo: equally spaced cuts between min and max of each column.
// Split rule: left if x < cut[c]  <=>  Xbin <= c, where
// Xbin = smallest c such that x < cuts[c] (or C if none).
static void build_bins_and_cuts(
    const NumericMatrix& X, int num_cutpoints,
    std::vector<std::vector<double> >& cuts,
    IntegerMatrix& Xbin
) {
  int n = X.nrow();
  int p = X.ncol();
  cuts.assign(p, std::vector<double>());
  Xbin = IntegerMatrix(n, p);

  for (int j = 0; j < p; ++j) {
    double xmin = X(0, j), xmax = X(0, j);
    for (int i = 1; i < n; ++i) {
      double xx = X(i, j);
      if (xx < xmin) xmin = xx;
      if (xx > xmax) xmax = xx;
    }
    std::vector<double> keep;
    keep.reserve(num_cutpoints);
    if (xmax > xmin && num_cutpoints > 0) {
      double xinc = (xmax - xmin) / (num_cutpoints + 1.0);
      for (int k = 0; k < num_cutpoints; ++k)
        keep.push_back(xmin + (k + 1) * xinc);
    }
    cuts[j] = keep;
    int C = static_cast<int>(keep.size());
    for (int i = 0; i < n; ++i) {
      int b = 0;
      while (b < C && !(X(i, j) < keep[b])) ++b;
      Xbin(i, j) = b;
    }
  }
}

// Chipman tree::rg — restrict cut index range for variable v using ancestors
static void rg(const Tree& tree, int node_id, int v, int* L, int* U) {
  int pid = tree.nodes[node_id].parent;
  if (pid < 0) return;
  const Node& par = tree.nodes[pid];
  if (par.split_var == v) {
    if (par.left == node_id) {
      if (par.split_cut <= *U) *U = par.split_cut - 1;
    } else {
      if (par.split_cut >= *L) *L = par.split_cut + 1;
    }
  }
  rg(tree, pid, v, L, U);
}

// Chipman cansplit: some variable still has U >= L in the cut grid
static bool cansplit_node(
    const Tree& tree, int node_id,
    const std::vector<std::vector<double> >& cuts
) {
  int p = static_cast<int>(cuts.size());
  for (int v = 0; v < p; ++v) {
    int C = static_cast<int>(cuts[v].size());
    if (C <= 0) continue;
    int L = 0, U = C - 1;
    rg(tree, node_id, v, &L, &U);
    if (U >= L) return true;
  }
  return false;
}

static void get_goodvars(
    const Tree& tree, int node_id,
    const std::vector<std::vector<double> >& cuts,
    std::vector<int>& goodvars
) {
  goodvars.clear();
  int p = static_cast<int>(cuts.size());
  for (int v = 0; v < p; ++v) {
    int C = static_cast<int>(cuts[v].size());
    if (C <= 0) continue;
    int L = 0, U = C - 1;
    rg(tree, node_id, v, &L, &U);
    if (U >= L) goodvars.push_back(v);
  }
}

// Return cut indices c that leave >= min_leaf on each side (count-based check)
static void valid_cut_ids(
    const IntegerMatrix& Xbin, const std::vector<int>& obs, int j, int C,
    int min_leaf, std::vector<int>& out, std::vector<int>& hist
) {
  out.clear();
  if (C <= 0 || obs.empty()) return;
  if (static_cast<int>(hist.size()) < static_cast<size_t>(C + 1)) hist.assign(C + 1, 0);
  else std::fill(hist.begin(), hist.begin() + C + 1, 0);

  for (size_t i = 0; i < obs.size(); ++i) ++hist[Xbin(obs[i], j)];

  int n_left = 0;
  int n_tot = static_cast<int>(obs.size());
  for (int c = 0; c < C; ++c) {
    n_left += hist[c];
    int n_right = n_tot - n_left;
    if (n_left >= min_leaf && n_right >= min_leaf) out.push_back(c);
  }
}

static void predict_tree_into(const Tree& tree, NumericVector& pred) {
  std::fill(pred.begin(), pred.end(), 0.0);
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    const Node& nd = tree.nodes[i];
    if (!nd.alive || !nd.is_leaf) continue;
    for (size_t k = 0; k < nd.obs.size(); ++k) pred[nd.obs[k]] = nd.mu;
  }
}

// Root is the unique alive node with parent < 0.
static int find_root(const Tree& tree) {
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    if (tree.nodes[i].alive && tree.nodes[i].parent < 0)
      return static_cast<int>(i);
  }
  return 0;
}

// OOS walk: left if x < cuts[v][split_cut]  (same rule as Xbin <= split_cut).
static double predict_tree_row(
    const Tree& tree, const NumericMatrix& X, int row,
    const std::vector<std::vector<double> >& cuts
) {
  int nid = find_root(tree);
  while (true) {
    const Node& nd = tree.nodes[nid];
    if (!nd.alive) stop("dead node encountered during prediction");
    if (nd.is_leaf) return nd.mu;
    if (nd.split_var < 0 || nd.split_cut < 0 || nd.left < 0 || nd.right < 0)
      stop("invalid split during prediction");
    double xv = X(row, nd.split_var);
    if (xv < cuts[nd.split_var][nd.split_cut]) nid = nd.left;
    else nid = nd.right;
  }
}

static double predict_ensemble_row(
    const std::vector<Tree>& ensemble, const NumericMatrix& X, int row,
    const std::vector<std::vector<double> >& cuts
) {
  double f = 0.0;
  for (size_t t = 0; t < ensemble.size(); ++t)
    f += predict_tree_row(ensemble[t], X, row, cuts);
  return f;
}

static void find_nogs(const Tree& tree, std::vector<int>& nogs) {
  nogs.clear();
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    const Node& nd = tree.nodes[i];
    if (!nd.alive || nd.is_leaf) continue;
    if (nd.left < 0 || nd.right < 0) continue;
    if (tree.nodes[nd.left].alive && tree.nodes[nd.left].is_leaf &&
        tree.nodes[nd.right].alive && tree.nodes[nd.right].is_leaf) {
      nogs.push_back(static_cast<int>(i));
    }
  }
}

static Tree make_root_tree(int n) {
  Tree t;
  int id = t.alloc_node();
  Node& root = t.nodes[id];
  root.is_leaf = true;
  root.depth = 0;
  root.parent = -1;
  root.mu = 0.0;
  root.obs.resize(n);
  for (int i = 0; i < n; ++i) root.obs[i] = i;
  return t;
}

static void update_leaf_mus(
    Tree& tree, const NumericVector& r_t, const NumericVector& w_star,
    double sigma2, double sigma_mu
) {
  double tau = sigma_mu * sigma_mu;
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    Node& nd = tree.nodes[i];
    if (!nd.alive || !nd.is_leaf) continue;
    Stats st = stats_from_obs(r_t, w_star, nd.obs);
    double post_var, post_mean;
    if (st.W <= 1e-12) {
      post_var = tau;
      post_mean = 0.0;
    } else {
      post_var = 1.0 / (st.W / sigma2 + 1.0 / tau);
      post_mean = post_var * (st.sum / sigma2);
    }
    nd.mu = R::rnorm(post_mean, std::sqrt(post_var));
  }
}

static bool mh_grow(
    Tree& tree, const NumericVector& r_t, const NumericVector& w_star,
    const IntegerMatrix& Xbin,
    const std::vector<std::vector<double> >& cuts,
    double sigma2, double sigma_mu, double alpha, double beta,
    double p_grow, double pb_default, int min_leaf,
    std::vector<int>& growable, std::vector<int>& goodvars
) {
  growable.clear();
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    Node& nd = tree.nodes[i];
    if (!nd.alive || !nd.is_leaf) continue;
    if (cansplit_node(tree, static_cast<int>(i), cuts))
      growable.push_back(static_cast<int>(i));
  }
  int n_grow = static_cast<int>(growable.size());
  if (n_grow == 0) return false;

  int parent_pos = growable[sample_int(n_grow)];
  int parent_depth = tree.nodes[parent_pos].depth;
  double parent_mu = tree.nodes[parent_pos].mu;
  // copy obs — may realloc on alloc_node
  std::vector<int> parent_obs = tree.nodes[parent_pos].obs;

  get_goodvars(tree, parent_pos, cuts, goodvars);
  if (goodvars.empty()) return false;

  int split_var = goodvars[sample_int(static_cast<int>(goodvars.size()))];
  int C = static_cast<int>(cuts[split_var].size());
  int L = 0, U = C - 1;
  rg(tree, parent_pos, split_var, &L, &U);
  if (U < L) return false;
  int split_cut = L + sample_int(U - L + 1);

  // Chipman: propose any (L,U) cut; reject if nl or nr < min_leaf (count-based)
  std::vector<int> left_obs, right_obs;
  left_obs.reserve(parent_obs.size());
  right_obs.reserve(parent_obs.size());
  for (size_t i = 0; i < parent_obs.size(); ++i) {
    int r = parent_obs[i];
    if (Xbin(r, split_var) <= split_cut) left_obs.push_back(r);
    else right_obs.push_back(r);
  }
  if (static_cast<int>(left_obs.size()) < min_leaf ||
      static_cast<int>(right_obs.size()) < min_leaf) {
    return false;
  }

  Stats parent_stats = stats_from_obs(r_t, w_star, parent_obs);
  Stats left_stats = stats_from_obs(r_t, w_star, left_obs);
  Stats right_stats = stats_from_obs(r_t, w_star, right_obs);

  Node& parent0 = tree.nodes[parent_pos];
  bool old_is_leaf = parent0.is_leaf;
  int old_left = parent0.left;
  int old_right = parent0.right;
  int old_svar = parent0.split_var;
  int old_scut = parent0.split_cut;
  double old_mu = parent0.mu;
  int old_parent = parent0.parent;

  int left_id = tree.alloc_node();
  int right_id = tree.alloc_node();
  Node& parent = tree.nodes[parent_pos];
  Node& left = tree.nodes[left_id];
  Node& right = tree.nodes[right_id];

  left.is_leaf = true;
  left.depth = parent_depth + 1;
  left.parent = parent_pos;
  left.mu = parent_mu;
  left.obs.swap(left_obs);

  right.is_leaf = true;
  right.depth = parent_depth + 1;
  right.parent = parent_pos;
  right.mu = parent_mu;
  right.obs.swap(right_obs);

  parent.is_leaf = false;
  parent.split_var = split_var;
  parent.split_cut = split_cut;
  parent.left = left_id;
  parent.right = right_id;

  int d = parent_depth;
  double sigma = std::sqrt(sigma2);
  double tau = sigma_mu;

  // prior grow probs for new children (Chipman: may be 0 if v exhausted)
  double pg_child = alpha * std::pow(1.0 + d + 1.0, -beta);
  double PGly, PGry;
  if (goodvars.size() > 1) {
    PGly = pg_child;
    PGry = pg_child;
  } else {
    PGly = ((split_cut - 1) < L) ? 0.0 : pg_child;
    PGry = (U < (split_cut + 1)) ? 0.0 : pg_child;
  }

  double PDy;
  if (n_grow > 1) {
    PDy = 1.0 - pb_default;
  } else {
    PDy = (PGly == 0.0 && PGry == 0.0) ? 1.0 : (1.0 - pb_default);
  }

  std::vector<int> nogs;
  find_nogs(tree, nogs);
  int n_prune_new = static_cast<int>(nogs.size());
  if (n_prune_new < 1) {
    parent.is_leaf = old_is_leaf;
    parent.parent = old_parent;
    parent.left = old_left;
    parent.right = old_right;
    parent.split_var = old_svar;
    parent.split_cut = old_scut;
    parent.mu = old_mu;
    tree.free_node(left_id);
    tree.free_node(right_id);
    return false;
  }

  double PGnx = alpha * std::pow(1.0 + d, -beta);
  double log_lik =
      lh(left_stats.W, left_stats.sum, sigma, tau) +
      lh(right_stats.W, right_stats.sum, sigma, tau) -
      lh(parent_stats.W, parent_stats.sum, sigma, tau) +
      std::log(sigma);
  // Chipman pr = PGnx*(1-PGly)*(1-PGry)*PDy*Pnogy / ((1-PGnx)*Pbotx*PBx)
  double log_pr =
      std::log(PGnx) + std::log(1.0 - PGly) + std::log(1.0 - PGry) +
      std::log(PDy) - std::log(static_cast<double>(n_prune_new)) -
      std::log(1.0 - PGnx) - std::log(static_cast<double>(n_grow)) -
      std::log(p_grow);

  double log_alpha = log_lik + log_pr;
  if (R_finite(log_alpha) && std::log(runif01()) < log_alpha) {
    left.mu = drawnodemu(left_stats.W, left_stats.sum, tau, sigma);
    right.mu = drawnodemu(right_stats.W, right_stats.sum, tau, sigma);
    return true;
  }

  parent.is_leaf = old_is_leaf;
  parent.parent = old_parent;
  parent.left = old_left;
  parent.right = old_right;
  parent.split_var = old_svar;
  parent.split_cut = old_scut;
  parent.mu = old_mu;
  tree.free_node(left_id);
  tree.free_node(right_id);
  return false;
}

static bool mh_prune(
    Tree& tree, const NumericVector& r_t, const NumericVector& w_star,
    const std::vector<std::vector<double> >& cuts,
    double sigma2, double sigma_mu, double alpha, double beta,
    double pb_default, double p_prune,
    std::vector<int>& nogs, std::vector<int>& growable
) {
  find_nogs(tree, nogs);
  int n_prune = static_cast<int>(nogs.size());
  if (n_prune == 0) return false;

  // goodbots on current tree (before death)
  growable.clear();
  for (size_t i = 0; i < tree.nodes.size(); ++i) {
    if (!tree.nodes[i].alive || !tree.nodes[i].is_leaf) continue;
    if (cansplit_node(tree, static_cast<int>(i), cuts))
      growable.push_back(static_cast<int>(i));
  }
  int n_good_before = static_cast<int>(growable.size());

  int parent_pos = nogs[sample_int(n_prune)];
  Node& parent = tree.nodes[parent_pos];
  int left_id = parent.left;
  int right_id = parent.right;
  Node& left = tree.nodes[left_id];
  Node& right = tree.nodes[right_id];

  int old_svar = parent.split_var;
  int old_scut = parent.split_cut;
  int d = parent.depth;

  bool left_cansplit = cansplit_node(tree, left_id, cuts);
  bool right_cansplit = cansplit_node(tree, right_id, cuts);

  Stats parent_stats = stats_from_obs(r_t, w_star, parent.obs);
  Stats left_stats = stats_from_obs(r_t, w_star, left.obs);
  Stats right_stats = stats_from_obs(r_t, w_star, right.obs);

  Node left_bak = left;
  Node right_bak = right;
  bool old_is_leaf = parent.is_leaf;
  int old_left = parent.left;
  int old_right = parent.right;
  double old_mu = parent.mu;

  double Wl = left_stats.W;
  double Wr = right_stats.W;
  double Wtot = Wl + Wr;
  double wmu = (Wtot > 1e-12)
      ? (Wl * left.mu + Wr * right.mu) / Wtot
      : 0.5 * (left.mu + right.mu);

  parent.is_leaf = true;
  parent.split_var = -1;
  parent.split_cut = -1;
  parent.left = -1;
  parent.right = -1;
  parent.mu = wmu;
  tree.free_node(left_id);
  tree.free_node(right_id);

  // Chipman Pboty = 1 / (goodbots - cansplit(L) - cansplit(R) + 1)
  int ngood = n_good_before;
  if (left_cansplit) --ngood;
  if (right_cansplit) --ngood;
  ++ngood;
  if (ngood < 1) ngood = 1;

  double PBy = (parent.parent < 0) ? 1.0 : pb_default;

  double sigma = std::sqrt(sigma2);
  double tau = sigma_mu;
  double PGny = alpha * std::pow(1.0 + d, -beta);
  double PGlx = left_cansplit ? (alpha * std::pow(1.0 + d + 1.0, -beta)) : 0.0;
  double PGrx = right_cansplit ? (alpha * std::pow(1.0 + d + 1.0, -beta)) : 0.0;

  double log_lik =
      lh(parent_stats.W, parent_stats.sum, sigma, tau) -
      lh(left_stats.W, left_stats.sum, sigma, tau) -
      lh(right_stats.W, right_stats.sum, sigma, tau) -
      std::log(sigma);
  // Chipman dprop: pr = ((1-PGny)*PBy*Pboty) / (PGny*(1-PGlx)*(1-PGrx)*PDx*(1/nnogs))
  double log_pr =
      std::log(1.0 - PGny) + std::log(PBy) - std::log(static_cast<double>(ngood)) -
      std::log(PGny) - std::log(1.0 - PGlx) - std::log(1.0 - PGrx) -
      std::log(p_prune) + std::log(static_cast<double>(n_prune));

  double log_alpha = log_lik + log_pr;
  if (R_finite(log_alpha) && std::log(runif01()) < log_alpha) {
    parent.mu = drawnodemu(parent_stats.W, parent_stats.sum, tau, sigma);
    return true;
  }

  // undo prune
  parent.is_leaf = old_is_leaf;
  parent.split_var = old_svar;
  parent.split_cut = old_scut;
  parent.left = old_left;
  parent.right = old_right;
  parent.mu = old_mu;
  std::vector<int> fl;
  fl.reserve(tree.free_list.size());
  for (size_t i = 0; i < tree.free_list.size(); ++i) {
    if (tree.free_list[i] != left_id && tree.free_list[i] != right_id)
      fl.push_back(tree.free_list[i]);
  }
  tree.free_list.swap(fl);
  tree.nodes[left_id] = left_bak;
  tree.nodes[right_id] = right_bak;
  tree.nodes[left_id].alive = true;
  tree.nodes[right_id].alive = true;
  return false;
}

static bool mh_change(
    Tree& tree, const NumericVector& r_t, const NumericVector& w_star,
    const IntegerMatrix& Xbin,
    const std::vector<std::vector<double> >& cuts,
    double sigma2, double sigma_mu, int min_leaf,
    std::vector<int>& nogs, std::vector<int>& tmp_cuts, std::vector<int>& hist,
    std::vector<int>& valid_vars, std::vector<std::vector<int> >& valid_cut_store
) {
  find_nogs(tree, nogs);
  if (nogs.empty()) return false;

  int parent_pos = nogs[sample_int(static_cast<int>(nogs.size()))];
  Node& node = tree.nodes[parent_pos];
  Node& left = tree.nodes[node.left];
  Node& right = tree.nodes[node.right];
  int p = Xbin.ncol();

  valid_vars.clear();
  valid_cut_store.clear();
  for (int j = 0; j < p; ++j) {
    valid_cut_ids(Xbin, node.obs, j, static_cast<int>(cuts[j].size()), min_leaf, tmp_cuts, hist);
    if (!tmp_cuts.empty()) {
      valid_vars.push_back(j);
      valid_cut_store.push_back(tmp_cuts);
    }
  }
  if (valid_vars.empty()) return false;

  int vi = sample_int(static_cast<int>(valid_vars.size()));
  int split_var = valid_vars[vi];
  int split_cut = valid_cut_store[vi][sample_int(static_cast<int>(valid_cut_store[vi].size()))];

  Stats old_left = stats_from_obs(r_t, w_star, left.obs);
  Stats old_right = stats_from_obs(r_t, w_star, right.obs);

  std::vector<int> bak_left = left.obs;
  std::vector<int> bak_right = right.obs;
  int bak_svar = node.split_var;
  int bak_scut = node.split_cut;

  left.obs.clear();
  right.obs.clear();
  left.obs.reserve(node.obs.size());
  right.obs.reserve(node.obs.size());
  for (size_t i = 0; i < node.obs.size(); ++i) {
    int r = node.obs[i];
    if (Xbin(r, split_var) <= split_cut) left.obs.push_back(r);
    else right.obs.push_back(r);
  }
  node.split_var = split_var;
  node.split_cut = split_cut;

  Stats new_left = stats_from_obs(r_t, w_star, left.obs);
  Stats new_right = stats_from_obs(r_t, w_star, right.obs);

  double sigma = std::sqrt(sigma2);
  double tau = sigma_mu;
  double log_lik =
      lh(new_left.W, new_left.sum, sigma, tau) +
      lh(new_right.W, new_right.sum, sigma, tau) -
      lh(old_left.W, old_left.sum, sigma, tau) -
      lh(old_right.W, old_right.sum, sigma, tau);

  if (R_finite(log_lik) && std::log(runif01()) < log_lik) {
    left.mu = drawnodemu(new_left.W, new_left.sum, tau, sigma);
    right.mu = drawnodemu(new_right.W, new_right.sum, tau, sigma);
    return true;
  }

  // undo
  node.split_var = bak_svar;
  node.split_cut = bak_scut;
  left.obs.swap(bak_left);
  right.obs.swap(bak_right);
  return false;
}

//' Run Kim-Rao weighted BART MCMC (one-step random-weight)
//'
//' Draws fresh Kim, Rao & Wang (2024) bootstrap weights each iteration and
//' uses them in leaf MH ratios and Gibbs updates. Tree prior / proposals
//' remain Chipman (unweighted). Leaf min size stays count-based.
//'
//' @param weights Design weights w_i (length n).
//' @param strata Stratum id per unit (any integers).
//' @param psu PSU id per unit within strata (any integers).
//' @param K Multiplicative scale for bootstrap weights (1 = raw Kim-Rao).
//' @param normalize_weights If TRUE, rescale each draw so sum(w*) = n.
//' @param X_test Optional test covariates (same p as X) for OOS predictions.
//' @export
// [[Rcpp::export]]
List run_bart_mcmc_cpp(
    NumericMatrix X,
    NumericVector y_scaled,
    NumericVector weights,
    IntegerVector strata,
    IntegerVector psu,
    int num_trees = 50,
    int num_iter = 500,
    int burn_in = 100,
    int min_leaf = 5,
    int num_cutpoints = 100,
    double alpha = 0.95,
    double beta = 2.0,
    double sigma_mu = 0.05,
    double sigma2_init = 0.1,
    double nu = 3.0,
    double lambda = 0.1,
    double p_grow = 0.28,
    double p_prune = 0.28,
    double p_change = 0.44,
    int verbose_every = 0,
    double y_center = 0.0,
    double y_range = 1.0,
    double K = 1.0,
    bool normalize_weights = false,
    Nullable<NumericMatrix> X_test = R_NilValue,
    bool return_weights = false
) {
  int n = X.nrow();
  if (weights.size() != n || strata.size() != n || psu.size() != n)
    stop("weights, strata, and psu must each have length nrow(X)");
  if (K <= 0.0) stop("K must be positive");

  bool has_test = X_test.isNotNull();
  NumericMatrix X_test_mat;
  int n_test = 0;
  if (has_test) {
    X_test_mat = X_test.get();
    n_test = X_test_mat.nrow();
    if (X_test_mat.ncol() != X.ncol())
      stop("X_test must have the same number of columns as X");
  }

  std::vector<std::vector<double> > cuts;
  IntegerMatrix Xbin;
  build_bins_and_cuts(X, num_cutpoints, cuts, Xbin);

  std::vector<StratumGroup> design;
  build_design(strata, psu, design);

  std::vector<Tree> ensemble(num_trees);
  for (int t = 0; t < num_trees; ++t) ensemble[t] = make_root_tree(n);

  double sigma2 = sigma2_init;
  int n_keep = std::max(num_iter - burn_in, 0);
  NumericMatrix pred_samples(n_keep, n);
  NumericMatrix pred_test_samples(n_keep, has_test ? n_test : 0);
  NumericMatrix weight_samples(return_weights ? n_keep : 0, return_weights ? n : 0);
  NumericVector sigma2_samples(n_keep);

  NumericMatrix tree_preds(n, num_trees);
  NumericVector total_pred(n);
  NumericVector tmp_pred(n);
  for (int t = 0; t < num_trees; ++t) {
    predict_tree_into(ensemble[t], tmp_pred);
    for (int i = 0; i < n; ++i) {
      tree_preds(i, t) = tmp_pred[i];
      total_pred[i] += tmp_pred[i];
    }
  }

  int att_grow = 0, acc_grow = 0;
  int att_prune = 0, acc_prune = 0;

  NumericVector r_t(n);
  NumericVector w_star(n);
  // Chipman defaults: birth/death only; pb = P(birth | can birth)
  const double pb_default = 0.5;
  (void)p_grow;
  (void)p_prune;
  (void)p_change;

  // reusable scratch buffers
  std::vector<int> growable, nogs, goodvars;

  for (int iter = 1; iter <= num_iter; ++iter) {
    // Fresh Kim-24 bootstrap weights for this iteration
    draw_kim24_weights(weights, design, K, normalize_weights, w_star);

    for (int t = 0; t < num_trees; ++t) {
      for (int i = 0; i < n; ++i) r_t[i] = y_scaled[i] - (total_pred[i] - tree_preds(i, t));

      // Adaptive P(birth) like Chipman getpb() using cansplit (ancestor L,U)
      growable.clear();
      for (size_t i = 0; i < ensemble[t].nodes.size(); ++i) {
        Node& nd = ensemble[t].nodes[i];
        if (!nd.alive || !nd.is_leaf) continue;
        if (cansplit_node(ensemble[t], static_cast<int>(i), cuts))
          growable.push_back(static_cast<int>(i));
      }
      find_nogs(ensemble[t], nogs);
      double pb_x;
      if (growable.empty()) pb_x = 0.0;
      else if (nogs.empty()) pb_x = 1.0; // stump
      else pb_x = pb_default;

      bool accepted = false;
      if (runif01() < pb_x) {
        ++att_grow;
        accepted = mh_grow(ensemble[t], r_t, w_star, Xbin, cuts, sigma2, sigma_mu,
                           alpha, beta, pb_x, pb_default, min_leaf,
                           growable, goodvars);
        if (accepted) ++acc_grow;
      } else {
        ++att_prune;
        double pd_x = 1.0 - pb_x;
        accepted = mh_prune(ensemble[t], r_t, w_star, cuts, sigma2, sigma_mu,
                            alpha, beta, pb_default, pd_x,
                            nogs, growable);
        if (accepted) ++acc_prune;
      }
      (void)accepted;

      // Chipman: always redraw all bottom-node mus after birth/death
      update_leaf_mus(ensemble[t], r_t, w_star, sigma2, sigma_mu);
      predict_tree_into(ensemble[t], tmp_pred);
      for (int i = 0; i < n; ++i) {
        total_pred[i] += tmp_pred[i] - tree_preds(i, t);
        tree_preds(i, t) = tmp_pred[i];
      }
    }

    // sigma^2: n-scaled weighted RSS for numerical stability (same w* draw)
    double Wtot = 0.0;
    for (int i = 0; i < n; ++i) Wtot += w_star[i];
    double ss = 0.0;
    if (Wtot > 1e-12) {
      double scale = static_cast<double>(n) / Wtot;
      for (int i = 0; i < n; ++i) {
        double e = y_scaled[i] - total_pred[i];
        ss += scale * w_star[i] * e * e;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        double e = y_scaled[i] - total_pred[i];
        ss += e * e;
      }
    }
    double shape = (nu + n) / 2.0;
    double rate = (nu * lambda + ss) / 2.0;
    sigma2 = 1.0 / R::rgamma(shape, 1.0 / rate);

    if (iter > burn_in) {
      int row = iter - burn_in - 1;
      for (int i = 0; i < n; ++i) pred_samples(row, i) = total_pred[i] * y_range + y_center;
      sigma2_samples[row] = sigma2 * y_range * y_range;
      if (return_weights) {
        for (int i = 0; i < n; ++i) weight_samples(row, i) = w_star[i];
      }
      if (has_test) {
        for (int j = 0; j < n_test; ++j) {
          double f = predict_ensemble_row(ensemble, X_test_mat, j, cuts);
          pred_test_samples(row, j) = f * y_range + y_center;
        }
      }
    }

    if (verbose_every > 0 && iter % verbose_every == 0) {
      Rprintf("iter %d / %d | sigma2(scaled)=%.5f\n", iter, num_iter, sigma2);
    }
  }

  NumericVector accept_rates = NumericVector::create(
      _["grow"] = att_grow ? static_cast<double>(acc_grow) / att_grow : NA_REAL,
      _["prune"] = att_prune ? static_cast<double>(acc_prune) / att_prune : NA_REAL
  );

  List out = List::create(
      _["pred_samples"] = pred_samples,
      _["sigma2_samples"] = sigma2_samples,
      _["accept_rates"] = accept_rates,
      _["y_center"] = y_center,
      _["y_range"] = y_range,
      _["K"] = K,
      _["normalize_weights"] = normalize_weights
  );
  if (has_test) out["pred_test_samples"] = pred_test_samples;
  else out["pred_test_samples"] = R_NilValue;
  if (return_weights) out["weight_samples"] = weight_samples;
  else out["weight_samples"] = R_NilValue;
  return out;
}
