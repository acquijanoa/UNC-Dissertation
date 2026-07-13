// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <random>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

std::mt19937 engine(42);

// Truncated Normal Sampler (Inverse CDF)
double r_trunc_norm(double mu, double lower, double upper) {
    double p_lower = R::pnorm(lower, mu, 1.0, 1, 0);
    double p_upper = R::pnorm(upper, mu, 1.0, 1, 0);
    double u = R::runif(0.0, 1.0);
    double p = p_lower + u * (p_upper - p_lower);
    if (p >= 1.0) p = 1.0 - 1e-15;
    if (p <= 0.0) p = 1e-15;
    return R::qnorm(p, mu, 1.0, 1, 0);
}

// ---------------------------------------------------------
// Node Structure
// ---------------------------------------------------------
struct Node {
    bool is_leaf;
    int feature;
    double threshold;
    double mu;
    int depth;
    int node_id;
    Node* left;
    Node* right;

    Node(int d = 0, int id = 0) : is_leaf(true), feature(-1), threshold(0.0), 
                                  mu(0.0), depth(d), node_id(id), left(nullptr), right(nullptr) {}

    ~Node() {
        delete left;
        delete right;
    }
    
    double predict(const rowvec& x) const {
        if (is_leaf) return mu;
        if (x(feature) <= threshold) {
            return left->predict(x);
        } else {
            return right->predict(x);
        }
    }
    
    void get_leaves(std::vector<Node*>& leaves) {
        if (is_leaf) {
            leaves.push_back(this);
        } else {
            left->get_leaves(leaves);
            right->get_leaves(leaves);
        }
    }
    
    void get_internals(std::vector<Node*>& internals) {
        if (!is_leaf) {
            internals.push_back(this);
            left->get_internals(internals);
            right->get_internals(internals);
        }
    }
    
    void route(const mat& X, const std::vector<int>& idx, std::vector<std::vector<int>>& leaf_idx) {
        if (is_leaf) {
            leaf_idx[node_id] = idx;
        } else {
            std::vector<int> left_idx, right_idx;
            for (int i : idx) {
                if (X(i, feature) <= threshold) left_idx.push_back(i);
                else right_idx.push_back(i);
            }
            if (!left_idx.empty()) left->route(X, left_idx, leaf_idx);
            if (!right_idx.empty()) right->route(X, right_idx, leaf_idx);
        }
    }
};

// ---------------------------------------------------------
// Prior config
// ---------------------------------------------------------
struct Prior {
    double alpha;
    double beta;
    double k;
    double sigma_mu;
    double nu;
    double lam;
};

// ---------------------------------------------------------
// Tree Class
// ---------------------------------------------------------
class Tree {
public:
    Node* root;
    int node_counter;
    Prior prior;

    Tree(Prior p) : prior(p) {
        node_counter = 0;
        root = new Node(0, next_id());
    }
    
    ~Tree() {
        delete root;
    }

    int next_id() { return node_counter++; }

    vec predict(const mat& X) const {
        int n = X.n_rows;
        vec preds(n);
        for(int i=0; i<n; ++i) {
            preds(i) = root->predict(X.row(i));
        }
        return preds;
    }
    
    double weighted_leaf_ll(const vec& R_l, const vec& w_l, double sigma2) {
        if (R_l.is_empty()) return 0.0;
        double W = sum(w_l);
        if (W <= 1e-8) return 0.0;
        double R_bar = sum(w_l % R_l) / W;
        double V = 1.0 / (W / sigma2 + 1.0 / std::pow(prior.sigma_mu, 2));
        double mu_hat = V * W * R_bar / sigma2;
        
        double ll = -0.5 * std::log(sigma2) * R_l.n_elem;
        for (uword i=0; i<R_l.n_elem; ++i) {
            ll -= 0.5 * w_l(i) * std::pow(R_l(i) - mu_hat, 2) / sigma2;
        }
        return ll;
    }
    
    void grow_proposal(const mat& X, const vec& R, const vec& w_star, double sigma2, std::vector<std::vector<int>>& leaf_idx) {
        std::vector<Node*> leaves;
        root->get_leaves(leaves);
        
        std::vector<Node*> eligible;
        for (Node* l : leaves) {
            double eff_n = 0;
            for (int i : leaf_idx[l->node_id]) if (w_star(i) > 0) eff_n++;
            if (eff_n >= 2) eligible.push_back(l);
        }
        if (eligible.empty()) return;
        
        std::uniform_int_distribution<int> unif_l(0, eligible.size() - 1);
        Node* leaf = eligible[unif_l(engine)];
        const std::vector<int>& idx = leaf_idx[leaf->node_id];
        
        std::uniform_int_distribution<int> unif_feat(0, X.n_cols - 1);
        int feat = unif_feat(engine);
        
        std::vector<double> vals;
        for (int i : idx) if (w_star(i) > 0) vals.push_back(X(i, feat));
        if (vals.size() < 2) return;
        
        std::sort(vals.begin(), vals.end());
        auto last = std::unique(vals.begin(), vals.end());
        vals.erase(last, vals.end());
        if (vals.size() < 2) return;
        
        std::uniform_int_distribution<int> unif_val(0, vals.size() - 2);
        double threshold = vals[unif_val(engine)];
        
        vec R_l(idx.size()), w_l(idx.size());
        for(size_t i=0; i<idx.size(); ++i) { R_l(i) = R(idx[i]); w_l(i) = w_star(idx[i]); }
        double ll_before = weighted_leaf_ll(R_l, w_l, sigma2);
        
        std::vector<double> R_left, w_left, R_right, w_right;
        for (int i : idx) {
            if (X(i, feat) <= threshold) {
                R_left.push_back(R(i)); w_left.push_back(w_star(i));
            } else {
                R_right.push_back(R(i)); w_right.push_back(w_star(i));
            }
        }
        
        if (R_left.empty() || R_right.empty()) return;
        
        double ll_left = weighted_leaf_ll(vec(R_left), vec(w_left), sigma2);
        double ll_right = weighted_leaf_ll(vec(R_right), vec(w_right), sigma2);
        double ll_after = ll_left + ll_right;
        
        double depth = leaf->depth;
        double prior_grow = prior.alpha / std::pow(1.0 + depth, prior.beta);
        double prior_ratio = std::log(prior_grow) + 2.0 * std::log(1.0 - prior.alpha / std::pow(1.0 + depth + 1.0, prior.beta)) - std::log(1.0 - prior_grow);
        
        double log_accept = ll_after - ll_before + prior_ratio;
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        if (std::log(unif(engine)) < log_accept) {
            leaf->is_leaf = false;
            leaf->feature = feat;
            leaf->threshold = threshold;
            leaf->left = new Node(leaf->depth + 1, next_id());
            leaf->right = new Node(leaf->depth + 1, next_id());
        }
    }
    
    void prune_proposal(const mat& X, const vec& R, const vec& w_star, double sigma2, std::vector<std::vector<int>>& leaf_idx) {
        std::vector<Node*> internals;
        root->get_internals(internals);
        
        std::vector<Node*> prunable;
        for (Node* n : internals) {
            if (n->left->is_leaf && n->right->is_leaf) prunable.push_back(n);
        }
        if (prunable.empty()) return;
        
        std::uniform_int_distribution<int> unif_n(0, prunable.size() - 1);
        Node* node = prunable[unif_n(engine)];
        const std::vector<int>& idx = leaf_idx[node->node_id];
        
        vec R_n(idx.size()), w_n(idx.size());
        for(size_t i=0; i<idx.size(); ++i) { R_n(i) = R(idx[i]); w_n(i) = w_star(idx[i]); }
        
        double ll_after = weighted_leaf_ll(R_n, w_n, sigma2);
        
        std::vector<double> R_left, w_left, R_right, w_right;
        for (int i : idx) {
            if (X(i, node->feature) <= node->threshold) {
                R_left.push_back(R(i)); w_left.push_back(w_star(i));
            } else {
                R_right.push_back(R(i)); w_right.push_back(w_star(i));
            }
        }
        
        double ll_left = weighted_leaf_ll(vec(R_left), vec(w_left), sigma2);
        double ll_right = weighted_leaf_ll(vec(R_right), vec(w_right), sigma2);
        double ll_before = ll_left + ll_right;
        
        double depth = node->depth;
        double prior_grow = prior.alpha / std::pow(1.0 + depth, prior.beta);
        double prior_ratio = std::log(1.0 - prior_grow) - std::log(prior_grow) - 2.0 * std::log(1.0 - prior.alpha / std::pow(1.0 + depth + 1.0, prior.beta));
        
        double log_accept = ll_after - ll_before + prior_ratio;
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        if (std::log(unif(engine)) < log_accept) {
            node->is_leaf = true;
            node->feature = -1;
            node->threshold = 0.0;
            delete node->left; node->left = nullptr;
            delete node->right; node->right = nullptr;
        }
    }
    
    void update_tree(const mat& X, const vec& R, const vec& w_star, double sigma2) {
        std::vector<std::vector<int>> leaf_idx(node_counter + 100); 
        std::vector<int> all_idx(X.n_rows);
        std::iota(all_idx.begin(), all_idx.end(), 0);
        root->route(X, all_idx, leaf_idx);
        
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        if (unif(engine) < 0.5) grow_proposal(X, R, w_star, sigma2, leaf_idx);
        else prune_proposal(X, R, w_star, sigma2, leaf_idx);
        
        gibbs_sample_leaves(X, R, w_star, sigma2);
    }
    
    void gibbs_sample_leaves(const mat& X, const vec& R, const vec& w_star, double sigma2) {
        std::vector<Node*> leaves;
        root->get_leaves(leaves);
        
        std::vector<std::vector<int>> leaf_idx(node_counter + 100);
        std::vector<int> all_idx(X.n_rows);
        std::iota(all_idx.begin(), all_idx.end(), 0);
        root->route(X, all_idx, leaf_idx);
        
        std::normal_distribution<double> norm(0.0, 1.0);
        
        for (Node* l : leaves) {
            const std::vector<int>& idx = leaf_idx[l->node_id];
            if (idx.empty()) {
                l->mu = norm(engine) * prior.sigma_mu;
                continue;
            }
            
            double W = 0.0;
            double WR = 0.0;
            for (int i : idx) {
                W += w_star(i);
                WR += w_star(i) * R(i);
            }
            
            if (W <= 1e-8) {
                l->mu = norm(engine) * prior.sigma_mu;
            } else {
                double R_bar = WR / W;
                double V = 1.0 / (W / sigma2 + 1.0 / std::pow(prior.sigma_mu, 2));
                double mu_hat = V * W * R_bar / sigma2;
                l->mu = mu_hat + std::sqrt(V) * norm(engine);
            }
        }
    }
};

// ---------------------------------------------------------
// MCMC Runner
// ---------------------------------------------------------
// [[Rcpp::export]]
List fit_bart_cpp(mat X, vec y, mat W_mat, int n_trees, int n_mcmc, int burn_in, String family, Nullable<NumericMatrix> X_test = R_NilValue) {
    int n = X.n_rows;
    int n_draws = n_mcmc - burn_in;
    mat y_hat_train = zeros<mat>(n_draws, n);
    
    bool has_test = X_test.isNotNull();
    mat X_test_mat;
    mat y_hat_test;
    int n_test = 0;
    if (has_test) {
        X_test_mat = as<mat>(X_test);
        n_test = X_test_mat.n_rows;
        y_hat_test = zeros<mat>(n_draws, n_test);
    }
    
    bool is_probit = (family == "probit");
    bool is_hurdle = (family == "hurdle");
    
    Prior p; p.alpha = 0.95; p.beta = 2.0; p.k = 2.0; p.nu = 3.0;
    double sigma_hat = std::sqrt(var(y) * 0.9);
    p.lam = (sigma_hat * sigma_hat * R::qchisq(0.1, p.nu, 1, 0)) / p.nu;
    
    std::vector<Tree*> trees_probit;
    vec f_probit = zeros<vec>(n);
    vec Z = zeros<vec>(n);
    if (is_probit || is_hurdle) {
        Prior p_probit = p;
        p_probit.sigma_mu = 3.0 / (p_probit.k * std::sqrt(n_trees));
        for(int i=0; i<n_trees; i++) trees_probit.push_back(new Tree(p_probit));
    }
    
    std::vector<Tree*> trees_cont;
    vec y_cont, f_cont;
    std::vector<int> pos_idx;
    if (!is_probit) {
        if (is_hurdle) {
            for(int i=0; i<n; ++i) if(y(i) > 0) pos_idx.push_back(i);
            y_cont = zeros<vec>(pos_idx.size());
            for(size_t i=0; i<pos_idx.size(); ++i) y_cont(i) = std::log(y(pos_idx[i]));
            f_cont = zeros<vec>(pos_idx.size());
        } else {
            y_cont = y;
            f_cont = zeros<vec>(n);
        }
        
        Prior p_cont = p;
        double y_range = y_cont.max() - y_cont.min();
        p_cont.sigma_mu = y_range / (2.0 * p_cont.k * std::sqrt(n_trees));
        for(int i=0; i<n_trees; i++) trees_cont.push_back(new Tree(p_cont));
    }
    
    double sigma2 = is_probit ? 1.0 : var(y_cont) * 0.5;
    vec sigma2_draws = zeros<vec>(n_draws);
    
    for (int iter = 0; iter < n_mcmc; ++iter) {
        vec w_star = W_mat.row(iter).t();
        
        if (is_probit || is_hurdle) {
            for(int i=0; i<n; ++i) {
                if (y(i) > 0) Z(i) = r_trunc_norm(f_probit(i), 0.0, R_PosInf);
                else Z(i) = r_trunc_norm(f_probit(i), R_NegInf, 0.0);
            }
            for (int t = 0; t < n_trees; ++t) {
                vec f_t = trees_probit[t]->predict(X);
                f_probit -= f_t;
                vec R = Z - f_probit;
                trees_probit[t]->update_tree(X, R, w_star, 1.0);
                f_t = trees_probit[t]->predict(X);
                f_probit += f_t;
            }
        }
        
        if (!is_probit) {
            vec w_cont; mat X_cont;
            if (is_hurdle) {
                w_cont = zeros<vec>(pos_idx.size());
                X_cont = zeros<mat>(pos_idx.size(), X.n_cols);
                for(size_t i=0; i<pos_idx.size(); ++i) {
                    w_cont(i) = w_star(pos_idx[i]);
                    X_cont.row(i) = X.row(pos_idx[i]);
                }
            } else {
                w_cont = w_star;
                X_cont = X;
            }
            
            for (int t = 0; t < n_trees; ++t) {
                vec f_t = trees_cont[t]->predict(X_cont);
                f_cont -= f_t;
                vec R = y_cont - f_cont;
                trees_cont[t]->update_tree(X_cont, R, w_cont, sigma2);
                f_t = trees_cont[t]->predict(X_cont);
                f_cont += f_t;
            }
            
            double sum_sq = sum(w_cont % square(y_cont - f_cont));
            double df = sum(w_cont);
            if (df > 0) sigma2 = (p.nu * p.lam + sum_sq) / R::rchisq(df + p.nu);
        }
        
        if (iter >= burn_in) {
            int d = iter - burn_in;
            sigma2_draws(d) = sigma2;
            
            // 1. Train predictions
            if (is_probit) {
                for(int i=0; i<n; ++i) y_hat_train(d, i) = R::pnorm(f_probit(i), 0.0, 1.0, 1, 0);
            } else if (is_hurdle) {
                vec f_cont_full = zeros<vec>(n);
                for(int t=0; t<n_trees; ++t) f_cont_full += trees_cont[t]->predict(X);
                for(int i=0; i<n; ++i) {
                    double p_pos = R::pnorm(f_probit(i), 0.0, 1.0, 1, 0);
                    double E_amt = std::exp(f_cont_full(i) + sigma2 / 2.0);
                    y_hat_train(d, i) = p_pos * E_amt;
                }
            } else {
                y_hat_train.row(d) = f_cont.t();
            }
            
            // 2. Test predictions (G-computation)
            if (has_test) {
                if (is_probit) {
                    vec f_test = zeros<vec>(n_test);
                    for(int t=0; t<n_trees; ++t) f_test += trees_probit[t]->predict(X_test_mat);
                    for(int i=0; i<n_test; ++i) y_hat_test(d, i) = R::pnorm(f_test(i), 0.0, 1.0, 1, 0);
                } else if (is_hurdle) {
                    vec f_test_prob = zeros<vec>(n_test);
                    for(int t=0; t<n_trees; ++t) f_test_prob += trees_probit[t]->predict(X_test_mat);
                    
                    vec f_test_cont = zeros<vec>(n_test);
                    for(int t=0; t<n_trees; ++t) f_test_cont += trees_cont[t]->predict(X_test_mat);
                    
                    for(int i=0; i<n_test; ++i) {
                        double p_pos = R::pnorm(f_test_prob(i), 0.0, 1.0, 1, 0);
                        double E_amt = std::exp(f_test_cont(i) + sigma2 / 2.0);
                        y_hat_test(d, i) = p_pos * E_amt;
                    }
                } else {
                    vec f_test = zeros<vec>(n_test);
                    for(int t=0; t<n_trees; ++t) f_test += trees_cont[t]->predict(X_test_mat);
                    y_hat_test.row(d) = f_test.t();
                }
            }
        }
    }
    
    if (is_probit || is_hurdle) for(int i=0; i<n_trees; i++) delete trees_probit[i];
    if (!is_probit) for(int i=0; i<n_trees; i++) delete trees_cont[i];
    
    if (has_test) {
        return List::create(Named("y_hat_train") = y_hat_train, Named("y_hat_test") = y_hat_test, Named("sigma2_draws") = sigma2_draws);
    } else {
        return List::create(Named("y_hat_train") = y_hat_train, Named("y_hat_test") = R_NilValue, Named("sigma2_draws") = sigma2_draws);
    }
}
