import pandas as pd

# read data
data = pd.read_csv('data/population.csv');

# Unique cluster identifier (since psuid is nested within strata)
data['cluster_id'] = data['strata'].astype(str) + '_' + data['psuid'].astype(str)

# number of observations, strata and clusters
n_obs = len(data)
n_strata = data['strata'].nunique()
n_clusters = data['cluster_id'].nunique()

# Print summary data
print('Total number of observations (N): ', n_obs);
print('Total number of strata: ', n_strata);
print('Total number of unique clusters (PSUs): ', n_clusters);

# --- Estimate R^2_alpha (Sharon Lohr's ICC Estimator) ---
# S^2 is the overall sample variance of y
s2 = data['y'].var()
s = data['y'].std()

# MSW is the mean square within-cluster variance
cluster_means = data.groupby('cluster_id')['y'].transform('mean');
ssw = ((data['y'] - cluster_means) ** 2).sum();
msw = ssw / (n_obs - n_clusters);

# Lohr's adjusted R^2 (R^2_alpha)
r2_alpha = 1 - (msw / s2)
print('\n--- Lohr Sample Size & Design Effect Parameters ---')
print('S^2 (Total Variance): ', s2)
print('MSW (Within-cluster Mean Square): ', msw)
print('Lohr R^2_alpha (adjusted R^2 for ICC): ', r2_alpha)

# Average cluster size (m_i)
m_i = n_obs / n_clusters
print('Average cluster size (m_i): ', m_i)

# Design Effect (DEF)
deff = 1 + r2_alpha * (m_i - 1)  # Equivalent to 1 - r2_alpha * (1 - m_i)
print('Estimated Design Effect (DEF): ', deff)

# Sample size estimation assuming SRS with e=0.05
n_srs = (1.96 * s / 0.5)**2 
print('Sample size assuming SRS: ', n_srs)
# Note: FPC is ignored as this is a complex sample design (Lohr book)

# Sample size estimation with clustering assuming conservative DEFF_w (~ 15%)
n = n_srs * deff * 1.15
print('Final sample size: ', round(n,-2)) 

