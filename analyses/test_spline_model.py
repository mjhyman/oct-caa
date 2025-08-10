import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from patsy import dmatrix

# Synthetic data creation
np.random.seed(42)
n_subjects = 5
n_points_per_subject = 20
subjects = np.repeat(np.arange(n_subjects), n_points_per_subject)
distance = np.linspace(0, 10, n_points_per_subject).reshape(1, -1).repeat(n_subjects, axis=0).flatten()
condition = np.random.binomial(1, 0.5, size=n_subjects * n_points_per_subject)
region = np.random.choice(['front', 'occip'], size=n_subjects * n_points_per_subject)
# Create true spline bases
spline_basis = dmatrix("bs(distance, df=7, degree=3, include_intercept=False) - 1",
                      {"distance": distance}, return_type='dataframe')
X_base = spline_basis.values
n_splines = X_base.shape[1]

# True betas for conditions (just for simulation)
beta_cond0_true = np.sin(np.linspace(0, 3, n_splines))
beta_cond1_true = np.cos(np.linspace(0, 3, n_splines))
# Simulate y (on log scale)
log_y = (X_base * (1 - condition)[:, None]) @ beta_cond0_true + \
        (X_base * condition[:, None]) @ beta_cond1_true
log_y += np.random.normal(0, 0.1, size=log_y.shape)  # noise
y = np.exp(log_y)

# Create DataFrame
data = pd.DataFrame({
    "y": y,
    "distance": distance,
    "condition": condition,
    "subject": subjects,
    "region": region
})

# Encode categorical
condition_code = pd.Categorical(data["condition"]).codes
subject_idx = pd.Categorical(data["subject"]).codes

with pm.Model() as model:
    # Priors
    mu_intercept = pm.Normal("mu_intercept", mu=np.log(np.median(y)), sigma=1)
    sigma_subject = pm.HalfNormal("sigma_subject", sigma=1)
    sigma_residual = pm.HalfNormal("sigma_residual", sigma=1)

    # Subject effects
    subject_effect = pm.Normal("subject_effect", mu=0, sigma=sigma_subject, shape=len(np.unique(subject_idx)))

    # Spline coefficients for condition 0 and 1
    beta_spline = pm.Normal("beta_spline", mu=0, sigma=1, shape=2 * n_splines)
    beta_diff = pm.Deterministic("beta_diff", beta_spline[n_splines:] - beta_spline[:n_splines])

    # Design matrix for splines duplicated per condition
    X_cond0 = X_base * (1 - condition_code)[:, None]
    X_cond1 = X_base * condition_code[:, None]
    X_spline = np.hstack([X_cond0, X_cond1])

    # Linear model
    mu = mu_intercept + pm.math.dot(X_spline, beta_spline) + subject_effect[subject_idx]

    # Likelihood
    y_obs = pm.LogNormal("y_obs", mu=mu, sigma=sigma_residual, observed=data["y"].values)

    # Sample
    trace = pm.sample(draws=1000, tune=1000, chains=2, target_accept=0.9, random_seed=42)

# Extract posterior beta_diff samples and posterior mean
distance_grid = np.linspace(data["distance"].min(), data["distance"].max(), 100)
spline_basis_grid = dmatrix(
    "bs(distance_grid, df=7, degree=3, include_intercept=False) - 1",
    {"distance_grid": distance_grid}, return_type='dataframe')
X_grid = spline_basis_grid.values

beta_diff_samples = trace.posterior["beta_diff"].stack(sample=("chain", "draw")).values.T  # shape(samples, n_splines)
diff_posterior = beta_diff_samples @ X_grid.T  # shape(samples, 100)

diff_mean = diff_posterior.mean(axis=0)
diff_hpd = az.hdi(diff_posterior, hdi_prob=0.95)
significant = (diff_hpd[:, 0] > 0) | (diff_hpd[:, 1] < 0)

# Plot difference spline with shaded credible interval + significant distances
plt.figure(figsize=(8, 5))
plt.plot(distance_grid, diff_mean, label="Difference (Cond1 - Cond0)", color='black')
plt.fill_between(distance_grid, diff_hpd[:, 0], diff_hpd[:, 1], color='gray', alpha=0.3, label="95% Credible Interval")

# Shade significant intervals
def highlight_intervals(x, sig_mask, color='yellow', alpha=0.2):
    ax = plt.gca()
    in_sig = False
    start = None
    for i, val in enumerate(sig_mask):
        if val and not in_sig:
            in_sig = True
            start = x[i]
        elif not val and in_sig:
            ax.axvspan(start, x[i-1], color=color, alpha=alpha)
            in_sig = False
    if in_sig:
        ax.axvspan(start, x[-1], color=color, alpha=alpha)

highlight_intervals(distance_grid, significant)

plt.xlabel("Distance")
plt.ylabel("Difference in log-scale spline fit")
plt.title("Posterior Difference Between Conditions (Synthetic Data)")
plt.legend()
plt.tight_layout()
plt.show()