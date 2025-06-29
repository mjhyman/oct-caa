import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# -----------------------
# Import spreadsheet of data
# -----------------------
# spreadsheet sheet names for each optical property
sheet_names = ["scattering", "retardance", "orientation"]

# Optional: store summaries for each response
summary_dict = {}

# Excel writer
output_excel_path = "/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/trace_summaries_29Jun2025.xlsx"
writer = pd.ExcelWriter(output_excel_path, engine="openpyxl")

# Booleans for viewing distributions or running MCMC
view_distro = True
run_model = True

for sheet in sheet_names:
    # Read in Excel file
    data = pd.read_excel('/autofs/cluster/octdata3/users/mjhyman/'
                         'oct_caa_analyses/optical_properties/lmm_test1_same_id_100um_donut_100um_outer.xlsx',
                         sheet_name=sheet)

    # Extract columns (subjectID, region, optical property)
    data.columns = data.columns.str.strip().str.lower()

    # Encode categorical IDs as integer indices
    data['condition'] = pd.Categorical(data['groups']).codes
    data['region'] = pd.Categorical(data['region']).codes
    data['subject'] = pd.Categorical(data['subid']).codes
    data['y'] = data['opticalproperty']

    # Extract necessary variables
    subject_idx = data["subject"].values
    tissue_idx = data["region"].values
    condition = data["condition"].values.astype(int)
    y_obs = data["y"]

    n_subjects = len(np.unique(subject_idx))
    n_tissues = len(np.unique(tissue_idx))

    # Automatically estimate prior for beta_condition
    ctrl_vals = data[data['condition'] == 0]['y']
    exp_vals = data[data['condition'] == 1]['y']
    mu_beta = exp_vals.mean() - ctrl_vals.mean()
    sigma_beta = np.sqrt(
        ctrl_vals.var(ddof=1) / len(ctrl_vals) +
        exp_vals.var(ddof=1) / len(exp_vals)
    )

    # Automatically estimate prior for mu_intercept
    mu_inter = y_obs.mean()
    sigma_inter = y_obs.std(ddof=1)

    if view_distro:
        # -------------------------
        # Quality control - check distribution of data
        # -------------------------
        sns.histplot(y_obs,kde=True)
        plt.title(f'Histogram of {sheet}')
        plt.show()
        print(f"\nSummary statistics for {sheet}:\n")
        print("Mean:", np.mean(y_obs))
        print("Std Dev:", np.std(y_obs, ddof=1))
        print("Min:", np.min(y_obs))
        print("Max:", np.max(y_obs))
        print("Skew:", stats.skew(y_obs))
        print("Kurtosis:", stats.kurtosis(y_obs))

    if run_model:
        # -------------------------
        # Bayesian Model with PyMC
        # -------------------------
        with pm.Model() as model:
            # Distributions of mean optical property
            mu_intercept = pm.Normal("mu_intercept", mu=mu_inter, sigma=sigma_inter)

            # Distribution of difference between means of groups (experimental - control)
            beta_condition = pm.Normal("beta_condition", mu=mu_beta, sigma=sigma_beta)

            # Standard Deviations of subject, tissue region, and residual (unexplained)
            sigma_subject = pm.HalfNormal("sigma_subject", sigma=1.0)
            sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=1.0)
            sigma_residual = pm.HalfNormal("sigma_residual", sigma=1.0)

            # Distributions of subjects and tissues
            subject_effect = pm.Normal("subject_effect", mu=0, sigma=sigma_subject, shape=n_subjects)
            tissue_effect = pm.Normal("tissue_effect", mu=0, sigma=sigma_tissue, shape=n_tissues)

            # Formula for modeling the data
            mu = mu_intercept + beta_condition * condition + subject_effect[subject_idx] + tissue_effect[tissue_idx]

            # LogNormal if strictly positive, else fallback to Normal
            if (y_obs > 0).all() and stats.skew(y_obs) > 1:
                y = pm.LogNormal("y", mu=mu, sigma=sigma_residual, observed=y_obs)
            else:
                y = pm.Normal("y", mu=mu, sigma=sigma_residual, observed=y_obs)

            # Sampling
            trace = pm.sample(1000, tune=1000, target_accept=0.9, cores=2, random_seed=42)
        # -------------------------
        # Posterior Summaries
        # -------------------------
        # Create summary line
        summary_df = az.summary(trace, var_names=["beta_condition", "mu_intercept",
                                                  "sigma_subject", "sigma_tissue", "sigma_residual"])
        # Save to Excel
        summary_df.to_excel(writer, sheet_name=sheet)
        # Plot posterior
        az.plot_posterior(trace, var_names=["beta_condition"])
        plt.title(f'Posterior of {sheet}')
        plt.show()

# Final save to disk
writer.close()