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
sheet_names = ["scattering", "retardance", "orientation"]
# sheet_names = ["scattering"]  # For debugging one sheet only

# Output file path
output_excel_path = (
    "/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/"
    "trace_summaries_16Jul2025__outliers_removed__non_centered_parameterization.xlsx"
)

# Booleans for viewing distributions or running MCMC
view_distro = False
run_model = True

# -----------------------
# Begin Excel writing context
# -----------------------
with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:

    for sheet in sheet_names:
        # Load data
        data = pd.read_excel(
            '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
            'lmm_test1_same_id_40um_donut_40um_outer.xlsx',
            sheet_name=sheet
        )
        data.columns = data.columns.str.strip().str.lower()

        # -------------------------
        # Outlier Removal (IQR Method)
        # -------------------------
        y_obs = data["opticalproperty"]
        Q1 = np.percentile(y_obs, 25)
        Q3 = np.percentile(y_obs, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        data = data[(y_obs >= lower_bound) & (y_obs <= upper_bound)]
        y_obs = data["opticalproperty"]  # Refresh after filtering

        print(f"\nSheet: {sheet} - {len(y_obs)} samples after outlier removal.")

        # Encode IDs
        data['condition'] = pd.Categorical(data['groups']).codes
        data['region'] = pd.Categorical(data['region']).codes
        data['subject'] = pd.Categorical(data['subid']).codes
        data['y'] = data['opticalproperty']

        subject_idx = data["subject"].values
        tissue_idx = data["region"].values
        condition = data["condition"].values.astype(int)
        y_obs = data["y"]

        n_subjects = len(np.unique(subject_idx))
        n_tissues = len(np.unique(tissue_idx))

        # Estimate priors
        ctrl_vals = data[data['condition'] == 0]['y']
        exp_vals = data[data['condition'] == 1]['y']
        mu_beta = exp_vals.mean() - ctrl_vals.mean()
        sigma_beta = np.sqrt(
            ctrl_vals.var(ddof=1) / len(ctrl_vals) +
            exp_vals.var(ddof=1) / len(exp_vals)
        )
        mu_inter = y_obs.mean()
        sigma_inter = y_obs.std(ddof=1)

        # Empirical priors for sigmas
        subject_std = data.groupby("subject")["y"].mean().std(ddof=1)
        tissue_std = data.groupby("region")["y"].mean().std(ddof=1)
        residual_std = np.sqrt(np.mean((y_obs - data.groupby("condition")["y"].transform("mean"))**2))

        # View distributions
        if view_distro:
            sns.histplot(y_obs, kde=True)
            plt.title(f'Histogram of {sheet} (Outliers Removed)')
            plt.show()

        # Run bayesian MCMC model
        with pm.Model() as model:
            mu_intercept = pm.Normal("mu_intercept", mu=mu_inter, sigma=sigma_inter)
            beta_condition = pm.Normal("beta_condition", mu=mu_beta, sigma=sigma_beta)

            sigma_subject = pm.HalfNormal("sigma_subject", sigma=subject_std)
            sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=tissue_std)
            sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

            if sheet == "scattering":
                z_subject = pm.Normal("z_subject", mu=0, sigma=1.0, shape=n_subjects)
                z_tissue = pm.Normal("z_tissue", mu=0, sigma=1.0, shape=n_tissues)
                subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)
                tissue_effect = pm.Deterministic("tissue_effect", z_tissue * sigma_tissue)
            else:
                subject_effect = pm.Normal("subject_effect", mu=0, sigma=sigma_subject, shape=n_subjects)
                tissue_effect = pm.Normal("tissue_effect", mu=0, sigma=sigma_tissue, shape=n_tissues)

            mu = mu_intercept + beta_condition * condition + subject_effect[subject_idx] + tissue_effect[tissue_idx]
            y = pm.LogNormal("y", mu=mu, sigma=sigma_residual, observed=y_obs)

            trace = pm.sample(draws=4000, tune=4000, target_accept=0.99,
                                chains=2, cores=1, random_seed=42, progressbar=False)

        summary_df = az.summary(trace, var_names=[
            "beta_condition", "mu_intercept", "sigma_subject", "sigma_tissue", "sigma_residual"
        ])
        summary_df.to_excel(writer, sheet_name=sheet)

        print(f"\nSummary statistics for sheet '{sheet}':")
        print(summary_df.to_string())

        if view_distro:
            az.plot_posterior(trace, var_names=["beta_condition"])
            plt.title(f'Posterior of {sheet} (Outliers Removed)')
            plt.show()