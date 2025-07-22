import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# -----------------------
# Configuration
# -----------------------
sheet_names = ["scattering", "retardance", "orientation"]
sheet_names = ["orientation"]  # Debug single sheet

output_excel_path = (
    "/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/"
    "trace_summaries_18Jul2025__outliers_removed__non_centered_parameterization.xlsx"
)

view_distro = False
run_model = True

# -----------------------
# Begin Excel writing context
# -----------------------
with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:

    for sheet in sheet_names:
        # -------------------------
        # Import and preprocess data
        # -------------------------
        data = pd.read_excel(
            '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
            'lmm_test1_same_id_40um_donut_40um_outer.xlsx',
            sheet_name=sheet
        )
        data.columns = data.columns.str.strip().str.lower()
        data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

        # Convert to categorical indices
        data['condition_code'] = pd.Categorical(data['condition']).codes
        data['region_code'] = pd.Categorical(data['region']).codes
        data['subject_code'] = pd.Categorical(data['subject']).codes

        # -------------------------
        # Outlier Removal (IQR Method)
        # -------------------------
        if sheet in ["scattering", "retardance"]:
            Q1 = np.percentile(data["y"], 25)
            Q3 = np.percentile(data["y"], 75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            data = data[(data["y"] >= lower_bound) & (data["y"] <= upper_bound)]
            print(f"\nSheet: {sheet} - {len(data)} samples after outlier removal.")

        # -------------------------
        # Remove nonpositive, NaN, Inf, and log-transform
        # -------------------------
        data = data[data["y"] > 0]
        data = data[np.isfinite(data["y"])]
        data = data.dropna(subset=["y"])
        assert len(data) > 10, "Too few observations left after cleaning!"
        assert (data["y"] > 0).all(), "All values must be > 0 for log transform"
        # Log transform
        data["log_y"] = np.log(data["y"])

        # -------------------------
        # Set up modeling variables
        # -------------------------
        y_obs = data["log_y"].values
        condition = data["condition_code"].values
        subject_idx = data["subject_code"].values
        tissue_idx = data["region_code"].values

        n_subjects = len(np.unique(subject_idx))
        n_tissues = len(np.unique(tissue_idx))

        # -------------------------
        # Estimate priors
        # -------------------------
        ctrl_vals = data[data["condition_code"] == 0]["log_y"]
        exp_vals = data[data["condition_code"] == 1]["log_y"]
        mu_beta = exp_vals.mean() - ctrl_vals.mean()
        sigma_beta = np.sqrt(
            ctrl_vals.var(ddof=1) / len(ctrl_vals) +
            exp_vals.var(ddof=1) / len(exp_vals)
        )
        mu_inter = y_obs.mean()
        sigma_inter = y_obs.std(ddof=1)

        subject_std = data.groupby("subject_code")["log_y"].mean().std(ddof=1)
        tissue_std = data.groupby("region_code")["log_y"].mean().std(ddof=1)
        group_means = data.groupby("condition_code")["log_y"].transform("mean")
        residual_std = np.sqrt(np.mean((y_obs - group_means) ** 2))

        # -------------------------
        # View Distributions
        # -------------------------
        if view_distro:
            sns.histplot(y_obs, kde=True)
            plt.title(f'Histogram of {sheet} (Outliers Removed)')
            plt.show()

        # -------------------------
        # Run PyMC Model
        # -------------------------
        if run_model:
            with pm.Model() as model:
                mu_intercept = pm.Normal("mu_intercept", mu=mu_inter, sigma=sigma_inter)
                beta_condition = pm.Normal("beta_condition", mu=mu_beta, sigma=sigma_beta)

                sigma_subject = pm.HalfNormal("sigma_subject", sigma=subject_std)
                sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=tissue_std)
                sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

                if sheet == "scattering":
                    z_subject = pm.Normal("z_subject", mu=0, sigma=1, shape=n_subjects)
                    z_tissue = pm.Normal("z_tissue", mu=0, sigma=1, shape=n_tissues)
                    subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)
                    tissue_effect = pm.Deterministic("tissue_effect", z_tissue * sigma_tissue)
                else:
                    subject_effect = pm.Normal("subject_effect", mu=0, sigma=sigma_subject, shape=n_subjects)
                    tissue_effect = pm.Normal("tissue_effect", mu=0, sigma=sigma_tissue, shape=n_tissues)

                mu = pm.Deterministic("mu", mu_intercept + beta_condition * condition +
                                      subject_effect[subject_idx] + tissue_effect[tissue_idx])

                y = pm.Normal("y", mu=mu, sigma=sigma_residual, observed=y_obs)

                trace = pm.sample(draws=4000, tune=4000, target_accept=0.99,
                                  chains=4, cores=2, random_seed=42, progressbar=False)

            # -------------------------
            # Save Trace Summary
            # -------------------------
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
