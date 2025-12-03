import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# -----------------------
# Configuration
# -----------------------
sheet_names = ["scattering", "retardance", "orientation"]

output_excel_path = (
    "/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/statistics/"
    "trace_summaries_20Jul2025__outliers_removed__non_centered_parameterization__regional.xlsx"
)

posterior_dir = "/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/statistics"
os.makedirs(posterior_dir, exist_ok=True)

view_distro = False
run_model = True

def run_pymc_model(data, sheet, region_str, use_tissue_effect=True):
    data.columns = data.columns.str.strip().str.lower()
    data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

    # Convert to categorical indices
    data['condition_code'] = pd.Categorical(data['condition']).codes
    data['region_code'] = pd.Categorical(data['region']).codes
    data['subject_code'] = pd.Categorical(data['subject']).codes

    y_obs = data["y"].values
    condition = data["condition_code"].values
    subject_idx = data["subject_code"].values
    tissue_idx = data["region_code"].values

    n_subjects = len(np.unique(subject_idx))
    n_tissues = len(np.unique(tissue_idx))

    # Priors
    ctrl_vals = data[data["condition_code"] == 0]["y"]
    exp_vals = data[data["condition_code"] == 1]["y"]
    mu_beta = np.log(exp_vals.mean() / ctrl_vals.mean())  # ratio effect
    sigma_beta = np.sqrt((ctrl_vals.std(ddof=1) / ctrl_vals.mean()) ** 2 +
                         (exp_vals.std(ddof=1) / exp_vals.mean()) ** 2)
    # mean and std dev of intercept
    mu_inter = np.log(y_obs.mean())
    sigma_inter = y_obs.std(ddof=1) / y_obs.mean()

    #
    subject_std = data.groupby("subject_code")["y"].mean().std(ddof=1)
    residuals = y_obs - data.groupby("condition_code")["y"].transform("mean")
    residual_std = residuals.std(ddof=1)

    with pm.Model() as model:
        mu_intercept = pm.Normal("mu_intercept", mu=mu_inter, sigma=sigma_inter)
        beta_condition = pm.Normal("beta_condition", mu=mu_beta, sigma=sigma_beta)

        sigma_subject = pm.HalfNormal("sigma_subject", sigma=subject_std)
        sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

        if sheet == "scattering":
            z_subject = pm.Normal("z_subject", mu=0, sigma=1, shape=n_subjects)
            subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)
        else:
            subject_effect = pm.Normal("subject_effect", mu=0, sigma=sigma_subject, shape=n_subjects)

        if use_tissue_effect:
            tissue_std = data.groupby("region_code")["y"].mean().std(ddof=1)
            sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=tissue_std)
            if sheet == "scattering":
                z_tissue = pm.Normal("z_tissue", mu=0, sigma=1, shape=n_tissues)
                tissue_effect = pm.Deterministic("tissue_effect", z_tissue * sigma_tissue)
            else:
                tissue_effect = pm.Normal("tissue_effect", mu=0, sigma=sigma_tissue, shape=n_tissues)

            mu = pm.Deterministic("mu", mu_intercept + beta_condition * condition +
                                  subject_effect[subject_idx] + tissue_effect[tissue_idx])
        else:
            mu = pm.Deterministic("mu", mu_intercept + beta_condition * condition +
                                  subject_effect[subject_idx])

        y = pm.LogNormal("y", mu=mu, sigma=sigma_residual, observed=y_obs)

        trace = pm.sample(draws=4000, tune=2000, target_accept=0.97,
                          chains=4, cores=2, random_seed=42, progressbar=False)

    summary_vars = ["beta_condition", "mu_intercept", "sigma_subject", "sigma_residual"]
    if use_tissue_effect:
        summary_vars.append("sigma_tissue")

    summary_df = az.summary(trace, var_names=summary_vars)
    print(f"\nSummary statistics for sheet '{sheet}' and region '{region_str}':")
    print(summary_df.to_string())

    return trace, summary_df


with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
    for sheet in sheet_names:
        full_data = pd.read_excel(
            '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
            'lmm_test1_same_id_40um_donut_40um_outer.xlsx',
            sheet_name=sheet
        )

        full_data.columns = full_data.columns.str.strip().str.lower()
        full_data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

        if sheet in ["scattering", "retardance"]:
            Q1 = np.percentile(full_data["y"], 25)
            Q3 = np.percentile(full_data["y"], 75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            full_data = full_data[(full_data["y"] >= lower_bound) & (full_data["y"] <= upper_bound)]
            print(f"\nSheet: {sheet} - {len(full_data)} samples after outlier removal.")

        full_data = full_data[full_data["y"] > 0]
        full_data = full_data[np.isfinite(full_data["y"])]
        full_data = full_data.dropna(subset=["y"])
        assert len(full_data) > 10, "Too few observations left after cleaning!"
        assert (full_data["y"] > 0).all(), "All values must be > 0 for LogNormal"

        posterior_samples_dict = {}

        # Frontal
        print(f"Running model for sheet '{sheet}' with frontal region only.")
        frontal_data = full_data[full_data["region"].str.lower() == "front"].copy()
        if len(frontal_data) > 10:
            trace_frontal, summary_frontal = run_pymc_model(frontal_data, sheet, 'frontal', False)
            summary_frontal.to_excel(writer, sheet_name=f"{sheet}_frontal")
            posterior_samples_frontal = trace_frontal.posterior["beta_condition"].values.flatten()
            posterior_samples_dict["Frontal"] = posterior_samples_frontal
            pd.DataFrame(posterior_samples_frontal, columns=["beta_condition"]).to_csv(
                os.path.join(posterior_dir, f"{sheet}_frontal_posterior.csv"), index=False)
            az.plot_posterior(trace_frontal, var_names=["beta_condition"])
            plt.title(f"Posterior of beta_condition - {sheet} (Frontal only)")
            plt.savefig(os.path.join(posterior_dir, f"{sheet}_frontal_posterior.png"))
            plt.close()

        # Occipital
        print(f"Running model for sheet '{sheet}' with occipital region only.")
        occipital_data = full_data[full_data["region"].str.lower() == "occip"].copy()
        if len(occipital_data) > 10:
            trace_occipital, summary_occipital = run_pymc_model(occipital_data, sheet, 'occipital', False)
            summary_occipital.to_excel(writer, sheet_name=f"{sheet}_occipital")
            posterior_samples_occipital = trace_occipital.posterior["beta_condition"].values.flatten()
            posterior_samples_dict["Occipital"] = posterior_samples_occipital
            pd.DataFrame(posterior_samples_occipital, columns=["beta_condition"]).to_csv(
                os.path.join(posterior_dir, f"{sheet}_occipital_posterior.csv"), index=False)
            az.plot_posterior(trace_occipital, var_names=["beta_condition"])
            plt.title(f"Posterior of beta_condition - {sheet} (Occipital only)")
            plt.savefig(os.path.join(posterior_dir, f"{sheet}_occipital_posterior.png"))
            plt.close()

        # All regions combined
        print(f"\nRunning model for sheet '{sheet}' with all regions combined.")
        trace_all, summary_all = run_pymc_model(full_data.copy(), sheet, 'combined', True)
        summary_all.to_excel(writer, sheet_name=f"{sheet}_all")
        posterior_samples_all = trace_all.posterior["beta_condition"].values.flatten()
        posterior_samples_dict["All regions"] = posterior_samples_all
        pd.DataFrame(posterior_samples_all, columns=["beta_condition"]).to_csv(
            os.path.join(posterior_dir, f"{sheet}_all_posterior.csv"), index=False)
        az.plot_posterior(trace_all, var_names=["beta_condition"])
        plt.title(f"Posterior of beta_condition - {sheet} (All regions)")
        plt.savefig(os.path.join(posterior_dir, f"{sheet}_all_posterior.png"))
        plt.close()

        # Overlay plot
        plt.figure(figsize=(10, 6))
        for label, samples in posterior_samples_dict.items():
            sns.kdeplot(samples, label=label, fill=True, alpha=0.4)
        plt.title(f"Overlay Posterior Distributions of beta_condition - {sheet}")
        plt.xlabel("beta_condition")
        plt.ylabel("Density")
        plt.legend()
        plt.tight_layout()
        overlay_path = os.path.join(posterior_dir, f"{sheet}_overlay_posterior.png")
        plt.savefig(overlay_path)
        plt.close()
        print(f"Overlay posterior plot saved to {overlay_path}")
