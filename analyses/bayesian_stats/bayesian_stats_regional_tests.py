import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
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

    # Extract the observation data
    y_obs = data["log_y"].values
    condition = data["condition_code"].values
    subject_idx = data["subject_code"].values
    tissue_idx = data["region_code"].values

    # Calculate number of subjects and tissues
    n_subjects = len(np.unique(subject_idx))
    n_tissues = len(np.unique(tissue_idx))

    # Priors
    ctrl_vals = data[data["condition_code"] == 0]["log_y"]
    exp_vals = data[data["condition_code"] == 1]["log_y"]
    mu_beta = exp_vals.mean() - ctrl_vals.mean()
    sigma_beta = np.sqrt(
        ctrl_vals.var(ddof=1) / len(ctrl_vals) +
        exp_vals.var(ddof=1) / len(exp_vals)
    )
    mu_inter = y_obs.mean()
    sigma_inter = y_obs.std(ddof=1)

    # Measure standard deviation from samples
    subject_std = data.groupby("subject_code")["log_y"].mean().std(ddof=1)
    residual_std = np.sqrt(np.mean((y_obs - data.groupby("condition_code")["log_y"].transform("mean")) ** 2))

    with pm.Model() as model:
        # Mu and beta
        mu_intercept = pm.Normal("mu_intercept", mu=mu_inter, sigma=sigma_inter)
        beta_condition = pm.Normal("beta_condition", mu=mu_beta, sigma=sigma_beta)
        # Sigma of subject and residual
        sigma_subject = pm.HalfNormal("sigma_subject", sigma=subject_std)
        sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

        # Set the subject effect deterministically for scattering to avoid errors
        if sheet == "scattering":
            z_subject = pm.Normal("z_subject", mu=0, sigma=1, shape=n_subjects)
            subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)
        else:
            subject_effect = pm.Normal("subject_effect", mu=0, sigma=sigma_subject, shape=n_subjects)

        # If combining frontal and occipital, then incorporate the tissue effect
        if use_tissue_effect:
            tissue_std = data.groupby("region_code")["log_y"].mean().std(ddof=1)
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

        y = pm.Normal("y", mu=mu, sigma=sigma_residual, observed=y_obs)

        trace = pm.sample(draws=4000, tune=2000, target_accept=0.97,
                          chains=4, cores=2, random_seed=42, progressbar=False)

    # Save summary to Excel sheet
    summary_vars = ["beta_condition", "mu_intercept", "sigma_subject", "sigma_residual"]
    if use_tissue_effect:
        summary_vars.append("sigma_tissue")

    summary_df = az.summary(trace, var_names=summary_vars)
    print(f"\nSummary statistics for sheet '{sheet}' and region '{region_str}':")
    print(summary_df.to_string())

    return trace, summary_df


with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
    for sheet in sheet_names:
        # Read full data
        full_data = pd.read_excel(
            '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties/'
            'lmm_test1_same_id_40um_donut_40um_outer.xlsx',
            sheet_name=sheet
        )

        full_data.columns = full_data.columns.str.strip().str.lower()
        full_data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

        # Outlier removal (only for scattering and retardance)
        if sheet in ["scattering", "retardance"]:
            Q1 = np.percentile(full_data["y"], 25)
            Q3 = np.percentile(full_data["y"], 75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            full_data = full_data[(full_data["y"] >= lower_bound) & (full_data["y"] <= upper_bound)]
            print(f"\nSheet: {sheet} - {len(full_data)} samples after outlier removal.")

        # Remove nonpositive, NaN, Inf, log transform
        full_data = full_data[full_data["y"] > 0]
        full_data = full_data[np.isfinite(full_data["y"])]
        full_data = full_data.dropna(subset=["y"])
        assert len(full_data) > 10, "Too few observations left after cleaning!"
        assert (full_data["y"] > 0).all(), "All values must be > 0 for log transform"
        full_data["log_y"] = np.log(full_data["y"])

        # Prepare a dictionary to hold posterior samples for overlay plot
        posterior_samples_dict = {}

        # 1) Frontal region only
        print(f"Running model for sheet '{sheet}' with frontal region only.")
        frontal_data = full_data[full_data["region"].str.lower() == "front"].copy()
        if len(frontal_data) > 10:
            frontal_data["condition_code"] = pd.Categorical(frontal_data["condition"]).codes
            frontal_data["region_code"] = pd.Categorical(frontal_data["region"]).codes
            frontal_data["subject_code"] = pd.Categorical(frontal_data["subject"]).codes
            frontal_data["log_y"] = np.log(frontal_data["y"])

            # Call PY MCMC model
            trace_frontal, summary_frontal = run_pymc_model(frontal_data, sheet,'frontal',False)

            # Summarize the data
            summary_frontal.to_excel(writer, sheet_name=f"{sheet}_frontal")
            posterior_samples_frontal = trace_frontal.posterior["beta_condition"].values.flatten()
            posterior_samples_dict["Frontal"] = posterior_samples_frontal
            # Save posterior samples
            pd.DataFrame(posterior_samples_frontal, columns=["beta_condition"]).to_csv(
                os.path.join(posterior_dir, f"{sheet}_frontal_posterior.csv"), index=False
            )
            # Plot posterior and save
            az.plot_posterior(trace_frontal, var_names=["beta_condition"])
            plt.title(f"Posterior of beta_condition - {sheet} (Frontal only)")
            plt.savefig(os.path.join(posterior_dir, f"{sheet}_frontal_posterior.png"))
            plt.close()
        else:
            print(f"Not enough frontal data points to run model for sheet '{sheet}'.")

        # 2) Occipital region only
        print(f"Running model for sheet '{sheet}' with occipital region only.")
        occipital_data = full_data[full_data["region"].str.lower() == "occip"].copy()
        if len(occipital_data) > 10:
            occipital_data["condition_code"] = pd.Categorical(occipital_data["condition"]).codes
            occipital_data["region_code"] = pd.Categorical(occipital_data["region"]).codes
            occipital_data["subject_code"] = pd.Categorical(occipital_data["subject"]).codes
            occipital_data["log_y"] = np.log(occipital_data["y"])

            # Call MCMC Model
            trace_occipital, summary_occipital = run_pymc_model(occipital_data, sheet,'occipital',False)

            # Summarize the results
            summary_occipital.to_excel(writer, sheet_name=f"{sheet}_occipital")
            posterior_samples_occipital = trace_occipital.posterior["beta_condition"].values.flatten()
            posterior_samples_dict["Occipital"] = posterior_samples_occipital
            # Save posterior samples
            pd.DataFrame(posterior_samples_occipital, columns=["beta_condition"]).to_csv(
                os.path.join(posterior_dir, f"{sheet}_occipital_posterior.csv"), index=False
            )
            # Plot posterior and save
            az.plot_posterior(trace_occipital, var_names=["beta_condition"])
            plt.title(f"Posterior of beta_condition - {sheet} (Occipital only)")
            plt.savefig(os.path.join(posterior_dir, f"{sheet}_occipital_posterior.png"))
            plt.close()
        else:
            print(f"Not enough occipital data points to run model for sheet '{sheet}'.")

        # 3) All regions combined
        print(f"\nRunning model for sheet '{sheet}' with all regions combined.")

        # Call the MCMC model
        trace_all, summary_all = run_pymc_model(full_data.copy(), sheet, 'combined',True)

        # Summarize the results
        summary_all.to_excel(writer, sheet_name=f"{sheet}_all")
        posterior_samples_all = trace_all.posterior["beta_condition"].values.flatten()
        posterior_samples_dict["All regions"] = posterior_samples_all
        # Save posterior samples
        pd.DataFrame(posterior_samples_all, columns=["beta_condition"]).to_csv(
            os.path.join(posterior_dir, f"{sheet}_all_posterior.csv"), index=False
        )
        # Plot posterior and save
        az.plot_posterior(trace_all, var_names=["beta_condition"])
        plt.title(f"Posterior of beta_condition - {sheet} (All regions)")
        plt.savefig(os.path.join(posterior_dir, f"{sheet}_all_posterior.png"))
        plt.close()

        ### Overlay plot of the three posteriors
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