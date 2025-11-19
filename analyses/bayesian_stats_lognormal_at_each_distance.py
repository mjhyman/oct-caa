import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import argparse as ap
from datetime import datetime

# -----------------------
# Configuration
# -----------------------

# ---- Import the distance (microns) from the console for the batch script
try:
    parser = ap.ArgumentParser(description="Run OCT Bayesian analysis for specified distance")
    parser.add_argument('--distance', type=float, required=True, help='Distance in microns to filter data on')
    args = parser.parse_args()
    # Access distance via args.distance
    distance = int(args.distance)
    print(f"Running model for distance = {distance} microns")
except:
    distance = 40
    print(f"Running model with default distance (40 microns) as debug")

sheet_names = ["scattering", "retardance"]

# Create output filepath
date = datetime.now().strftime("%Y_%m_%d")
output_excel_path = (
    "/projectnb/npbssmic/ns/CAA/beta_stats/posteriors/"
    f"bayes_posterior_distance_{distance}_{date}.xlsx")

posterior_dir = "/projectnb/npbssmic/ns/CAA/beta_stats/posteriors/"
os.makedirs(posterior_dir, exist_ok=True)

def run_pymc_model(data: pd.DataFrame, sheet: str, region_str: str, use_tissue_effect=True, distance=None):
    data.columns = data.columns.str.strip().str.lower()
    data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

    # Filter by distance if specified
    if distance is not None:
        if "distance" not in data.columns:
            raise ValueError("Column 'distance' not found in data")
        data = data[data["distance"] == distance]
        if len(data) == 0:
            raise ValueError(f"No data found for distance={distance}")

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

        # Standard deviation of subject and residual
        sigma_subject = pm.HalfNormal("sigma_subject", sigma=subject_std)
        sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

        # Subject-level effect
        z_subject = pm.Normal("z_subject", mu=0, sigma=1, shape=n_subjects)
        subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)

        if use_tissue_effect:
            tissue_std = data.groupby("region_code")["y"].mean().std(ddof=1)
            sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=tissue_std)
            # Tissue Effect
            z_tissue = pm.Normal("z_tissue", mu=0, sigma=1, shape=n_tissues)
            tissue_effect = pm.Deterministic("tissue_effect", z_tissue * sigma_tissue)
            # Linear predictor
            mu = pm.Deterministic("mu", mu_intercept + beta_condition * condition +
                                  subject_effect[subject_idx] + tissue_effect[tissue_idx])
        else:
            mu = pm.Deterministic("mu", mu_intercept + beta_condition * condition +
                                  subject_effect[subject_idx])

        # ---- Likelihood Distribution from log normal (only positive values for each optical property) ----
        y = pm.LogNormal("y", mu=mu, sigma=sigma_residual, observed=y_obs)

        trace = pm.sample(draws=4000, tune=2000, target_accept=0.97, init="jitter+adapt_diag",
                          random_seed=42, progressbar=False)

    summary_vars = ["beta_condition", "mu_intercept", "sigma_subject", "sigma_residual"]
    if use_tissue_effect:
        summary_vars.append("sigma_tissue")

    # print the summary to the console
    summary_df = az.summary(trace, var_names=summary_vars, hdi_prob=0.95)
    print(f"\nSummary statistics for sheet '{sheet}', region '{region_str}', distance '{distance}':")
    print(summary_df.to_string())
    return trace, summary_df


def run_and_save_model(data, sheet, region_str, use_tissue_effect, writer, posterior_dir, distance):
    """Run PyMC model, save summary and posterior samples/plot for one sheet-region combination."""

    trace, summary_df = run_pymc_model(data, sheet, region_str, use_tissue_effect,distance)

    # Save summary to Excel writer
    summary_df.to_excel(writer, sheet_name=f"{sheet}_{region_str}")

    # Save trace to pickle
    trace_path = os.path.join(posterior_dir, f"{sheet}_{region_str}_dist_{distance}_trace.pkl")
    with open(trace_path, "wb") as f:
        pickle.dump(trace, f)
    print(f"Trace saved to {trace_path}")

    # Extract posterior samples of beta_condition as 1D array
    beta_samples = trace.posterior["beta_condition"].values.flatten()

    # Save posterior samples to CSV
    csv_path = os.path.join(posterior_dir, f"{sheet}_{region_str}_{distance}_posterior.csv")
    pd.DataFrame(beta_samples, columns=["beta_condition"]).to_csv(csv_path, index=False)

    # Plot posterior distribution
    az.plot_posterior(trace, var_names=["beta_condition"])
    plt.title(f"Posterior of beta_condition - {sheet} ({region_str})")
    png_path = os.path.join(posterior_dir, f"{sheet}_{region_str}_{distance}_posterior.png")
    plt.savefig(png_path)
    plt.close()

    print(f"Finished model for {sheet}, {region_str}, {distance} um."
          f"Results saved to:\n  {csv_path}\n  {png_path}")

    return beta_samples


def plot_overlay_posteriors(beta_samples_dict, sheet, posterior_dir, distance):
    """Plot overlay KDE plots for posterior samples from multiple regions."""

    plt.figure(figsize=(10, 6))
    for label, samples in beta_samples_dict.items():
        sns.kdeplot(samples, label=label, fill=True, alpha=0.4)

    plt.title(f"Overlay Posterior Distributions of beta_condition - {sheet} distance = {distance}")
    plt.xlabel("beta_condition")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()

    overlay_path = os.path.join(posterior_dir, f"{sheet}_distance_{distance}_overlay_posterior.png")
    plt.savefig(overlay_path)
    plt.close()
    print(f"\nOverlay posterior plot saved to {overlay_path}")

# ---- Main Code for Running Bayesian Stats ----
with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
    for sheet in sheet_names:
        full_data = pd.read_excel(
            '/projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_03Nov2025.xlsx',sheet_name=sheet)

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

        # Regions to analyze, with use_tissue_effect = False for single regions, True for combined
        region_list = [
            ("frontal", False),
            ("occipital", False),
            ("all", True)  # combined analysis over all
        ]

        for region, use_tissue in region_list:
            if region == "frontal":
                region_data = full_data[full_data["region"].str.lower() == "front"].copy()
            elif region == "occipital":
                region_data = full_data[full_data["region"].str.lower() == "occip"].copy()
            else:  # combined
                region_data = full_data.copy()

            if len(region_data) > 10:
                key_label = region.capitalize() if region != "all" else "All regions"
                print(f"\nRunning model for sheet '{sheet}' with region '{region_data['region'].unique()}'")
                beta_samples = run_and_save_model(region_data, sheet, region if region != "all" else "all",
                                                  use_tissue, writer, posterior_dir, distance)
                posterior_samples_dict[key_label] = beta_samples

        # Overlay plots comparing posteriors
        if posterior_samples_dict:
            plot_overlay_posteriors(posterior_samples_dict, sheet, posterior_dir, distance)