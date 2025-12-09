import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
from datetime import datetime as dt

# -----------------------
# Configuration
# -----------------------
# Optical properties to parse
sheet_names = ["scattering", "retardance"]
date = dt.now().strftime('%Y_%m_%d_%H%M')

# Directory for storing the beta spline model
posterior_dir = (f"/projectnb/npbssmic/ns/CAA/beta_stats_{date}/")
os.makedirs(posterior_dir, exist_ok=True)

# Outliers for each optical property
ulimit = {'scattering': 25,
          'retardance': 45}

# Initialize dictionary for storing summary stats
dfs_dict = {}

# Python Monte-Carlo (pymc) Sampling parameters
pymc_params = {
    'draws': 2000,          # typically 1000 - 2000
    'tune': 2000,           # typically 1000 - 2000
    'target_accept': 0.95,  # typically 0.9 - 0.99 for hierarchical/difficult parameters
    'max_treedepth': 15,    # up to 15 for difficult posteriors
    'random_seed': 42,      # any integer for reproducibility
    'progressbar': False}   # Set to false if running in script

def run_pymc_per_distance(data, sheet, region_str, posterior_dir, pymc_params, use_tissue_effect):
    """
    Run PyMC model at each discrete distance, comparing the two conditions.
    Returns a dataframe with posterior mean, HPD, significant flag for each distance.
    """
    # Print the optical property that's being tested
    print(f'\n##### STARTING {sheet} for {region_str} #####\n')

    # Get unique distances (sorted for nice plotting)
    distances = np.sort(data["distance"].unique())

    # Initialize lists
    results = []
    trace_paths = []

    for d in distances:
        # Print the distance
        print(f'\n## Distance {d} ##\n')

        # Create local copy at specific distance
        ddata = data[data['distance'] == d].copy()
        # Copy condition code
        ddata['condition_code'] = (ddata['condition'] == 'experimental').astype(int)

        # ---- Convert Categories to Codes ----
        # Condition (1 = EPVS. 0 = vessels)
        # condition_code = pd.Categorical(data["condition"]).codes

        # Identify subject and tissue indices
        subject_names = ddata['subject'].unique()
        tissue_names = ddata['region'].unique()
        subject_idx = ddata['subject'].map({name: i for i, name in enumerate(subject_names)}).values
        tissue_idx = ddata['region'].map({name: i for i, name in enumerate(tissue_names)}).values

        ### Measure standard deviations
        y_obs = ddata["y"].values
        subject_std = ddata.groupby('subject')["y"].mean().std(ddof=1)
        residual_std = (y_obs - ddata.groupby('condition')["y"].transform("mean")).std(ddof=1)

        with pm.Model() as model:
            # -----------------------
            # Distribution of priors for intercept, subject, residual
            # -----------------------
            mu_intercept = pm.Normal('mu_intercept', mu=np.mean(np.log(ddata['y'])), sigma=2)
            sigma_subject = pm.HalfNormal('sigma_subject', sigma=subject_std)
            sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

            # -----------------------
            # subject-level effects (non-centered parameterization)
            # -----------------------
            # Random effects of each subject (deviation from group average [units = Group Std. Dev.])
            z_subject = pm.Normal("z_subject", mu=0, sigma=1, shape=len(subject_names))
            # Scale each subject's deviation by the group-level Std. Dev. - more realistic effect sizes
            subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)

            # -----------------------
            # Beta conditions
            # -----------------------
            beta_condition = pm.Normal('beta_condition', mu=0, sigma=2)

            # -----------------------
            # Model predictor
            # -----------------------
            # use_tissue_effect is true when combining frontal and occipital into single vector
            if use_tissue_effect:
                tissue_std = ddata.groupby('region')["y"].mean().std(ddof=1)
                sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=tissue_std)
                z_tissue = pm.Normal("z_tissue", mu=0, sigma=1, shape=len(tissue_names))
                tissue_effect = pm.Deterministic("tissue_effect", z_tissue * sigma_tissue)
                # Linear Predictor
                mu = mu_intercept + beta_condition * ddata['condition_code'].values + subject_effect[subject_idx] + tissue_effect[tissue_idx]
            else:
                # Linear Predictor
                mu = mu_intercept + beta_condition * ddata['condition_code'].values + subject_effect[subject_idx]

            # -----------------------
            # Model the likelihood
            # -----------------------
            # degrees of freedom (normality) for StudentT distributions
            mu_inter = np.mean(np.log(y_obs))
            nu = pm.Exponential("nu", 1 / mu_inter)
            # Likelihood: StudentT on log-scale
            log_y_obs = np.log(y_obs)
            y = pm.StudentT("y", mu=mu, sigma=sigma_residual, nu=nu, observed=log_y_obs)

            # -----------------------
            # Sample the distribution "y" (time-consuming step)
            # -----------------------
            trace = pm.sample(draws=pymc_params['draws'],
                              tune=pymc_params['tune'],
                              target_accept=pymc_params['target_accept'],
                              max_treedepth=pymc_params['max_treedepth'],
                              random_seed=pymc_params['random_seed'],
                              progressbar=pymc_params['progressbar'])

            trace_path = os.path.join(posterior_dir, f'trace_distance_{d}.nc')
            az.to_netcdf(trace, trace_path)
            trace_paths.append(trace_path)
            posterior_summary = az.summary(trace, var_names=["beta_condition"], hdi_prob=0.95)
            results.append({
                "distance": d,
                "beta_mean": posterior_summary.loc["beta_condition", "mean"],
                "beta_sd": posterior_summary.loc["beta_condition", "sd"],
                "beta_hdi_lower": posterior_summary.loc["beta_condition", "hdi_2.5%"],
                "beta_hdi_upper": posterior_summary.loc["beta_condition", "hdi_97.5%"]})

    # Convert result list to dataframe
    result_df = pd.DataFrame(results)
    result_df.sort_values("distance", inplace=True)
    result_df.reset_index(drop=True, inplace=True)
    # Add significant to data fram
    result_df["significant"] = ((result_df["beta_hdi_lower"] > 0) | (result_df["beta_hdi_upper"] < 0))
    return result_df

def plot_per_distance_effect(result_df, label, savepath=None):
    plt.figure(figsize=(8, 5))
    plt.plot(result_df["distance"], result_df["beta_mean"], label=f"{label} mean", color="black")
    plt.fill_between(result_df["distance"], result_df["beta_hdi_lower"], result_df["beta_hdi_upper"],
                     color="gray", alpha=0.4)
    # Highlight significant intervals
    if "significant" in result_df.columns:
        sig_mask = result_df["significant"].values
        highlight_significant_intervals(result_df["distance"].values, sig_mask)
    plt.xlabel("Distance (micron)")
    plt.ylabel("Posterior mean effect (log(EPVS) - log(Vessels))")
    plt.title(f"{label}: Bayesian condition effect by distance")
    plt.legend()
    plt.tight_layout()
    if savepath is not None:
        plt.savefig(savepath, dpi=300)
        plt.close()
    else:
        plt.show()

def highlight_significant_intervals(x_vals, sig_mask, color='yellow', alpha=0.3):
    ax = plt.gca()
    in_sig = False
    start_idx = None
    for i, sig in enumerate(sig_mask):
        if sig and not in_sig:
            in_sig = True
            start_idx = i
        elif not sig and in_sig:
            in_sig = False
            ax.axvspan(x_vals[start_idx], x_vals[i-1], color=color, alpha=alpha)
    if in_sig:
        ax.axvspan(x_vals[start_idx], x_vals[-1], color=color, alpha=alpha)

def iqr_outlier_mask(g):
    """
    Perform 1.5*IQR outlier removal
    :param g: the group (region/subject) of data for 1.5*IQR removal
    :return: g after performing the 1.5*IQR removal
    """
    q1 = np.percentile(g['y'], 25)
    q3 = np.percentile(g['y'], 75)
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return (g['y'] >= lower) & (g['y'] <= upper)

# ---- Iterate over optical properties ----
for sheet in sheet_names:
    full_data = pd.read_excel(
        '/projectnb/npbssmic/ns/CAA/'
        'caa_all_radii_40um_donut_03Nov2025.xlsx',sheet_name=sheet)

    full_data.columns = full_data.columns.str.strip().str.lower()
    full_data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

    # -----------------------
    # Remove upper outliers, apply lower threshold, remove NAN ###
    # -----------------------
    # Remove upper limit outliers
    full_data = full_data[full_data["y"] < ulimit[sheet]]
    # Apply lower limit (0 for both mus and retardance)
    full_data = full_data[full_data["y"] > 0]
    # Ensure all are finite
    full_data = full_data[np.isfinite(full_data["y"])]
    # Remove NaNs
    full_data = full_data.dropna(subset=["y"])
    assert len(full_data) > 10, "Too few observations left after cleaning!"
    assert (full_data["y"] > 0).all(), "All values must be > 0 for LogNormal"

    # Apply 1.5*IQR outlier removal for each subject in each region
    mask = full_data.groupby(['region', 'subject'], group_keys=False).apply(iqr_outlier_mask)
    full_data = full_data[mask]
    print(f"\nSheet: {sheet} - {len(full_data)} samples after subject+region IQR outlier removal.")
    outlier_counts = full_data.groupby(['region', 'subject'])['y'].count()
    print(outlier_counts)

    ### Run MCMC model at each distance for each region
    # Frontal
    print(f"Running Bayesian model for each distance, frontal region.")
    frontal_data = full_data[full_data["region"].str.lower() == "front"].copy()
    frontal_result_df = run_pymc_per_distance(frontal_data, sheet,
                                              'frontal', posterior_dir, pymc_params,
                                              False)
    # Occipital
    print(f"Running Bayesian model for each distance, occipital region.")
    occipital_data = full_data[full_data["region"].str.lower() == "occip"].copy()
    occipital_result_df = run_pymc_per_distance(occipital_data, sheet,
                                                'occipital', posterior_dir, pymc_params,
                                                False)
    # Combined
    print(f"Running Bayesian model for each distance, all regions combined.")
    combined_result_df = run_pymc_per_distance(full_data, sheet,
                                               'combined', posterior_dir, pymc_params,
                                               True)

    # Save summary stats to dictionary
    dfs_dict[f"frontal_{sheet}"] = frontal_result_df
    dfs_dict[f"occipital_{sheet}"] = occipital_result_df
    dfs_dict[f"combined_{sheet}"] = combined_result_df

    ### Plot data
    plot_per_distance_effect(frontal_result_df, "Frontal",
                                 savepath=os.path.join(posterior_dir, f"{sheet}_frontal_per_distance.png"))
    plot_per_distance_effect(occipital_result_df, "Occipital",
                                 savepath=os.path.join(posterior_dir, f"{sheet}_occipital_per_distance.png"))
    plot_per_distance_effect(combined_result_df, "Combined",
                                 savepath=os.path.join(posterior_dir, f"{sheet}_combined_per_distance.png"))

# ----------------------
# Save summary statistics to spreadsheet
# ----------------------
# Excel filepath for summary statistics
excel_path = os.path.join(posterior_dir, "all_properties_per_distance.xlsx")
with pd.ExcelWriter(excel_path) as writer:
    for sheet_name, df in dfs_dict.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)

# ----------------------
# Save a metadata text file to posterior_dir
# ----------------------
# Filepath for meta data
metadata_path = os.path.join(posterior_dir, "metadata.txt")
# Write meta data
with open(metadata_path, "w") as metafile:
    metafile.write("Bayesian Condition Effect Model Metadata\n")
    metafile.write(f"Run date: {date}\n")
    metafile.write(f"Input sheets: {sheet_names}\n")
    metafile.write("Input file: /projectnb/npbssmic/ns/CAA/caa_all_radii_40um_donut_03Nov2025.xlsx\n")
    metafile.write("Outlier upper limits:\n")
    for sheet, limit in ulimit.items():
        metafile.write(f"  {sheet}: {limit}\n")
    metafile.write("PyMC sampling parameters:\n")
    for k, v in pymc_params.items():
        metafile.write(f"  {k}: {v}\n")
    metafile.write("Analysis regions: frontal, occipital, combined\n")
    metafile.write(f"Summary Excel saved to: {excel_path}\n")
    for key in dfs_dict:
        png_name = f"{key.replace('_', '_')}_per_distance.png"
        metafile.write(f"Plot PNG: {os.path.join(posterior_dir, png_name)}\n")
    # Optionally, list traces and any other outputs as needed
metafile.close()
print(f"Metadata file saved to {metadata_path}")