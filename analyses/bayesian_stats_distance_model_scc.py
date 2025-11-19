import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
from datetime import datetime
from patsy import dmatrix

# -----------------------
# Configuration
# -----------------------
# Optical properties to parse
sheet_names = ["scattering", "retardance"]

# Directory for storing the beta spline model
posterior_dir = "/projectnb/npbssmic/ns/CAA/beta_spline/"
os.makedirs(posterior_dir, exist_ok=True)

# Output file path
date = datetime.now().strftime("%Y_%m_%d")
output_excel_path = ("/projectnb/npbssmic/ns/CAA/beta_spline/"
                     f"beta_spline_{date}.xlsx")

# Outliers for each optical property
ulimit = {'scattering': 25,
          'retardance': 45}

# Python Monte-Carlo (pymc) Sampling parameters
pymc_params = {
    'draws': 1000,          # typically 1000 - 2000
    'tune': 1000,           # typically 1000 - 2000
    'target_accept': 0.8,  # typically 0.9 - 0.99 for hierarchical/difficult parameters
    'max_treedepth': 15,    # up to 15 for difficult posterios
    'random_seed': 42,      # any integer for reproducibility
    'progressbar': False}   # Set to false if running in script

def run_pymc_model(data, sheet, region_str, stats_dir, pymc_params, use_tissue_effect=False):
    """
    Define  PYMC model based on data, perform random sampling, generate statistics
    :param data: Dataframe subset of spreadsheet
    :param sheet: the optical property
    :param region_str: the brain region (front or occip)
    :param stats_dir: directory to output the PYMC model as a pkl file
    :param pymc_params: dictionary of pymc values
    :param use_tissue_effect: boolean for combining frontal and occipital
    :return: Does not return any variables. saves the model to the stats_dir
    """
    # --- Preprocess ---
    data.columns = data.columns.str.strip().str.lower()
    data.rename(columns={"groups": "condition", "subid": "subject", "opticalproperty": "y"}, inplace=True)

    # Extract Values from dataframes
    y_obs = data["y"].values
    distance_vals = data["distance"].values
    subject_idx = data["subject"].values

    # ---- Scale/standardize distance to improve parameter scale + numerical stability
    distance_mean = data["distance"].mean()
    distance_std = data["distance"].std()
    distance_scaled = (data["distance"].values - distance_mean) / distance_std

    # ---- Convert Categories to Codes ----
    # Condition (1 = EPVS. 0 = vessels)
    condition_code = pd.Categorical(data["condition"]).codes
    tissue_idx = pd.Categorical(data["region"]).codes

    # ----  Retrieve number of subjects and tissue regions ----
    n_subjects = len(np.unique(subject_idx))
    n_tissues = len(np.unique(tissue_idx))

    # --- Spline basis for nonlinear distance effect ---
    DF_SPLINE = 4
    spline_basis = dmatrix(f"bs(distance_scaled, df={DF_SPLINE}, degree=3, include_intercept=False) - 1",
                           {"distance_scaled": distance_scaled}, return_type='dataframe')
    X_base = spline_basis.values
    n_splines = X_base.shape[1]

    # Duplicate basis for two conditions (interaction model)
    X_cond0 = X_base * (1 - condition_code)[:, None]
    X_cond1 = X_base * condition_code[:, None]
    X_spline = np.hstack([X_cond0, X_cond1])  # shape: (n_samples, 2 * n_splines)

    # --- Priors based on observed data ---
    mu_inter = np.mean(np.log(y_obs))
    sigma_inter = np.std(np.log(y_obs))

    subject_std = data.groupby(subject_idx)["y"].mean().std(ddof=1)
    residual_std = (y_obs - data.groupby(condition_code)["y"].transform("mean")).std(ddof=1)

    with pm.Model() as model:
        # -----------------------
        # Distribution of priors for intercept, subject, residual
        # -----------------------
        mu_intercept = pm.Normal("mu_intercept", mu=mu_inter, sigma=sigma_inter)
        sigma_subject = pm.HalfNormal("sigma_subject", sigma=subject_std)
        sigma_residual = pm.HalfNormal("sigma_residual", sigma=residual_std)

        # -----------------------
        # subject-level effects (non-centered parameterization)
        # -----------------------
        # Random effects of each subject (deviation from group average [units = Group Std. Dev.])
        z_subject = pm.Normal("z_subject", mu=0, sigma=1, shape=n_subjects)
        # Scale each subject's deviation by the group-level Std. Dev. - more realistic effect sizes
        subject_effect = pm.Deterministic("subject_effect", z_subject * sigma_subject)

        # -----------------------
        # Spline coefficients (2 × n_splines for 2 conditions)
        # -----------------------
        # Regression coefficients (n_splines is the degrees of freedom). The first n_splines values are the coefficients
        #  for EPVS, and the subsequent n_splines values are the coefficients for the vessels.
        beta_spline = pm.Normal("beta_spline", mu=0, sigma=0.5, shape=2 * n_splines)
        # Difference in spline weights [experimental (EPVS, condition 1) - control (vessels, condition 0)]
        # larger difference b/w weights indicates larger difference b/w groups
        beta_diff = pm.Deterministic("beta_diff", beta_spline[n_splines:] - beta_spline[:n_splines])
        # Non-linear effect for each coefficient (reflects actual value of spline curve on log scale)
        spline_term = pm.math.dot(X_spline, beta_spline)

        # -----------------------
        # Account for combining frontal and occipital
        # -----------------------
        # use_tissue_effect is true when combining frontal and occipital into single vector
        if use_tissue_effect:
            tissue_std = data.groupby(tissue_idx)["y"].mean().std(ddof=1)
            sigma_tissue = pm.HalfNormal("sigma_tissue", sigma=tissue_std)
            z_tissue = pm.Normal("z_tissue", mu=0, sigma=1, shape=n_tissues)
            tissue_effect = pm.Deterministic("tissue_effect", z_tissue * sigma_tissue)
            # Linear Predictor
            mu = mu_intercept + spline_term + subject_effect[subject_idx] + tissue_effect[tissue_idx]
        else:
            # Linear Predictor
            mu = mu_intercept + spline_term + subject_effect[subject_idx]

        # -----------------------
        # Model the likelihood
        # -----------------------
        # degrees of freedom (normality) for StudentT distributions
        nu = pm.Exponential("nu", 1 / mu_inter)
        # Likelihood: StudentT on log-scale
        log_y_obs = np.log(y_obs)
        y = pm.StudentT("y", mu=mu, sigma=sigma_residual, nu=nu, observed=log_y_obs)

        # -----------------------
        # Sample the distribution "y" (time-consuming step)
        # -----------------------
        print(f"\nSampling the MCMC distribution.")
        trace = pm.sample(draws=pymc_params['draws'], tune=pymc_params['tune'],
                          target_accept=pymc_params['target_accept'],
                          max_treedepth=pymc_params['max_treedepth'],
                          random_seed=pymc_params['random_seed'],
                          progressbar=pymc_params['progressbar'])

    # ---- Print the summary statistics to console (these are saved in log) ----
    var_names = ["beta_spline", "beta_diff", "mu_intercept", "sigma_subject", "sigma_residual"]
    if use_tissue_effect:
        var_names.append("sigma_tissue")
    # Create printable summary of results
    summary_df = az.summary(trace, var_names=var_names)
    print(f"\nSummary statistics for '{sheet}' and region '{region_str}':")
    print(summary_df.to_string())

    # ---- Save outputs to spreadsheets and pickle objects----
    safe_sheet = sheet.replace(" ", "_").lower()
    safe_region = region_str.replace(" ", "_").lower()
    suffix = f"{safe_sheet}_{safe_region}"
    outdir = os.path.join(stats_dir,f"spline_model_outputs/{suffix}")
    os.makedirs(outdir, exist_ok=True)
    # Save summary to spreadsheet
    with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name=f"{suffix}")
    # Save trace
    with open(os.path.join(outdir, f"trace_{suffix}.pkl"), "wb") as f:
        pickle.dump(trace, f)
    # Save summary as pickle
    with open(os.path.join(outdir, f"summary_{suffix}.pkl"), "wb") as f:
        pickle.dump(summary_df, f)
    # Save X_spline as pickle
    with open(os.path.join(outdir, f"xspline_{suffix}.pkl"), "wb") as f:
        pickle.dump(X_spline, f)
    # Save spline_basis as pickle (preserves DataFrame structure and column names)
    with open(os.path.join(outdir, f"spline_basis_{suffix}.pkl"), "wb") as f:
        pickle.dump(spline_basis, f)

    # ---- Plot the spline fit ----
    # This can also be called after generating the model from a batch script and using the separate plotting script
    try:
        # Extract posterior means for spline coefficients
        beta_spline_mean = trace.posterior["beta_spline"].mean(dim=["chain", "draw"]).values
        n_splines = spline_basis.shape[1]

        # Reconstruct splines for each condition
        spline_vals_cond0 = spline_basis.values @ beta_spline_mean[:n_splines]
        spline_vals_cond1 = spline_basis.values @ beta_spline_mean[n_splines:]

        # Sort the distance values (X-axis) for plotting
        sort_idx = np.argsort(distance_vals)

        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(distance_vals[sort_idx], spline_vals_cond0[sort_idx], label="Condition 0", color="blue")
        plt.plot(distance_vals[sort_idx], spline_vals_cond1[sort_idx], label="Condition 1", color="red")
        plt.xlabel("Distance (micron)")
        plt.ylabel("Spline fit (a.u.)")
        plt.title(f"Spline Fit Comparison: {sheet}, {region_str}")
        plt.legend()
        plt.tight_layout()

        # Save plot
        plot_path = os.path.join(outdir, f"spline_plot_{suffix}.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
    except Exception as e:
        print(f"Plotting failed for {sheet}, {region_str}: {e}")

    # ---- Plot difference spline curve and highlight significant distances ----
    try:
        # Create fine grid for distance
        distance_grid = np.linspace(distance_vals.min(), distance_vals.max(), 100)
        spline_basis_grid = dmatrix(
            f"bs(distance_grid, df={DF_SPLINE}, degree=3, include_intercept=False) - 1",
            {"distance_grid": distance_grid}, return_type='dataframe')
        X_grid = spline_basis_grid.values  # (100, n_splines)

        # Extract posterior beta_diff samples: combine chains and draws
        beta_diff_samples = trace.posterior["beta_diff"].stack(sample=("chain", "draw")).values.T  # (samples, n_splines)

        # Compute difference function posterior samples at grid points: (samples, 100)
        diff_posterior = beta_diff_samples @ X_grid.T

        # Mean and 95% credible interval
        diff_mean = np.mean(diff_posterior, axis=0)
        diff_hpd = az.hdi(diff_posterior.T, hdi_prob=0.95)  # shape (100, 2)

        # Identify significant bins where HPD excludes zero
        significant = (diff_hpd[:, 0] > 0) | (diff_hpd[:, 1] < 0)

        # Plot difference curve with credible interval
        plt.figure(figsize=(8, 5))
        plt.plot(distance_grid, diff_mean, color='black', label="Difference (Cond1 - Cond0)")
        plt.fill_between(distance_grid, diff_hpd[:, 0], diff_hpd[:, 1], color='gray', alpha=0.3, label='95% Credible Interval')

        # Highlight significant intervals
        highlight_significant_intervals(distance_grid, significant)

        plt.xlabel("Distance (micron)")
        plt.ylabel("Difference in spline fit (a.u.)")
        plt.title(f"Difference Between Conditions with Significant Regions\n{sheet}, {region_str}")
        plt.legend()
        plt.tight_layout()

        diff_plot_path = os.path.join(outdir, f"difference_spline_{suffix}.png")
        plt.savefig(diff_plot_path, dpi=300)
        plt.close()

        print(f"Saved difference spline plot with significant regions to {diff_plot_path}")
    except Exception as e:
        print(f"Failed to plot difference spline for {sheet}, {region_str}: {e}")

    return diff_mean, diff_hpd, significant, distance_grid

def plot_overlay_difference_splines(results_dict,
                                    distance_grid=None,
                                    title="Overlay of Difference Spline Posteriors",
                                    xlabel="Distance (micron)",
                                    ylabel="Difference (Cond1 - Cond0)",
                                    savepath=None):
    """
    Plot overlay of difference spline mean and credible intervals for multiple regions.

    Parameters
    ----------
    results_dict : dict
        Dictionary keyed by region name, each with a dict containing:
          - "mean": 1D numpy array of mean difference curve over distance_grid
          - "hpd": 2D numpy array shape (len(distance_grid), 2) with lower and upper HPD bounds
          - "significant": boolean 1D numpy array (len(distance_grid)) indicating significant points

    distance_grid : 1D numpy array, optional
        Common x-axis values (distance). If None, keys will be taken from first region.

    title : str, optional
        Plot title.

    xlabel, ylabel : str, optional
        Axis labels.

    savepath : str or None, optional
        If provided, save the plot PNG to this path.

    Returns
    -------
    None
    """

    colors = plt.cm.tab10.colors
    plt.figure(figsize=(10, 6))

    if distance_grid is None:
        # take from first region mean
        first_region = next(iter(results_dict))
        distance_grid = np.arange(len(results_dict[first_region]["mean"]))  # fallback

    for i, (region, region_data) in enumerate(results_dict.items()):
        mean = region_data["mean"]
        hpd = region_data["hpd"]
        sig = region_data["significant"]
        color = colors[i % len(colors)]
        # Create Plot
        plt.plot(distance_grid, mean, label=f"{region} mean", color=color)
        plt.fill_between(distance_grid, hpd[:, 0], hpd[:, 1], color=color, alpha=0.25)
        highlight_significant_intervals(distance_grid, sig, color=color, alpha=0.1)
    # Configure figure
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=300)
        plt.close()
        print(f"Overlay difference plot saved to: {savepath}")
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

    # Apply mask for each subject-region combo
    mask = full_data.groupby(['region', 'subject'], group_keys=False).apply(iqr_outlier_mask)
    full_data = full_data[mask]
    print(f"\nSheet: {sheet} - {len(full_data)} samples after subject+region IQR outlier removal.")
    outlier_counts = full_data.groupby(['region', 'subject'])['y'].count()
    print(outlier_counts)

    # Frontal
    print(f"Running model for sheet '{sheet}' with frontal region only.")
    frontal_data = full_data[full_data["region"].str.lower() == "front"].copy()
    if len(frontal_data) > 10:
        diff_mean_frontal, diff_hpd_frontal, significant_frontal, _ = \
            run_pymc_model(frontal_data, sheet, 'frontal',
                           posterior_dir, pymc_params, False)

    # Occipital
    print(f"Running model for sheet '{sheet}' with occipital region only.")
    occipital_data = full_data[full_data["region"].str.lower() == "occip"].copy()
    if len(occipital_data) > 10:
        diff_mean_occipital, diff_hpd_occipital, significant_occipital, _ = \
            run_pymc_model(occipital_data, sheet, 'occipital',
                           posterior_dir, pymc_params,False)

    # All regions combined
    print(f"\nRunning model for sheet '{sheet}' with all regions combined.")
    diff_mean_combined, diff_hpd_combined, significant_combined, distance_grid = \
        run_pymc_model(full_data.copy(), sheet, 'combined',
                       posterior_dir, pymc_params, True)

    print(f"\nCreating figure with all regions combined.")
    # ---- Overlay plot ----
    results_to_plot = {
        "Frontal": {
            "mean": diff_mean_frontal,
            "hpd": diff_hpd_frontal,
            "significant": significant_frontal
        },
        "Occipital": {
            "mean": diff_mean_occipital,
            "hpd": diff_hpd_occipital,
            "significant": significant_occipital
        },
        "Combined": {
            "mean": diff_mean_combined,
            "hpd": diff_hpd_combined,
            "significant": significant_combined
        }
    }

    plot_overlay_difference_splines(results_to_plot,
                                    distance_grid=distance_grid,
                                    title="Difference Between Conditions by Region",
                                    savepath=os.path.join(posterior_dir,
                                                          f"spline_model_outputs/{sheet}_regions_posteriors.png"))
