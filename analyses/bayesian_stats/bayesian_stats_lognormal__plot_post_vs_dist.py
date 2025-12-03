import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})

def parse_sheet_name(sheet_name):
    """
    Parse sheet/tab name like 'scattering_frontal' into (sheet, region).
    """
    parts = sheet_name.split('_')
    if len(parts) < 2:
        return None, None
    sheet = parts[0]
    region = parts[1]
    return sheet, region


def load_posteriors_from_excels(posterior_dir):
    """
    Loads mean & HDI for beta_condition from each Excel file (one distance),
    each with multiple sheet tabs named as '<sheet>_<region>'.
    Returns combined DataFrame with columns:
    sheet, region, distance, median, lower, upper
    """
    files = glob.glob(os.path.join(posterior_dir, "*.xlsx"))
    records = []
    for f in files:
        # Extract distance from filename, e.g. 'posterior_dist_40.xlsx'
        m = re.search(r'(\d+)', os.path.basename(f))
        if not m:
            print(f"Warning: Could not parse distance from filename {f}")
            continue
        distance = int(m.group(1))

        # Open all sheets in this Excel file
        xls = pd.ExcelFile(f)
        for sheetname in xls.sheet_names:
            sheet, region = parse_sheet_name(sheetname)
            if sheet is None or region is None:
                continue

            df = xls.parse(sheetname, index_col=0)
            # Confirm the variable 'beta_condition' is in the index
            if 'beta_condition' not in df.index:
                print(f"Warning: 'beta_condition' missing in {sheetname} in {f}")
                continue

            # Extract mean and HDI bounds: adjust these column names if needed
            try:
                mean_val = df.loc['beta_condition', 'mean']
                lower = df.loc['beta_condition', 'hdi_2.5%']
                upper = df.loc['beta_condition', 'hdi_97.5%']
            except KeyError:
                # Try alternative HDI names
                lower = df.loc['beta_condition', 'hdi_lower'] if 'hdi_lower' in df.columns else None
                upper = df.loc['beta_condition', 'hdi_upper'] if 'hdi_upper' in df.columns else None
                mean_val = df.loc['beta_condition', 'mean'] if 'mean' in df.columns else None
                if None in (mean_val, lower, upper):
                    print(f"Warning: Missing mean or HDI columns in {sheetname} in {f}")
                    continue

            records.append({
                "sheet": sheet,
                "region": region,
                "distance": distance,
                "median": mean_val,
                "lower": lower,
                "upper": upper
            })
    return pd.DataFrame(records)


# The plotting functions remain basically the same as provided before, just slightly adapted for column 'sheet' and 'region'

def plot_region_fill(df_region, region_name, sheet_name, output_dir, ylims):
    plt.figure(figsize=(8, 5))
    df_region = df_region.sort_values("distance")
    plt.fill_between(df_region["distance"], df_region["lower"], df_region["upper"],
                     color='lightblue', alpha=0.5, label="95% credible interval")
    plt.plot(df_region["distance"], df_region["median"], 'o-', color='blue', label="Mean estimate")
    plt.axhline(0, linestyle='--', color='gray', linewidth=0.7)
    plt.ylim(ylims)
    plt.xlabel("Distance ($\\mu$m)")
    plt.ylabel("Percentage Change (%)")
    if region_name == 'all':
        tstr = 'Frontal + Occipital'
        plt.title(f"Effect Size - {sheet_name.capitalize()} - {tstr}")
    else:
        plt.title(f"Effect Size - {sheet_name.capitalize()} - {region_name.capitalize()}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    filename = os.path.join(output_dir, f"{sheet_name}_{region_name}_fill_between.png")
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved {filename}")


def plot_region_errorbar(df_region, region_name, sheet_name, output_dir, ylims):
    plt.figure(figsize=(8, 5))
    df_region = df_region.sort_values("distance")
    yerr_lower = df_region["median"] - df_region["lower"]
    yerr_upper = df_region["upper"] - df_region["median"]
    plt.errorbar(df_region["distance"], df_region["median"],
                 yerr=[yerr_lower, yerr_upper], fmt='o-', capsize=5, color='green',
                 label="Mean with 95% CI")
    plt.axhline(0, linestyle='--', color='gray', linewidth=0.7)
    plt.ylim(ylims)
    plt.xlabel("Distance ($\\mu$m)")
    plt.ylabel("Percentage Change (%)")
    if region_name == 'all':
        tstr = 'Frontal + Occipital'
        plt.title(f"Effect Size - {sheet_name.capitalize()} - {tstr}")
    else:
        plt.title(f"Effect Size - {sheet_name.capitalize()} - {region_name.capitalize()}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    filename = os.path.join(output_dir, f"{sheet_name}_{region_name}_errorbar.png")
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved {filename}")


def plot_overlay_fill(df, sheet_name, output_dir, ylims):
    plt.figure(figsize=(10, 6))
    colors = {'frontal': '#1f77b4', 'occipital': '#ff7f0e', 'all': '#2ca02c'}
    for region in ['frontal', 'occipital', 'all']:
        df_region = df[df["region"] == region].sort_values("distance")
        if df_region.empty:
            continue
        plt.fill_between(df_region["distance"], df_region["lower"], df_region["upper"],
                         color=colors.get(region, 'gray'), alpha=0.3, label=f"{region.capitalize()} 95% CI")
        plt.plot(df_region["distance"], df_region["median"], 'o-', label=f"{region.capitalize()} Mean",
                 color=colors.get(region, 'gray'))
    plt.axhline(0, color='black', linestyle='--', linewidth=0.7)
    plt.ylim(ylims)
    plt.xlabel("Distance ($\\mu$m)")
    plt.ylabel("Percentage Change (%)")
    plt.title(f"Effect Size - {sheet_name.capitalize()} - Overlay Regions")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    filename = os.path.join(output_dir, f"{sheet_name}_overlay_fill_between.png")
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved {filename}")


def plot_overlay_errorbar(df, sheet_name, output_dir, ylims):
    plt.figure(figsize=(10, 6))
    colors = {'frontal': '#1f77b4', 'occipital': '#ff7f0e', 'all': '#2ca02c'}
    for region in ['frontal', 'occipital', 'all']:
        df_region = df[df["region"] == region].sort_values("distance")
        if df_region.empty:
            continue
        yerr_lower = df_region["median"] - df_region["lower"]
        yerr_upper = df_region["upper"] - df_region["median"]
        plt.errorbar(df_region["distance"], df_region["median"],
                     yerr=[yerr_lower, yerr_upper], fmt='o-', capsize=5,
                     label=f"{region.capitalize()}", color=colors.get(region, 'gray'))
    plt.axhline(0, color='black', linestyle='--', linewidth=0.7)
    plt.ylim(ylims)
    plt.xlabel("Distance ($\\mu$m)")
    plt.ylabel("Percentage Change (%)")
    plt.title(f"Effect Size - {sheet_name.capitalize()} - Overlay Regions")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    filename = os.path.join(output_dir, f"{sheet_name}_overlay_errorbar.png")
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved {filename}")


def main():
    posterior_dir = "/projectnb/npbssmic/ns/CAA/beta_spline/"
    output_dir = "/projectnb/npbssmic/ns/CAA/beta_spline/figures/"
    os.makedirs(output_dir, exist_ok=True)
    sheet_names = ["scattering", "retardance"]

    for sheet_name in sheet_names:
        print(f"Processing sheet: {sheet_name}")
        df = load_posteriors_from_excels(posterior_dir)

        if df.empty:
            print(f"No data found for sheet '{sheet_name}', skipping.")
            continue

        # Filter for current sheet
        df_sheet = df[df['sheet'] == sheet_name].copy()

        if df_sheet.empty:
            print(f"No data for sheet {sheet_name}, skipping.")
            continue

        # Convert from log scale to % change
        df_sheet["median"] = 100 * (np.exp(df_sheet["median"]) - 1)
        df_sheet["lower"] = 100 * (np.exp(df_sheet["lower"]) - 1)
        df_sheet["upper"] = 100 * (np.exp(df_sheet["upper"]) - 1)

        ylims_dict = {
            "scattering": [-2.5, 4],
            "retardance": [-7.5, 7],
            "orientation": [-60, 125]
        }
        ylims = ylims_dict.get(sheet_name, None)

        for region in ['all', 'frontal', 'occipital']:
            df_region = df_sheet[df_sheet['region'] == region]
            if df_region.empty:
                print(f"No data for region '{region}' in sheet '{sheet_name}'")
                continue

            plot_region_fill(df_region, region, sheet_name, output_dir, ylims)
            plot_region_errorbar(df_region, region, sheet_name, output_dir, ylims)

        plot_overlay_fill(df_sheet, sheet_name, output_dir, ylims)
        plot_overlay_errorbar(df_sheet, sheet_name, output_dir, ylims)


if __name__ == "__main__":
    main()