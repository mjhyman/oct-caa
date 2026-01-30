import matplotlib.pyplot as plt

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
    import numpy as np

    colors = plt.cm.tab10.colors
    plt.figure(figsize=(10, 6))

    if distance_grid is None:
        # take from first region mean
        first_region = next(iter(results_dict))
        distance_grid = np.arange(len(results_dict[first_region]["mean"]))  # fallback

    def highlight_significant_intervals(x_vals, sig_mask, color='yellow', alpha=0.15):
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

    for i, (region, region_data) in enumerate(results_dict.items()):
        mean = region_data["mean"]
        hpd = region_data["hpd"]
        sig = region_data["significant"]
        color = colors[i % len(colors)]

        plt.plot(distance_grid, mean, label=f"{region} mean", color=color)
        plt.fill_between(distance_grid, hpd[:, 0], hpd[:, 1], color=color, alpha=0.25)
        highlight_significant_intervals(distance_grid, sig, color=color, alpha=0.1)

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