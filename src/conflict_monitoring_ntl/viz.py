import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def plot_xarray_time_comparison(ds: xr.Dataset, band_name: str, cmap: str = "inferno"):
    """
    Plots all images for the specified band in the xarray dataset side by side or in a
    grid, automatically scaling each image using the 99th percentile of valid values.

    Args:
        ds: An xarray Dataset containing image data, with a 'time' dimension and bands.
        band_name: Name of the band variable in the Dataset to plot.
        cmap: Colormap for displaying the images. Defaults to "inferno".
    Returns:
        None. Displays the plot in matplotlib context.
    """
    date_dim = "time"
    max_cols = 4

    times = ds[date_dim].values
    images = ds[band_name]
    n = len(times)
    n_cols = min(n, max_cols)
    n_rows = int(np.ceil(n / max_cols))
    _, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = np.array(axes).reshape(-1)

    for i, ax in enumerate(axes[:n]):
        img = images.isel({date_dim: i}).values.squeeze()
        img = img / np.nanmax(img)  # normalize

        valid_vals = img[img > 0]
        vmax = np.percentile(valid_vals, 99) if valid_vals.size > 0 else 1

        if img.ndim == 3 and img.shape[0] == 3 and img.shape[-1] != 3:
            img = img.transpose(1, 2, 0)  # converts (3, h, w) → (h, w, 3)

        ax.imshow(img, cmap=cmap, vmin=img.min(), vmax=vmax)

        if isinstance(times[i], np.datetime64):
            label = times[i].astype("M8[D]").astype(object)
        else:
            label = str(times[i])
        ax.set_title(str(label))
        ax.axis("off")

    for ax in axes[n:]:
        ax.axis("off")

    plt.tight_layout()
    plt.show()
