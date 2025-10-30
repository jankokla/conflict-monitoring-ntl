import warnings

import folium
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from conflict_monitoring_ntl.utils import binarize_xarray


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


def plot_map_with_shape(
    gdf: gpd.GeoDataFrame, zoom_start: int = 10, is_layer_control: bool = True
) -> folium.Map:
    """
    Plots geographic shapes from a GeoDataFrame on an interactive Folium map
        with ESRI satellite and label layers.

    Args:
        gdf: A GeoDataFrame containing the shape(s) to display.
        zoom_start: The initial zoom level for the map. Defaults to 10.
        is_layer_control: If True, adds a layer control widget so users can
            toggle visibility of basemaps and overlays (default True).

    Returns:
        A folium.Map instance with the shape overlay, ESRI base layers,
            and (optionally) interactive layer controls.
    """
    warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS.*")

    center = [gdf.geometry.centroid.y.mean(), gdf.geometry.centroid.x.mean()]

    # tiles="CartoDB Positron"
    m = folium.Map(location=center, zoom_start=zoom_start, tiles=None)

    # add base layer from ESRI
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",  # noqa: E501
        attr="Esri",
        name="Esri Satellite",
        overlay=False,
        control=True,
    ).add_to(m)

    # add labels
    folium.TileLayer(
        tiles="https://services.arcgisonline.com/arcgis/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile/{z}/{y}/{x}",  # noqa: E501
        attr="Esri",
        name="Labels",
        overlay=True,
        control=True,
    ).add_to(m)

    folium.GeoJson(
        gdf.to_json(),
        name="Boundary",
        style_function=lambda _: {
            "fillColor": "none",
            "color": "red",
            "weight": 3,
            "fillOpacity": 0,
        },
    ).add_to(m)

    if is_layer_control:
        folium.LayerControl().add_to(m)

    return m


def plot_admin_map_with_tiles(
    country_gdf: gpd.GeoDataFrame,
    raster_gdf: gpd.GeoDataFrame,
    admin_gdf: gpd.GeoDataFrame,
    base_map: str = "CartoDB positron",
    zoom_start: int = 7,
) -> folium.Map:
    # TODO: add docstring
    centroid = country_gdf.geometry[0].centroid
    m = folium.Map(
        location=[centroid.y, centroid.x], zoom_start=zoom_start, tiles=base_map
    )

    folium.GeoJson(
        data=raster_gdf.geometry,
        name="Raster coverage",
        style_function=lambda x: {"fillColor": "black", "color": "black", "weight": 1},
    ).add_to(m)

    for _, row in admin_gdf.iterrows():
        color = "#07407B" if row.is_within_raster else "#F7931F"
        folium.GeoJson(
            data=row.geometry,
            name="Districts",
            style_function=lambda x, color=color: {
                "fillColor": color,
                "color": color,
                "weight": 1,
            },
        ).add_to(m)

    folium.GeoJson(
        data=country_gdf.iloc[0].geometry,
        style_function=lambda x: {"fillColor": "none", "color": "grey", "weight": 2},
    ).add_to(m)

    return m


def plot_binary(arr: xr.DataArray, satellite: str, threshold: float):
    """Helper function for plotting binarized tiles."""
    arr_binary = binarize_xarray(arr, threshold)

    return arr_binary.hvplot.image(
        x="lon",
        y="lat",
        geo=True,
        crs=arr.rio.crs,
        cmap="inferno",
        clim=(0, 1),
        title=f"{satellite} Binarized",
        colorbar=False,
    )


def plot_tile_comparison(
    arr_left,
    arr_right,
    title_left="Left",
    title_right="Right",
    clim_left=(0, 1),
    clim_right=(0, 1),
    cmap_left="inferno",
    cmap_right="inferno",
    colorbar=True,
    extra_plot_kwargs={},
):
    import hvplot.xarray

    left_plot = arr_left.hvplot.image(
        x="lon",
        y="lat",
        crs=arr_left.rio.crs,
        clim=clim_left,
        cmap=cmap_left,
        title=title_left,
        geo=True,
        colorbar=colorbar,
        **extra_plot_kwargs,
    )
    right_plot = arr_right.hvplot.image(
        x="lon",
        y="lat",
        crs=arr_right.rio.crs,
        clim=clim_right,
        cmap=cmap_right,
        title=title_right,
        geo=True,
        colorbar=colorbar,
        **extra_plot_kwargs,
    )

    return left_plot + right_plot
