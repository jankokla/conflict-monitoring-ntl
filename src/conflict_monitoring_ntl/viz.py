import warnings

import folium
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pycountry
import seaborn as sns
import xarray as xr
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.palettes import Set2
from bokeh.plotting import figure, show
from bokeh.transform import factor_cmap
from matplotlib.patches import Patch

from conflict_monitoring_ntl.config import PALETTE, SMOD_CLASS_ORDER
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


def plot_lighted_country_comparison_by_urbanisation(
    df,
    colors: tuple = ("#006666", "#B2D8D8"),
    labels: tuple = ("lighted", "non-lighted"),
    height: int = 4,
    aspect: float = 1.5,
):
    with sns.axes_style("whitegrid", {"grid.color": ".9", "grid.linestyle": ":"}):
        g = sns.FacetGrid(
            df, col="country_code", height=height, aspect=aspect, sharex=False
        )

        g.map_dataframe(
            sns.histplot,
            y="smod_class",
            hue="bm_binary",
            multiple="stack",
            stat="count",
            shrink=0.8,
            palette={boolean: color for boolean, color in zip([True, False], colors)},
        )

        for ax in g.axes.flat:
            ax.set_xlabel("Count of Pixels", labelpad=10)
            ax.set_ylabel("Degree of Urbanisation", labelpad=20)
            axis_country_code = ax.title._text.split("= ")[-1]
            ax.set_title(pycountry.countries.get(alpha_3=axis_country_code).name)

        labels_dict = {label: color for label, color in zip(labels, colors)}
        patches = [Patch(color=v, label=k) for k, v in labels_dict.items()]
        g.figure.legend(
            handles=patches,
            loc="upper left",
            ncol=len(patches),
            frameon=False,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.93])


def plot_lighted_proportion_comparison_by_urbanisation(
    df,
    col: str = "continent",
    colors: tuple = ("#383060", "#dad7eb"),
    labels: tuple = ("lighted", "non-lighted"),
    height: int = 4,
    aspect: float = 1.5,
):
    with sns.axes_style("whitegrid", {"grid.color": ".9", "grid.linestyle": ":"}):
        g = sns.FacetGrid(df, col=col, height=height, aspect=aspect)

        g.map_dataframe(
            sns.histplot,
            y="smod_class",
            hue="bm_binary",
            multiple="fill",
            stat="proportion",
            discrete=True,
            palette={boolean: color for boolean, color in zip([True, False], colors)},
            shrink=0.8,
        )

        for ax in g.axes.flat:
            ax.set_ylabel("Degree of Urbanisation", labelpad=20)

        labels_dict = {label: color for label, color in zip(labels, colors)}
        patches = [Patch(color=v, label=k) for k, v in labels_dict.items()]
        g.figure.legend(
            handles=patches,
            loc="upper left",
            ncol=len(patches),
            frameon=False,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.93])


def plot_world_coverage(country_codes: list[str]):
    from matplotlib.colors import ListedColormap

    world = gpd.read_file("https://datahub.io/core/geo-countries/r/countries.geojson")
    world["is_good"] = world["ISO3166-1-Alpha-3"].isin(country_codes)
    cmap = ListedColormap(["#bdbdbd", "#377eb8"])  # e.g.: red, blue

    _ = world.plot(
        column="is_good",
        cmap=cmap,
        legend=False,
        figsize=(17, 9),
        edgecolor="black",
        linewidth=0.3,
    )

    plt.axis("off")
    plt.show()


def plot_radiance_boxenplots(
    df: pd.DataFrame,
    vertical_line: int | None = None,
    height: int = 6,
    aspect: float = 1.8,
    show_labels: bool = True,
    show_legend: bool = True,
):
    with sns.axes_style("whitegrid", {"grid.color": ".9", "grid.linestyle": ":"}):
        g = sns.catplot(
            data=df,
            x="black_marble_radiance_monthly",
            y="smod_class",
            hue="continent",
            kind="boxen",
            order=SMOD_CLASS_ORDER,
            height=height,
            aspect=aspect,
            palette=PALETTE,
        )

        g._legend.remove()  # type: ignore
        if show_legend:
            g.add_legend(
                title="Continent",
                loc="upper center",
                bbox_to_anchor=(0.5, 1.05),
                ncol=len(df["continent"].unique()),
                frameon=False,
            )

        if show_labels:
            g.set_axis_labels(
                "Radiance - Black Marble (nW/cm²/sr) ",
                "Degree of Urbanisation",
                labelpad=20,
            )
        else:
            g.set_axis_labels("", "")

        if vertical_line:
            ax = g.axes.flat[0]
            ax.axvline(x=1, color="#DB8029", linestyle="--", linewidth=1)


def plot_scatter(
    df: pd.DataFrame,
    countries_to_plot: list[str],
    y: str = "f1",
    figsize: tuple = (12, 8),
    legend_loc: str = "upper left",
):
    df["is_plotted"] = False
    df.loc[df["country"].isin(countries_to_plot), "is_plotted"] = True

    plt.figure(figsize=figsize)
    with sns.axes_style("whitegrid", {"grid.color": ".9", "grid.linestyle": ":"}):
        ax = sns.scatterplot(
            data=df,
            x="hdi",
            y=y,
            hue="continent",
            palette=PALETTE,
            size="pixel_count",
            sizes=(20, 500),
            alpha=0.8,
        )

        # Axis labels
        ax.set_xlabel("Human Development Index", labelpad=20)
        ax.set_ylabel("F1", labelpad=20)

        leg = ax.get_legend()
        if leg:
            leg.remove()

        handles, labels = ax.get_legend_handles_labels()
        n_classes = df["continent"].nunique()
        ax.legend(
            handles[1 : n_classes + 1],
            labels[1 : n_classes + 1],
            title="Continent",
            loc=legend_loc,
            frameon=False,
        )

        for i, row in df.iterrows():
            if row.is_plotted:
                ax.annotate(
                    row["country"],
                    (row["hdi"], row["f1"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    arrowprops=dict(arrowstyle="->", lw=0.5, color="black"),
                    fontsize=9,
                )

        plt.tight_layout()
        plt.show()


def plot_scatter_bokeh(df: pd.DataFrame, y: str = "f1", legend_loc: str = "top_left"):
    continents = df["continent"].unique().tolist()
    palette = Set2[max(3, len(continents))]

    source = ColumnDataSource.from_df(df)

    fig = figure(x_axis_label="Human Development Index", y_axis_label=y, width=880)
    fig.scatter(
        x="hdi",
        y=y,
        source=source,
        alpha=0.8,
        color=factor_cmap("continent", palette=palette, factors=continents),
        legend_field="continent",
        size="size",
    )

    hover = HoverTool(
        tooltips=[
            ("Country", "@country"),
            ("F1", "@f1"),
            ("precision", "@precision"),
            ("recall", "@recall"),
            ("Pixel Count", "@pixel_count"),
        ]
    )

    fig.legend.location = legend_loc  # type: ignore

    fig.add_tools(hover)
    show(fig)


def plot_regression_radiance_by_country(df: pd.DataFrame, color: str):
    plt.figure(figsize=(12, 8))
    with sns.axes_style("whitegrid", {"grid.color": ".9", "grid.linestyle": ":"}):
        ax = sns.regplot(
            x="hdi_2020",
            y="black_marble_radiance_monthly",
            data=df,
            color=color,
        )

        ax.set_xlabel("Human Development Index (2020)", labelpad=20)
        ax.set_ylabel(
            "Mean Monthly Nighttime Radiance (nW/cm$^2$/sr) - (Jan 2020)", labelpad=20
        )

        for _, row in df.iterrows():
            if row.black_marble_radiance_monthly > 0.3:
                ax.annotate(
                    row["country_code"],
                    (row["hdi_2020"], row["black_marble_radiance_monthly"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    arrowprops=dict(arrowstyle="->", lw=0.5, color="black"),
                    fontsize=9,
                )

        plt.show()
