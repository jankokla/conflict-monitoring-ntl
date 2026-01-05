import warnings

import folium
import geopandas as gpd
import matplotlib.pyplot as plt
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


def plot_map_with_shape(
    gdf: gpd.GeoDataFrame, zoom_start: int = 10, is_layer_control: bool = True
) -> folium.Map:
    """
    Plot geographic shapes on an interactive Folium map with ESRI satellite layers.

    Args:
        gdf: A GeoDataFrame containing the shape(s) to display.
        zoom_start: The initial zoom level for the map.
        is_layer_control: If True, adds a widget to toggle basemaps and overlays.

    Returns:
        A folium.Map instance with shape overlay and ESRI base layers.
    """
    warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS.*")

    center = [gdf.geometry.centroid.y.mean(), gdf.geometry.centroid.x.mean()]

    m = folium.Map(location=center, zoom_start=zoom_start, tiles=None)

    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Esri",
        name="Esri Satellite",
        overlay=False,
        control=True,
    ).add_to(m)

    folium.TileLayer(
        tiles="https://services.arcgisonline.com/arcgis/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile/{z}/{y}/{x}",
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
    """
    Visualize administrative districts and their intersection with raster coverage.

    Districts are color-coded based on whether they fall within the bounds of 
    the available satellite raster data.

    Args:
        country_gdf: GeoDataFrame representing the national boundary.
        raster_gdf: GeoDataFrame representing the footprint of satellite data.
        admin_gdf: GeoDataFrame of administrative districts with 'is_within_raster' column.
        base_map: The folium/leaflet tile provider string.
        zoom_start: Initial zoom level.

    Returns:
        A folium.Map with country, raster, and color-coded district overlays.
    """
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
    """
    Helper function to visualize binarized satellite tiles using hvPlot.

    Args:
        arr: The continuous data array (e.g., radiance).
        satellite: Name of the satellite sensor for the plot title.
        threshold: The value used to separate 'lit' from 'unlit' pixels.

    Returns:
        An interactive HoloViews image plot.
    """
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
    """
    Create a side-by-side interactive comparison of two xarray tiles.

    Args:
        arr_left: Data array for the left-hand plot.
        arr_right: Data array for the right-hand plot.
        title_left: Title for the left plot.
        title_right: Title for the right plot.
        clim_left: Color limit range for the left plot.
        clim_right: Color limit range for the right plot.
        cmap_left: Colormap for the left plot.
        cmap_right: Colormap for the right plot.
        colorbar: Whether to display the colorbar.
        extra_plot_kwargs: Dictionary of additional arguments passed to hvplot.

    Returns:
        A HoloViews Layout containing two synchronized image plots.
    """
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


def plot_lit_country_comparison_by_urbanisation(
    df,
    colors: tuple = ("#006666", "#B2D8D8"),
    labels: tuple = ("lit", "nonlit"),
    height: int = 4,
    aspect: float = 1.5,
):
    """
    Plot stacked histograms of pixel counts by urbanization degree across countries.

    Args:
        df: DataFrame containing 'country_code', 'smod_class', and 'bm_binary'.
        colors: Tuple of two hex colors for 'lit' and 'non-lit' bars.
        labels: Tuple of strings for the legend.
        height: Height of each facet.
        aspect: Aspect ratio of each facet.
    """
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
        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))


def plot_lit_proportion_comparison_by_urbanisation(
    df,
    col: str = "continent",
    colors: tuple = ("#383060", "#dad7eb"),
    labels: tuple = ("lit", "unlit"),
    height: int = 4,
    aspect: float = 1.5,
):
    """
    Plot the proportion of lit vs unlit pixels as a function of urbanization degree.

    Useful for comparing urbanization light profiles across continents or regions.

    Args:
        df: DataFrame with geographic and classification data.
        col: The column to facet the plots by (e.g., 'continent').
        colors: Tuple of two hex colors for proportions.
        labels: Legend labels.
        height: Facet height.
        aspect: Facet aspect ratio.
    """
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
        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))


def plot_world_coverage(country_codes: list[str], country_names: list[str]):
    """
    Highlight specific countries on a global map.

    Args:
        country_codes: List of ISO3166-1-Alpha-3 codes to highlight.
        country_names: List of country names to highlight.
    """
    from matplotlib.colors import ListedColormap

    world = gpd.read_file("https://datahub.io/core/geo-countries/r/countries.geojson")

    country_mask = world["name"].isin(country_names)
    code_mask = world["ISO3166-1-Alpha-3"].isin(country_codes)
    world["is_good"] = country_mask | code_mask

    cmap = ListedColormap(["#bdbdbd", "#377eb8"])

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
    """
    Plot boxen plots (letter-value plots) of radiance distribution by urbanization degree.

    Args:
        df: DataFrame containing radiance, smod_class, and continent.
        vertical_line: If provided, draws a vertical dashed line at this x-value (threshold).
        height: Plot height.
        aspect: Plot aspect ratio.
        show_labels: Whether to show axis labels.
        show_legend: Whether to show the continent legend.
    """
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
            ax.axvline(x=1, color="#F9CB40", linestyle="--", linewidth=3)


def plot_scatter(
    df: pd.DataFrame,
    countries_to_plot: list[str],
    x: str = "hdi",
    y: str = "f1",
    xlabel: str = "Human Development Index",
    ylabel: str = "f1",
    figsize: tuple = (12, 7),
    legend_loc: str = "upper left",
):
    """
    Plot a scatter plot of metrics (e.g., F1 vs HDI) with country annotations.

    Args:
        df: DataFrame with country metrics.
        countries_to_plot: List of country names to label with arrows.
        x: Column name for the x-axis.
        y: Column name for the y-axis.
        xlabel: Custom label for x-axis.
        ylabel: Custom label for y-axis.
        figsize: Matplotlib figure size.
        legend_loc: Location of the legend.
    """
    df["is_plotted"] = False
    df.loc[df["country"].isin(countries_to_plot), "is_plotted"] = True

    plt.figure(figsize=figsize)
    with sns.axes_style("whitegrid", {"grid.color": ".9", "grid.linestyle": ":"}):
        ax = sns.scatterplot(
            data=df,
            x=x,
            y=y,
            hue="continent",
            palette=PALETTE,
            size="pixel_count",
            sizes=(20, 500),
            alpha=0.8,
        )

        ax.set_xlabel(xlabel, fontsize=13, labelpad=20)
        ax.set_ylabel(ylabel, fontsize=13, labelpad=20)

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
            fontsize=12,
            title_fontsize=13,
        )

        for i, row in df.iterrows():
            if row.is_plotted:
                ax.annotate(
                    row["country"],
                    (row[x], row[y]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    arrowprops=dict(arrowstyle="->", lw=0.5, color="black"),
                    fontsize=10,
                )

        plt.tight_layout()
        plt.show()


def plot_scatter_bokeh(
    df: pd.DataFrame,
    x: str = "hdi",
    y: str = "f1",
    xlabel: str = "Human Development Index",
    ylabel: str = "f1",
    legend_loc: str = "top_left",
):
    """
    Create an interactive Bokeh scatter plot with hover tooltips for country metrics.

    Args:
        df: DataFrame with country metrics.
        x: X-axis column name.
        y: Y-axis column name.
        xlabel: Label for x-axis.
        ylabel: Label for y-axis.
        legend_loc: Location of the interactive legend.
    """
    continents = df["continent"].unique().tolist()
    palette = Set2[max(3, len(continents))]

    source = ColumnDataSource.from_df(df)

    fig = figure(x_axis_label=xlabel, y_axis_label=ylabel, width=880)
    fig.scatter(
        x=x,
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
    """
    Plot a linear regression between HDI and mean monthly radiance.

    Args:
        df: DataFrame containing 'hdi_2020' and 'black_marble_radiance_monthly'.
        color: Color for the regression line and points.
    """
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
