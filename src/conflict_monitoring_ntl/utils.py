from functools import reduce
from typing import Union

import geemap
import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry
import pycountry_convert as pc
import pyproj
import rasterio
import xarray as xr
from ee.featurecollection import FeatureCollection
from pygadm import Items
from rasterio.crs import CRS
from shapely import Geometry
from shapely.geometry import box
from shapely.ops import transform, unary_union
from sklearn.metrics import precision_score, recall_score
from tqdm import tqdm


def get_raster_polygon(*raster_paths) -> tuple[Geometry, CRS]:
    """
    Computes the union of bounding polygons from any number of raster files.

    Args:
        *raster_paths: Variable number of raster file paths to evaluate.

    Returns:
        A tuple containing the CRS of the first raster and a single merged 
        geometry representing the total spatial extent of all inputs.
    """
    boxes = []
    crs = None

    for idx, path in enumerate(raster_paths):
        with rasterio.open(path) as src:
            if idx == 0:
                crs = src.crs
            boxes.append(box(*src.bounds))

    return unary_union(boxes), crs


def binarize_xarray(
    arr: xr.DataArray, threshold: float, crs: str = "EPSG:4326"
) -> xr.DataArray:
    """
    Convert a continuous data array into a binary mask based on a threshold.

    Args:
        arr: The input data array to process.
        threshold: The value above which pixels are set to 1.
        crs: The coordinate reference system string to assign to the output.

    Returns:
        A binarized data array where values > threshold are 1, others are 0, 
        maintaining original NaN locations.
    """
    new_arr = xr.where(arr.notnull(), (arr > threshold).astype(int), np.nan)
    return new_arr.rio.write_crs(crs)


def get_non_nan_flat_array(arr: xr.DataArray) -> np.ndarray:
    """
    Flatten an xarray and remove all NaN values.

    Args:
        arr: The input data array.

    Returns:
        A 1D numpy array containing only valid numeric data.
    """
    np_arr = arr.to_numpy().flatten()
    return np_arr[~np.isnan(np_arr)]


def get_combined_mask(ds: xr.Dataset) -> xr.DataArray:
    """
    Generate a bitwise-AND mask across all variables in a dataset.

    Identifies spatial locations where all variables in the dataset contain 
    finite (non-NaN) values. Useful for finding the valid intersection 
    across multiple satellite layers.

    Args:
        ds: The input dataset containing multiple data variables.

    Returns:
        A boolean data array where True represents valid pixels in all layers.
    """
    masks = [np.isfinite(ds[var]) for var in ds.data_vars]
    combined = reduce(lambda x, y: x & y, masks)

    assert isinstance(combined, xr.DataArray)

    return combined


def get_raster_gdf(*raster_paths, crs: str | None = None) -> gpd.GeoDataFrame:
    """
    Create a GeoDataFrame representing the spatial extent of one or more rasters.

    Args:
        *raster_paths: Paths to the raster files.
        crs: Optional target CRS for the resulting GeoDataFrame.

    Returns:
        A GeoDataFrame containing a single geometry entry for the merged extent.
    """
    polygon, tif_crs = get_raster_polygon(*raster_paths)
    gdf = gpd.GeoDataFrame(geometry=[polygon], crs=tif_crs)
    return gdf.to_crs(crs) if crs else gdf


def reproject_gdf(
    gdf: gpd.GeoDataFrame,
    to_crs: Union[str, CRS],
    from_crs: Union[str, CRS] = "EPSG:4326",
) -> gpd.GeoDataFrame:
    """
    Reproject a GeoDataFrame and store the result in a new geometry column.

    Args:
        gdf: Input GeoDataFrame.
        to_crs: The target coordinate reference system.
        from_crs: The source coordinate reference system.

    Returns:
        The GeoDataFrame with an additional 'geometry_proj' column.
    """
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform
    gdf["geometry_proj"] = gdf["geometry"].apply(lambda geom: transform(project, geom))
    return gdf


def get_gdf_for_admin(admin_id: str, content_level: int = 2) -> gpd.GeoDataFrame:
    """
    Retrieve administrative boundaries as a GeoDataFrame using GADM.

    Args:
        admin_id: The GADM administrative ID (e.g., 'UKR').
        content_level: The administrative depth (0=country, 1=state, 2=county).

    Returns:
        A GeoDataFrame containing the requested administrative shapes.
    """
    county_gdf = Items(admin=admin_id, content_level=content_level)  # type: ignore
    new_gdf = gpd.GeoDataFrame(geometry=county_gdf.geometry)
    new_gdf.set_crs("EPSG:4326", inplace=True)
    return new_gdf


def get_ee_for_admin(admin_id: str, content_level: int = 2) -> FeatureCollection:
    """
    Retrieve administrative boundaries as a Google Earth Engine FeatureCollection.

    Args:
        admin_id: The GADM administrative ID.
        content_level: The administrative depth.

    Returns:
        An Earth Engine FeatureCollection for the specified region.
    """
    county_gdf = get_gdf_for_admin(admin_id, content_level)
    county_ee = geemap.gdf_to_ee(county_gdf)

    assert isinstance(county_ee, FeatureCollection)

    return county_ee


def get_precision_recall(
    arr: xr.DataArray, y_true: np.ndarray, thresholds: list[float]
) -> pd.DataFrame:
    """
    Calculate precision and recall metrics across a range of binarization thresholds.

    Args:
        arr: Predicted radiance or intensity data array.
        y_true: Flattened ground truth binary array.
        thresholds: List of numeric thresholds to evaluate.

    Returns:
        A DataFrame indexed by threshold with precision and recall scores.
    """
    results = []

    for threshold in tqdm(thresholds, desc="Computing precision/recall"):
        bm_binary = binarize_xarray(arr, threshold)
        y_pred = get_non_nan_flat_array(bm_binary)

        results.append(
            {
                "precision": precision_score(y_true, y_pred),
                "recall": recall_score(y_true, y_pred),
                "threshold": round(threshold, 2),
            }
        )

    return pd.DataFrame(results)


def country_to_continent(alpha3):
    """
    Convert an ISO Alpha-3 country code to its corresponding continent name.

    Args:
        alpha3: Three-letter ISO country code (e.g., 'USA', 'FRA').

    Returns:
        String containing the continent name, or None if the code is invalid.
    """
    try:
        country = pycountry.countries.get(alpha_3=alpha3)
        assert isinstance(country, pycountry.db.Country)
        continent_code = pc.country_alpha2_to_continent_code(country.alpha_2)
        continent_name = pc.convert_continent_code_to_continent_name(continent_code)
        return continent_name
    except KeyError:
        return None
