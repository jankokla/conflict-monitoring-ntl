from functools import reduce
from typing import Union

import geemap
import geopandas as gpd
import numpy as np
import pandas as pd
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
        *raster_paths: Variable number of raster file paths.

    Returns:
        Tuple containing CRS (from the first raster) and the unary union
            geometry of all bounding boxes.
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
    new_arr = xr.where(arr.notnull(), (arr > threshold).astype(int), np.nan)
    return new_arr.rio.write_crs(crs)


def get_non_nan_flat_array(arr: xr.DataArray) -> np.ndarray:
    np_arr = arr.to_numpy().flatten()
    return np_arr[~np.isnan(np_arr)]


def get_combined_mask(ds: xr.Dataset) -> xr.DataArray:
    """TODO: update"""
    masks = [np.isfinite(ds[var]) for var in ds.data_vars]
    combined = reduce(lambda x, y: x & y, masks)

    assert isinstance(combined, xr.DataArray)

    return combined


def get_raster_gdf(*raster_paths, crs: str | None = None) -> gpd.GeoDataFrame:
    polygon, tif_crs = get_raster_polygon(*raster_paths)
    gdf = gpd.GeoDataFrame(geometry=[polygon], crs=tif_crs)
    return gdf.to_crs(crs) if crs else gdf


def reproject_gdf(
    gdf: gpd.GeoDataFrame,
    to_crs: Union[str, CRS],
    from_crs: Union[str, CRS] = "EPSG:4326",
) -> gpd.GeoDataFrame:
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform
    gdf["geometry_proj"] = gdf["geometry"].apply(lambda geom: transform(project, geom))
    return gdf


def get_gdf_for_admin(admin_id: str, content_level: int = 2) -> gpd.GeoDataFrame:
    # TODO: add docstring
    county_gdf = Items(admin=admin_id, content_level=content_level)  # type: ignore
    county_gdf.set_crs("EPSG:4326", inplace=True)
    return county_gdf


def get_ee_for_admin(admin_id: str, content_level: int = 2) -> FeatureCollection:
    # TODO: add docstring
    county_gdf = get_gdf_for_admin(admin_id, content_level)
    county_ee = geemap.gdf_to_ee(county_gdf)

    assert isinstance(county_ee, FeatureCollection)

    return county_ee


def get_precision_recall(
    arr: xr.DataArray, y_true: np.ndarray, thresholds: list[float]
) -> pd.DataFrame:
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
