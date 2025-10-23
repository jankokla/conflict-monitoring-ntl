from typing import Union

import geopandas as gpd
import pyproj
import rasterio
from rasterio.crs import CRS
from shapely import Geometry
from shapely.geometry import box
from shapely.ops import transform, unary_union


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
