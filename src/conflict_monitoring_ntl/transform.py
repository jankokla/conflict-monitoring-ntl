import datetime

import geopandas as gpd
import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from conflict_monitoring_ntl.satellites import BaseRaster


def binarization(
    ds: xr.Dataset, satellite: str, threshold_var: str, thresholds: list[float]
):
    """
    Create binarized thresholded variables for a dataset.

    Takes a continuous raster variable and converts it into a series of binary 
    (0/1) masks based on a list of provided thresholds. Values greater than 
    the threshold are set to 1, while others are set to 0 (preserving NaNs).

    Args:
        ds: The input dataset to augment with binarized variables.
        satellite: Prefix for output variable names; typically the satellite name.
        threshold_var: The variable name in the dataset to apply the threshold to.
        thresholds: List of thresholds to apply for binarization.

    Returns:
        The dataset with new binarized variables added for each threshold level.
    """
    t_var = ds[threshold_var]
    satellite_str = satellite.lower().replace(" ", "_")
    for thresh in thresholds:
        d_var = f"{satellite_str}_binary_{str(thresh).replace('.', '_')}"
        ds[d_var] = xr.where(t_var.notnull(), (t_var > thresh).astype(int), np.nan)

    return ds.rio.write_crs(ds.rio.crs)


def reproject_match(ds: xr.Dataset, ref_ds: xr.Dataset, resampling: Resampling):
    """
    Reproject and match an xarray.Dataset to a reference dataset.

    Aligns the input dataset's grid, CRS, and resolution to match the 
    reference dataset. This is a critical step before merging datasets 
    from different satellite sensors.

    Args:
        ds: The input dataset to be re-projected and matched.
        ref_ds: The reference dataset providing the target grid and CRS.
        resampling: The resampling method to use (e.g., nearest, bilinear).

    Returns:
        The re-projected dataset with dimensions standardized to 'lon' and 'lat'.
    """

    def rename_if_necessary(ds: xr.Dataset):
        """
        Standardize coordinate names for rioxarray processing.

        Rioxarray's reproject_match logic often expects 'x' and 'y' naming 
        conventions rather than 'lon' and 'lat'.
        """
        dims = list(ds.sizes)
        if "lon" in dims and "lat" in dims:
            ds = ds.rename({"lon": "x", "lat": "y"})
        return ds

    ds = rename_if_necessary(ds)
    ref_ds = rename_if_necessary(ref_ds)

    ds = ds.transpose("y", "x", missing_dims="ignore")
    ref_ds = ref_ds.transpose("y", "x", missing_dims="ignore")

    ds = ds.rio.reproject_match(ref_ds, resampling=resampling)

    return ds.rename({"x": "lon", "y": "lat"})


TRANSFORMS = {"binarization": binarization, "reproject_match": reproject_match}


class RasterPipeline:
    """
    Orchestrates the loading and transformation of multiple satellite rasters.

    This class manages a sequence of satellite data objects, fetches their 
    respective rasters for a specific date and area, and applies a pipeline 
    of transformations (like binarization and alignment) to each.
    """
    def __init__(
        self,
        shape: gpd.GeoDataFrame,
        date: datetime.date,
        rasters: list[BaseRaster],
        configs: list[list[dict]],
        base_index: int = -1,
    ):
        """
        Initialize the RasterPipeline.

        Args:
            shape: GeoDataFrame defining the spatial area of interest.
            date: The target date for data extraction.
            rasters: List of BaseRaster implementations for different sensors.
            configs: Nested list of dictionaries defining transformations per raster.
            base_index: Index of the raster to use as the spatial reference for alignment.
        """
        self.shape = shape
        self.date = date
        self.rasters = rasters
        self.xarrays = []
        self.configs = configs
        self.base_index = base_index

    def run(self):
        """
        Process all rasters, apply pipelines, and merge results.

        Executes the full workflow: downloading/loading rasters, applying 
        user-defined transformations, and merging all results into a 
        single multi-variable dataset.

        Returns:
            A consolidated xarray Dataset containing all processed satellite data.
        """
        processed = []

        self.xarrays = [r.raster(self.shape, self.date) for r in self.rasters]  # type: ignore

        for ds, cfg in zip(self.xarrays, self.configs):
            ds = self.apply_transformations(ds, cfg)  # type: ignore
            processed.append(ds)

        return xr.merge(processed).rio.write_crs("EPSG:4326")

    def apply_transformations(self, ds, pipeline):
        """
        Apply a series of transformation functions to a single dataset.

        Args:
            ds: The input xarray Dataset to transform.
            pipeline: A dictionary mapping transformation names to their kwargs.

        Returns:
            The transformed xarray Dataset.
        """
        for func_str, kwargs in pipeline.items():
            fn = TRANSFORMS[func_str]

            if func_str == "reproject_match":
                ds = fn(ds, ref_ds=self.xarrays[self.base_index], **kwargs)
            else:
                ds = fn(ds, **kwargs)

        return ds
