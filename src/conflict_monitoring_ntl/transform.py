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

    Args:
        ds: The input dataset to augment with binarized variables.
        satellite: Prefix for output variable names; typically the satellite name.
        threshold_var: The variable name in `ds` to threshold.
        thresholds: List of thresholds to apply for binarization.

    Returns:
        The dataset with new binarized variables added for each threshold.
    """
    t_var = ds[threshold_var]
    satellite_str = satellite.lower().replace(" ", "_")
    for thresh in thresholds:
        d_var = f"{satellite_str}_binary_{str(thresh).replace(".", "_")}"
        ds[d_var] = xr.where(t_var.notnull(), (t_var > thresh).astype(int), np.nan)

    return ds.rio.write_crs(ds.rio.crs)  # write crs to all data vars


def reproject_match(ds: xr.Dataset, ref_ds: xr.Dataset, resampling: Resampling):
    """
    Reproject and match an xarray.Dataset to the grid, CRS, and resolution of a
        reference Dataset.

    Args:
        ds: The input dataset to be re-projected and matched.
        ref_ds: The reference dataset to match grid, CRS, and resolution.
        resampling: The resampling method to use during re-projection.

    Returns:
        The re-projected and rematched dataset, with dimensions renamed to
            'lon' and 'lat'.
    """

    def rename_if_necessary(ds: xr.Dataset):
        "Since reproject assumes 'x' and 'y' instead of 'lon' and 'lat'."
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
    def __init__(
        self,
        shape: gpd.GeoDataFrame,
        date: datetime.date,
        rasters: list[BaseRaster],
        configs: list[list[dict]],
        base_index: int = -1,
    ):
        self.shape = shape
        self.date = date
        self.rasters = rasters
        self.xarrays = []
        self.configs = configs
        self.base_index = base_index

    def run(self):
        """Process all rasters, apply pipelines, reproject and concatenate."""
        processed = []

        self.xarrays = [r.raster(self.shape, self.date) for r in self.rasters]  # type: ignore

        for ds, cfg in zip(self.xarrays, self.configs):
            ds = self.apply_transformations(ds, cfg)  # type: ignore
            processed.append(ds)

        return xr.merge(processed).rio.write_crs("EPSG:4326")

    def apply_transformations(self, ds, pipeline):
        """Apply a list of transformations to an xarray object."""
        for func_str, kwargs in pipeline.items():
            fn = TRANSFORMS[func_str]

            if func_str == "reproject_match":
                ds = fn(ds, ref_ds=self.xarrays[self.base_index], **kwargs)
            else:
                ds = fn(ds, **kwargs)

        return ds
