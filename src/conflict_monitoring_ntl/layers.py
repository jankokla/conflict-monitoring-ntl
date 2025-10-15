import datetime

import ee
import geemap
import xarray as xr
from ee.featurecollection import FeatureCollection
from ee.filter import Filter
from ee.image import Image
from rasterio.enums import Resampling

from conflict_monitoring_ntl.satellites import SDGSat

try:
    ee.Initialize()
except ee.ee_exception.EEException:
    ee.Authenticate()
    ee.Initialize()


ADMIN_UNIT_DS = "FAO/GAUL/2015/level2"


def get_rasters_for_admin(admin_name: str, date: datetime.date) -> xr.Dataset:
    admin = _get_feature_collection(admin_name)
    gdf = geemap.ee_to_gdf(admin)

    ghs_xds = _get_layer_from_ee("JRC/GHSL/P2023A/GHS_POP/2020", admin)
    bm_xds = _get_layer_from_ee("NOAA/VIIRS/DNB/ANNUAL_V22/20240101", admin)

    sdgsat = SDGSat()
    sdgsat_xds = sdgsat.raster(gdf, [date], variable="PH")

    assert isinstance(sdgsat_xds, xr.Dataset)

    ghs_xds = _standardize_ghs_pop(ghs_xds, bm_xds)
    sdgsat_xds = _standardize_sdgsat(sdgsat_xds, bm_xds)
    bm_xds = _standardize_black_marble(bm_xds)

    return xr.merge([ghs_xds, sdgsat_xds, bm_xds])


def _get_layer_from_ee(product: str, clip: FeatureCollection) -> xr.Dataset:
    image = Image(product)
    pop = image.clip(clip)

    return geemap.ee_to_xarray(
        pop, geometry=clip.geometry(), projection=pop.projection()
    )


def _standardize_black_marble(xds: xr.Dataset) -> xr.Dataset:
    return (
        xds.squeeze("time", drop=True)
        .drop_attrs()[["average"]]
        .rename({"average": "black_marble"})
        .transpose("lat", "lon")
    )


def _standardize_ghs_pop(xds: xr.Dataset, ref_xds: xr.Dataset) -> xr.Dataset:
    xds = xds.rename({"X": "x", "Y": "y"}).transpose("time", "y", "x")

    xds = xds.rio.reproject_match(
        ref_xds.rename({"lon": "x", "lat": "y"}),  # since func can only read x-y
        resampling=Resampling.sum,
    )

    return (
        xds.squeeze("time", drop=True)
        .drop_attrs()
        .drop_vars("spatial_ref")
        .rename({"x": "lon", "y": "lat"})
    )


def _standardize_sdgsat(xds: xr.Dataset, ref_xds: xr.Dataset) -> xr.Dataset:

    bandwidth = 0.675  # in micrometers
    gain = 0.00008757
    bias = 0.0000183897

    DN = xds["PH"]

    # Only apply to non-NaN and nonzero values
    cond = (~DN.isnull()) & (DN != 0)
    L = DN * gain + bias

    xds["radiance"] = xr.where(cond, L * 1e5 * bandwidth, DN)

    xds = xds.rio.reproject_match(
        ref_xds.rename({"lon": "x", "lat": "y"}),  # since func can only read x-y
        resampling=Resampling.average,
    )

    return (
        xds.squeeze("time", drop=True)
        .drop_attrs()
        .drop_vars(["spatial_ref", "PH"])
        .rename({"x": "lon", "y": "lat", "radiance": "sdgsat"})
    )


def _get_feature_collection(
    value: str, feature: str = "ADM2_NAME"
) -> FeatureCollection:
    countries = FeatureCollection(ADMIN_UNIT_DS)
    return countries.filter(Filter.eq(feature, value))
