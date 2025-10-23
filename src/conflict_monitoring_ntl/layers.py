import datetime
from typing import Literal

import ee
import geemap
import geopandas as gpd
import numpy as np
import pygadm
import xarray as xr
from ee.featurecollection import FeatureCollection
from ee.image import Image
from rasterio.enums import Resampling

from conflict_monitoring_ntl.satellites import SDGSat

try:
    ee.Initialize()
except ee.ee_exception.EEException:
    ee.Authenticate()
    ee.Initialize()


# gain, bias, bandwidth (in micrometers)
SDGSAT_COEF_MAPPING = {
    "PH": (0.00008757, 0.0000183897, 0.675),
    "PL": (0.00008832, 0.0000167808, 0.675),
}


def get_rasters_for_admin(
    admin_id: str,
    date: datetime.date,
    is_sdgsat: bool = True,
    is_ghspop: bool = False,
    sdgsat_dn: Literal["PL", "PH", "HDR"] = "PL",
    **kwargs,
) -> xr.Dataset:
    # TODO: add docstring
    county_gdf, county_ee = _get_gdf_ee_for_admin(admin_id)

    date_str = date.strftime("%Y_%m_%d")
    bm_xds = _get_layer_from_ee(f"NASA/VIIRS/002/VNP46A2/{date_str}", county_ee)

    arrays = []

    if is_ghspop:
        ghs_xds = _get_layer_from_ee("JRC/GHSL/P2023A/GHS_POP/2020", county_ee)
        ghs_xds = _standardize_ghs_pop(ghs_xds, bm_xds)
        arrays.append(ghs_xds)

    if is_sdgsat:
        sdgsat = SDGSat()
        sdgsat_xds = sdgsat.raster(county_gdf, [date], variable=sdgsat_dn)

        assert isinstance(sdgsat_xds, xr.Dataset)
        sdgsat_xds = _standardize_sdgsat(
            sdgsat_xds, bm_xds, sdgsat_dn=sdgsat_dn, **kwargs
        )
        arrays.append(sdgsat_xds)

    bm_xds = _standardize_black_marble(bm_xds, **kwargs)
    arrays.append(bm_xds)

    return xr.merge(arrays).rio.write_crs("EPSG:4326")


def _get_layer_from_ee(product: str, clip: FeatureCollection) -> xr.Dataset:
    image = Image(product)
    pop = image.clip(clip)

    return geemap.ee_to_xarray(
        pop, geometry=clip.geometry(), projection=pop.projection()
    )


def _standardize_black_marble(
    xds: xr.Dataset,
    bm_is_binary: bool = True,
    bm_binary_thresholds: list[float] = [0.5, 1.0],
    **kwargs,
) -> xr.Dataset:
    # TODO: add docstring
    band = "Gap_Filled_DNB_BRDF_Corrected_NTL"
    xds = xds[[band]]

    if bm_is_binary:
        DN = xds[band]
        for thresh in bm_binary_thresholds:
            threshold_str = str(thresh).replace(".", "_")
            xds[f"black_marble_binary_{threshold_str}"] = xr.where(
                DN.notnull(), (DN > thresh).astype(int), np.nan
            )

    return (
        xds.squeeze("time", drop=True)
        .rename({band: "black_marble_radiance"})
        .transpose("lat", "lon")
    )


def _standardize_ghs_pop(xds: xr.Dataset, ref_xds: xr.Dataset) -> xr.Dataset:
    # TODO: add docstring
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


def _standardize_sdgsat(
    xds: xr.Dataset,
    ref_xds: xr.Dataset,
    sdgsat_dn: str = "PL",
    is_radiance: bool = True,
    is_binary: bool = True,
    sdgsat_binary_thresholds: list[float] = [4.0],
    **kwargs,
) -> xr.Dataset:
    # TODO: add docstring
    DN = xds[sdgsat_dn]

    if is_radiance:

        gain, bias, bw = SDGSAT_COEF_MAPPING[sdgsat_dn]

        # only apply to non-NaN and nonzero values
        cond = (~DN.isnull()) & (DN != 0)
        L = DN * gain + bias

        xds["sdgsat_radiance"] = xr.where(cond, L * 1e5 * bw, DN)

    xds = xds.rio.reproject_match(
        ref_xds.rename({"lon": "x", "lat": "y"}),  # since func can only read x-y
        resampling=Resampling.average,
    )

    if is_binary:
        DN = xds[sdgsat_dn]
        for thresh in sdgsat_binary_thresholds:
            threshold_str = str(thresh).replace(".", "_")
            xds[f"sdgsat_binary_{threshold_str}"] = xr.where(
                DN.notnull(), (DN > thresh).astype(int), np.nan
            )

    return (
        xds.squeeze("time", drop=True)
        .drop_attrs()
        .drop_vars(["spatial_ref"])
        .rename({"x": "lon", "y": "lat", sdgsat_dn: "sdgsat_dn"})
    )


def _get_gdf_ee_for_admin(
    admin_id: str, content_level: int = 2
) -> tuple[gpd.GeoDataFrame, FeatureCollection]:
    # TODO: add docstring
    county_gdf = pygadm.Items(admin=admin_id, content_level=content_level)  # type: ignore
    county_gdf.set_crs("EPSG:4326", inplace=True)
    county_ee = geemap.gdf_to_ee(county_gdf)

    assert isinstance(county_ee, FeatureCollection)

    return county_gdf, county_ee
