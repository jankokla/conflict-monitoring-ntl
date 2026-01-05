import datetime
import logging
import os
import re
import time
import warnings
from abc import ABC, abstractmethod
from pathlib import Path, PosixPath
from typing import Literal
from zoneinfo import ZoneInfo

import ee
import geemap
import geopandas as gpd
import pandas as pd
import rasterio
import requests
import rioxarray
import xarray as xr
from blackmarble import BlackMarble
from ee.featurecollection import FeatureCollection
from ee.image import Image
from pydantic import ConfigDict, validate_call
from rioxarray.merge import merge_arrays
from shapely import box
from tqdm import tqdm

from conflict_monitoring_ntl.logger import logger

logging.getLogger("pyogrio._io").setLevel(logging.WARNING)

warnings.filterwarnings(
    "ignore",
    message="Connection pool is full, discarding connection*",
    module="urllib3.connectionpool",
)

try:
    ee.Initialize(project="ee-jankokla-ntl")
except ee.ee_exception.EEException:
    ee.Authenticate()
    ee.Initialize(project="ee-jankokla-ntl")


class BaseRaster(ABC):
    """
    Abstract base class for handling satellite raster data operations.
    
    Provides core utilities for spatial intersection, local file discovery,
    and Earth Engine integration.
    """
    DATA_FOLDER = ""
    FILE_PATTERN = ""
    DATE_REGEX = ""

    @abstractmethod
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable,
    ) -> xr.Dataset | None:
        """
        Abstract method to retrieve raster data as an xarray Dataset.

        Args:
            gdf: Input geospatial features to define the area of interest.
            date_range: Single date or list of dates for data retrieval.
            variable: The specific variable or band to extract.
        """
        pass

    def _get_matching_tiles(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
    ) -> dict[datetime.date, list[PosixPath]]:
        """
        Find local raster filepaths overlapping features within a specific date range.

        Args:
            gdf: Input features to test for spatial intersection.
            date_range: Date or list of dates for filtering tiles.

        Returns:
            Dictionary mapping dates to lists of matching raster filepaths.

        Raises:
            ValueError: If date_range is empty or invalid.
        """
        if isinstance(date_range, datetime.date):
            date_range = [date_range]

        data_path = Path(__file__).parents[2] / "data" / self.DATA_FOLDER
        filepaths = list(data_path.rglob(self.FILE_PATTERN))

        matching_files = {}

        for filepath in filepaths:
            try:
                date_str = re.findall(self.DATE_REGEX, filepath.stem)[0]  # type: ignore
            except AttributeError:
                continue

            date = datetime.datetime.strptime(date_str, "%Y%m%d").date()

            if date in date_range:
                with rasterio.open(filepath) as src:
                    bounds = src.bounds
                    raster_poly = box(
                        bounds.left, bounds.bottom, bounds.right, bounds.top
                    )
                    if gdf.to_crs(src.crs).intersects(raster_poly).any():
                        if date not in matching_files:
                            matching_files[date] = []
                        matching_files[date].append(filepath)

        return matching_files

    @staticmethod
    def _get_layer_from_ee(product, clip: FeatureCollection) -> xr.Dataset:
        """
        Fetch a specific product layer from Google Earth Engine and convert to xarray.

        Args:
            product: Earth Engine asset ID or image string.
            clip: Feature collection used to spatially clip the image.

        Returns:
            Xarray Dataset containing the clipped image data.
        """
        image = Image(product)
        pop = image.clip(clip)

        return geemap.ee_to_xarray(
            pop,
            geometry=clip.geometry(),
            projection=pop.projection(),
            ee_initialize=False,
        )

    @staticmethod
    def _get_patch(gdf: gpd.GeoDataFrame, src: rasterio.DatasetReader) -> xr.DataArray:
        """
        Extract a raster patch matching specific geometry and CRS from a local file.

        Args:
            gdf: GeoDataFrame holding the clipping geometries.
            src: Open rasterio dataset reader to sample from.

        Returns:
            Xarray DataArray containing the clipped and loaded raster patch.
        """
        gdf = gdf.to_crs(src.crs)
        patch = (
            rioxarray.open_rasterio(src, masked=True)
            .rio.clip(gdf.geometry.values, gdf.crs, from_disk=True)
            .load()
        )

        assert isinstance(patch, xr.DataArray)

        return patch


class GHSLPopulation(BaseRaster):
    """Handler for Global Human Settlement Layer (GHSL) Population data."""
    PRODUCT = "JRC/GHSL/P2023A/GHS_POP/2020"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "population_count",
    ) -> xr.Dataset | None:
        """
        Retrieve GHSL population count data from Earth Engine.

        Args:
            gdf: Input geospatial features for clipping.
            date_range: Target date (expected to be a single date for this static product).
            variable: Band name to extract.

        Returns:
            Xarray Dataset with population data reprojected to EPSG:4326.
        """
        if isinstance(date_range, list):
            assert len(date_range) == 1
            date_range = date_range[0]

        ee_clip = geemap.gdf_to_ee(gdf)
        assert isinstance(ee_clip, FeatureCollection)

        ds = self._get_layer_from_ee(self.PRODUCT, ee_clip)

        ds = (
            ds[[variable]]
            .rename({"X": "x", "Y": "y"})
            .transpose("time", "y", "x")
            .rio.reproject("EPSG:4326")
        )
        ds = ds.rename({"x": "lon", "y": "lat", variable: "ghsl_population"})
        ds.attrs["crs"] = "EPSG:4326"

        return ds.squeeze("time", drop=True)


class GHSLUrbanization(BaseRaster):
    """Handler for GHSL Degree of Urbanization (SMOD) data."""
    PRODUCT = "JRC/GHSL/P2023A/GHS_SMOD_V2-0/2020"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "smod_code",
    ) -> xr.Dataset | None:
        """
        Retrieve GHSL urbanization settlement classification from Earth Engine.

        Args:
            gdf: Input geospatial features for clipping.
            date_range: Target date for filtering.
            variable: Settlement classification band name.

        Returns:
            Xarray Dataset with SMOD codes.
        """
        if isinstance(date_range, list):
            assert len(date_range) == 1
            date_range = date_range[0]

        ee_clip = geemap.gdf_to_ee(gdf)
        assert isinstance(ee_clip, FeatureCollection)

        ds = self._get_layer_from_ee(self.PRODUCT, ee_clip)

        ds = (
            ds[[variable]]
            .rename({"X": "x", "Y": "y"})
            .transpose("time", "y", "x")
            .rio.reproject("EPSG:4326")
        )
        ds = ds.rename({"x": "lon", "y": "lat", variable: "smod_code"})
        ds.attrs["crs"] = "EPSG:4326"

        return ds.squeeze("time", drop=True)


class GHSLSurface(BaseRaster):
    """Handler for GHSL Built-up Surface characteristics."""
    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "built_surface",
    ) -> xr.Dataset:
        """
        Retrieve GHSL built-up surface area from Earth Engine.

        Args:
            gdf: Input geospatial features for clipping.
            date_range: Target date for filtering.
            variable: Built surface band name.

        Returns:
            Xarray Dataset with built-up surface density.
        """
        if isinstance(date_range, list):
            assert len(date_range) == 1
            date_range = date_range[0]

        ee_clip = geemap.gdf_to_ee(gdf)
        assert isinstance(ee_clip, FeatureCollection)

        product = "JRC/GHSL/P2023A/GHS_BUILT_S/2020"

        ds = self._get_layer_from_ee(product, ee_clip)

        ds = (
            ds[[variable]]
            .rename({"X": "x", "Y": "y"})
            .transpose("time", "y", "x")
            .rio.reproject("EPSG:4326")
        )
        ds = ds.rename({"x": "lon", "y": "lat", variable: "ghsl_surface"})
        ds.attrs["crs"] = "EPSG:4326"

        return ds.squeeze("time", drop=True)


class BlackMarbleEE(BaseRaster):
    """Handler for NASA Black Marble (VNP46A2) data via Earth Engine."""
    PRODUCT_BASE = "NASA/VIIRS/002/VNP46A2"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "Gap_Filled_DNB_BRDF_Corrected_NTL",
    ) -> xr.Dataset:
        """
        Fetch Black Marble Nighttime Lights data from Earth Engine for a specific date.

        Args:
            gdf: Input geospatial features for clipping.
            date_range: Target date to fetch.
            variable: The specific NTL band to retrieve.

        Returns:
            Xarray Dataset containing radiance values.
        """
        if isinstance(date_range, list):
            assert len(date_range) == 1
            date_range = date_range[0]

        ee_clip = geemap.gdf_to_ee(gdf)
        assert isinstance(ee_clip, FeatureCollection)

        date_str = date_range.strftime("%Y_%m_%d")
        product = f"NASA/VIIRS/002/VNP46A2/{date_str}"

        ds = self._get_layer_from_ee(product, ee_clip)
        ds = ds[[variable]].rename({variable: "black_marble_radiance"})

        return ds.squeeze("time", drop=True)


class BlackMarblePy(BaseRaster):
    """Handler for NASA Black Marble data using the local blackmarble library."""
    
    FREQUENCY_TO_VARIABLE = {
        "daily": "Gap_Filled_DNB_BRDF-Corrected_NTL",
        "monthly": "NearNadir_Composite_Snow_Free",
        "annual": "NearNadir_Composite_Snow_Free",
    }

    FREQUENCY_TO_PRODUCT = {
        "daily": "VNP46A2",
        "monthly": "VNP46A3",
        "annual": "VNP46A4",
    }

    def __init__(
        self,
        frequency: Literal["daily", "monthly", "annual"] = "daily",
        drop_values_by_quality_flag: list[int] = [],
    ) -> None:
        """
        Initialize the BlackMarblePy loader.

        Args:
            frequency: Temporal resolution (daily, monthly, or annual).
            drop_values_by_quality_flag: List of quality flags to mask out.
        """
        super().__init__()
        self.frequency = frequency

        output_dir = Path(__file__).parents[2] / "data" / "black_marble"
        self.bm = BlackMarble(
            output_directory=output_dir,
            drop_values_by_quality_flag=drop_values_by_quality_flag,
        )

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "Gap_Filled_DNB_BRDF_Corrected_NTL",
    ) -> xr.Dataset:
        """
        Fetch Black Marble data using the NASA API/local cache.

        Args:
            gdf: Input geospatial features.
            date_range: Single date or list of dates.
            variable: Radiance variable name.

        Returns:
            Xarray Dataset containing nighttime light radiance.
        """
        variable = self.FREQUENCY_TO_VARIABLE[self.frequency]

        if not isinstance(date_range, list):
            date_range = [date_range]

        ds = self.bm.raster(
            gdf,
            product_id=self.FREQUENCY_TO_PRODUCT[self.frequency],  # type: ignore
            date_range=date_range,
            variable=variable,
        ).drop_attrs()

        ds.attrs["crs"] = "EPSG:4326"
        ds = ds.rio.write_crs("EPSG:4326")

        ds = ds[[variable]].rename(
            {
                variable: f"black_marble_radiance_{self.frequency}",
                "x": "lon",
                "y": "lat",
            }
        )

        if ds.sizes["time"] == 1:
            return ds.squeeze("time", drop=True)

        return ds


class EnMAP(BaseRaster):
    """EnMAP satellite hyperspectral raster data loader."""

    BANDS = [10, 26, 44]

    @property
    def DATA_FOLDER(self) -> str:
        """Directory name for EnMAP data."""
        return "enmap"

    @property
    def FILE_PATTERN(self) -> str:
        """Glob pattern for VNIR spectral images."""
        return "*SPECTRAL_IMAGE_VNIR.TIF"

    @property
    def DATE_REGEX(self) -> str:
        """Regex to extract ISO timestamp from filename."""
        return r"(\d{8})T\d{6}Z"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "RGB",
    ) -> xr.Dataset | None:
        """
        Extract EnMAP VNIR bands for features over a date range.

        Args:
            gdf: Input geospatial features.
            date_range: Single date or list of dates to match acquisition.
            variable: Name to assign to the resulting dataset variable.

        Returns:
            Xarray Dataset containing concatenated bands for each matched date.
        """
        matching_files = self._get_matching_tiles(gdf, date_range)

        data_arrays = []
        for date, file_list in matching_files.items():
            data_arrays_per_date = []

            for tif_file in file_list:
                with rasterio.open(tif_file) as src:
                    patch = self._get_patch(gdf, src)

                    bands_patch = patch.isel(band=self.BANDS)
                    bands_patch = bands_patch.rio.reproject("EPSG:4326")
                    data_arrays_per_date.append(bands_patch)

            merged = merge_arrays(data_arrays_per_date)
            merged["time"] = date

            data_arrays.append(merged)

        if not data_arrays:
            cls_name = self.__class__.__name__
            logger.warning(f"[{cls_name}] No images for given region and date range.")
            return None

        combined = (
            xr.concat(data_arrays, dim="time", combine_attrs="drop_conflicts")
            .to_dataset(name=variable, promote_attrs=True)
            .sortby("time")
        )

        return combined


class SDGSat(BaseRaster):
    """SDGSat satellite (Glimmer) raster data loader and radiance calculator."""

    BAND_MAPPING = {"PL": 1, "PH": 2, "RGB": 3}
    SDGSAT_COEF_MAPPING = {
        "PH": (0.00008757, 0.0000183897, 0.675),
        "PL": (0.00008832, 0.0000167808, 0.675),
    }

    @property
    def DATA_FOLDER(self) -> str:
        """Directory name for SDGSat data."""
        return "sdgsat"

    @property
    def FILE_PATTERN(self) -> str:
        """Glob pattern for SDGSat Level 1H files."""
        return "*_LH.tif"

    @property
    def DATE_REGEX(self) -> str:
        """Regex to extract date from filename."""
        return r"\d{8}"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: Literal["PL", "PH", "HDR"] = "PH",
    ) -> xr.Dataset:
        """
        Extract SDGSat bands and calculate physical radiance.

        Args:
            gdf: Input geospatial features.
            date_range: Single date or list of dates.
            variable: The sensor product to extract (PL or PH).

        Returns:
            Xarray Dataset containing both raw DN and calculated radiance.
        """
        matching_files = self._get_matching_tiles(gdf, date_range)

        data_arrays = []
        for date, file_list in matching_files.items():
            data_arrays_per_date = []

            for tif_file in file_list:
                with rasterio.open(tif_file) as src:
                    patch = self._get_patch(gdf, src)

                    band_patch = patch.sel(band=self.BAND_MAPPING[variable])
                    band_patch = band_patch.drop_vars("band")
                    band_patch.attrs["long_name"] = variable

                    data_arrays_per_date.append(band_patch)

            merged = merge_arrays(data_arrays_per_date)
            merged["time"] = date

            data_arrays.append(merged)

        if not data_arrays:
            cls_name = self.__class__.__name__
            logger.warning(f"[{cls_name}] No images for given region and date range.")
            return xr.Dataset()

        combined = (
            xr.concat(data_arrays, dim="time", combine_attrs="drop_conflicts")
            .to_dataset(name=variable, promote_attrs=True)
            .sortby("time")
        )

        ds = combined.rio.reproject("EPSG:4326")

        DN = ds[variable]
        gain, bias, bw = self.SDGSAT_COEF_MAPPING[variable]

        cond = (~DN.isnull()) & (DN != 0)
        L = DN * gain + bias

        ds["sdgsat_radiance"] = xr.where(cond, L * 1e5 * bw, DN)

        ds = ds.drop_vars("spatial_ref").rename({variable: "sdgsat_dn"}).drop_attrs()
        ds.attrs["crs"] = "EPSG:4326"
        ds = ds.rio.write_crs("EPSG:4326")

        return ds.rename({"y": "lat", "x": "lon"}).squeeze("time", drop=True)


class Landsat:
    """Interface for USGS Landsat 8/9 Machine-to-Machine (M2M) API."""
    SERVICE_URL = "https://m2m.cr.usgs.gov/api/api/json/stable"
    DATASET = "landsat_ot_c2_l2"

    def __init__(self) -> None:
        """Initialize connection and default filters for Landsat search."""
        self._api_key = self._prompt_ers_login()
        self._spatial_filter = None
        self._metadata_filter = {
            "filterType": "value",
            "filterId": "61af9273566bb9a8",
            "value": "8",
        }
        self._cloud_cover_filter = {"min": 0, "max": 20}

    @property
    def DATA_FOLDER(self) -> str:
        """Directory name for Landsat downloads."""
        return "landsat"

    @property
    def FILE_PATTERN(self) -> str:
        """File pattern for Landsat products."""
        return "*_LH.tif"

    @property
    def DATE_REGEX(self) -> str:
        """Regex for Landsat date extraction."""
        return r"\d{8}"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str | None = None,
        timezone: str = "UTC",
    ) -> xr.Dataset | None:
        """
        Search for, download, and process Landsat scenes.

        Args:
            gdf: Input geospatial features for search area.
            date_range: Dates to search.
            variable: Specific band or variable to process.
            timezone: Local timezone for evening pass filtering.

        Returns:
            Xarray Dataset containing processed Landsat data.
        """
        self.spatial_filter = gdf

        scenes_df = self._scene_search(date_range)
        scene_ids = self._get_scenes_between_time(scenes_df, timezone)

        if not scene_ids:
            cls_name = self.__class__.__name__
            logger.warning(f"[{cls_name}] No images for given region and date range.")
            return None

        products = self._get_available_products(scene_ids)

        self._download_products(products)

        return xr.Dataset()

    def _download_products(self, products: list[dict]):
        """
        Request and handle downloads from the M2M API.

        Args:
            products: List of product dictionaries containing entity and product IDs.
        """
        label = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        download_req_payload = {"downloads": products, "label": label}

        download_request_results = self._send_request(
            os.path.join(self.SERVICE_URL, "download-request"),
            download_req_payload,
            self._api_key,
        )

        assert isinstance(download_request_results, dict)

        if len(download_request_results["availableDownloads"]) > 0:
            for result in download_request_results["availableDownloads"]:
                self._download_file(result["url"])
        elif len(download_request_results["preparingDownloads"]) > 0:
            preparingDownloadIds = []

            for result in download_request_results["preparingDownloads"]:
                preparingDownloadIds.append(result["downloadId"])

            download_ret_payload = {"label": label}
            print("Retrieving download urls...\n")
            download_retrieve_results = self._send_request(
                os.path.join(self.SERVICE_URL, "download-retrieve"),
                download_ret_payload,
                self._api_key,
            )

            assert isinstance(download_retrieve_results, dict)

            if download_retrieve_results:
                for result in download_retrieve_results["available"]:
                    if result["downloadId"] in preparingDownloadIds:
                        preparingDownloadIds.remove(result["downloadId"])
                        self._download_file(result["url"])

                for result in download_retrieve_results["requested"]:
                    if result["downloadId"] in preparingDownloadIds:
                        preparingDownloadIds.remove(result["downloadId"])
                        self._download_file(result["url"])

            while len(preparingDownloadIds) > 0:
                print(
                    f"{len(preparingDownloadIds)} downloads are not available yet. "
                    "Waiting for 30s to retrieve again\n"
                )

                time.sleep(30)
                download_retrieve_results = self._send_request(
                    os.path.join(self.SERVICE_URL, "download-retrieve"),
                    download_ret_payload,
                    self._api_key,
                )

                assert isinstance(download_retrieve_results, dict)

                if download_retrieve_results:
                    for result in download_retrieve_results["available"]:
                        if result["downloadId"] in preparingDownloadIds:
                            preparingDownloadIds.remove(result["downloadId"])
                            self._download_file(result["url"])

    def _get_available_products(self, scene_ids: list[str]) -> list[dict]:
        """
        Identify downloadable product IDs for a list of scene entity IDs.

        Args:
            scene_ids: List of USGS entity IDs.

        Returns:
            List of dictionaries with entityId and productId.
        """
        download_payload = {"datasetName": self.DATASET, "entityIds": scene_ids}
        download_options = self._send_request(
            os.path.join(self.SERVICE_URL, "download-options"),
            download_payload,
            self._api_key,
        )

        assert isinstance(download_options, dict)

        df = pd.json_normalize(download_options)
        df = df[(df.available) & (df.downloadSystem != "folder")]

        return (
            df[["entityId", "id"]]
            .rename(columns={"id": "productId"})
            .to_dict("records")
        )

    def _scene_search(
        self, date_range: datetime.date | list[datetime.date]
    ) -> pd.DataFrame:
        """
        Perform a scene search via the M2M API.

        Args:
            date_range: Dates to search.

        Returns:
            DataFrame containing search results.
        """
        if isinstance(date_range, datetime.date):
            date_range = [date_range]

        payload = self._get_search_payload(min(date_range), max(date_range))

        scenes = self._send_request(
            url=os.path.join(self.SERVICE_URL, "scene-search"),
            data=payload,
            api_key=self._api_key,
        )

        assert isinstance(scenes, dict)

        return pd.json_normalize(scenes["results"])

    def _get_scenes_between_time(
        self,
        scenes_df: pd.DataFrame,
        timezone: str,
        start_time: datetime.time = datetime.time(21, 0, 0),
        end_time: datetime.time = datetime.time(3, 0, 0),
    ) -> list[str]:
        """
        Filter scenes to find those acquired during a specific evening window.

        Args:
            scenes_df: DataFrame of scenes from search.
            timezone: Local timezone to calculate evening window.
            start_time: Start of the evening window.
            end_time: End of the evening window (next day).

        Returns:
            List of entity IDs matching the time criteria.
        """
        def get_datetime(meta_list, field):
            time_pattern = "%Y-%m-%d %H:%M:%S"
            for d in meta_list:
                if d["fieldName"] == field:
                    dt_str = d["value"][:19]
                    dt = datetime.datetime.strptime(dt_str, time_pattern)
                    dt_utc = dt.replace(tzinfo=ZoneInfo("UTC"))
                    dt_local = dt_utc.astimezone(ZoneInfo(timezone))
                    return dt_local

        def get_evening_anchor(dt):
            if dt.time() < end_time:
                return (dt - datetime.timedelta(days=1)).date()
            else:
                return dt.date()

        def is_within_evening(dt):
            return start_time <= dt.time() or dt.time() < end_time

        def is_in_time_window(start, end):
            return (
                get_evening_anchor(start) == get_evening_anchor(end)
                and is_within_evening(start)
                and is_within_evening(end)
                and start <= end
            )

        scenes_df["start_time"] = scenes_df["metadata"].apply(
            lambda meta: get_datetime(meta, "Start Time")  # type: ignore
        )
        scenes_df["stop_time"] = scenes_df["metadata"].apply(
            lambda meta: get_datetime(meta, "Stop Time")  # type: ignore
        )

        filtered_df = scenes_df[
            scenes_df.apply(
                lambda row: is_in_time_window(row["start_time"], row["stop_time"]),
                axis=1,
            )
        ]

        return filtered_df.entityId.tolist()

    def _get_search_payload(
        self, start_date: datetime.date, end_date: datetime.date
    ) -> dict:
        """
        Construct the JSON payload for the M2M scene-search request.

        Args:
            start_date: Earliest acquisition date.
            end_date: Latest acquisition date.

        Returns:
            Dictionary payload for the API request.
        """
        return {
            "datasetName": "landsat_ot_c2_l2",
            "sceneFilter": {
                "metadataFilter": self.metadata_filter,
                "spatialFilter": self.spatial_filter,
                "acquisitionFilter": {
                    "start": start_date.strftime("%Y-%m-%d"),
                    "end": end_date.strftime("%Y-%m-%d"),
                },
                "cloudCoverFilter": self.cloud_cover_filter,
            },
        }

    @property
    def cloud_cover_filter(self) -> dict:
        """Get the current cloud cover filter range."""
        return self._cloud_cover_filter

    @cloud_cover_filter.setter
    def cloud_cover_filter(self, min: int, max: int) -> None:
        """Set the cloud cover filter range."""
        self._cloud_cover_filter = {"min": min, "max": max}

    @property
    def metadata_filter(self) -> dict:
        """Get current metadata filters."""
        if self._metadata_filter:
            return self._metadata_filter
        raise ValueError("You need to set `metadata_filter` first.")

    @metadata_filter.setter
    def metadata_filter(self, satellite: Literal["8", "9"]) -> None:
        self._metadata_filter["value"] = satellite

    @property
    def spatial_filter(self) -> dict:
        if self._spatial_filter:
            return self._spatial_filter
        raise ValueError("You need to set `spatial_filter` first.")

    @spatial_filter.setter
    def spatial_filter(self, gdf: gpd.GeoDataFrame) -> None:
        self._spatial_filter = {
            "filterType": "mbr",
            "lowerLeft": {
                "latitude": gdf.bounds.miny[0],
                "longitude": gdf.bounds.minx[0],
            },
            "upperRight": {
                "latitude": gdf.bounds.maxy[0],
                "longitude": gdf.bounds.maxx[0],
            },
        }

    def _prompt_ers_login(self) -> str | None:
        """
        Sends a POST request to the EarthExplorer service to obtain a login token.

        Args:
            self: The instance reference.

        Returns:
            The 'data' field from the JSON response containing login information.

        Raises:
            requests.HTTPError: If the HTTP request fails.
            KeyError: If the 'data' field is missing in the response.
            ValueError: If the response is not valid JSON.
        """
        response = requests.post(
            os.path.join(self.SERVICE_URL, "login-token"),
            json={
                "username": os.getenv("EARTHEXPLORER_USERNAME"),
                "token": os.getenv("EARTHEXPLORER_TOKEN"),
            },
        )
        response.raise_for_status()
        return response.json()["data"]

    @staticmethod
    def _send_request(
        url: str,
        data: dict,
        api_key: str | None = None,
    ) -> dict | str | None:
        """
        Sends a POST request to a specified M2M endpoint and returns the parsed
            'data' field from the JSON response.

        Parameters:
            url: The endpoint URL to which the request is sent.
            data: Dictionary payload to send as the JSON body of the request.
            api_key: Authentication token added to the 'X-Auth-Token' header.

        Returns:
            dict or None: The value of the 'data' field from the JSON response,
                or None if the field is absent.

        Raises:
            requests.HTTPError: If the HTTP response status code indicates an error.
            ValueError: If the response body is not valid JSON.
        """
        headers = {"X-Auth-Token": api_key} if api_key else {}

        response = requests.post(url, json=data, headers=headers)
        response.raise_for_status()
        output = response.json()

        return output.get("data")

    def _download_file(self, url: str):
        output_dir = Path(__file__).parents[2] / "data" / self.DATA_FOLDER
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()

            # Attempt to fetch filename from headers, fallback if needed
            disposition = response.headers.get("content-disposition")
            if disposition:
                filenames = re.findall('filename="?([^"]+)"?', disposition)
                filename = filenames[0] if filenames else url.split("/")[-1]
            else:
                filename = url.split("/")[-1]

            # ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            file_path = os.path.join(output_dir, filename)

            total = int(response.headers.get("content-length", 0))
            chunk_size = 1024

            with (
                open(file_path, "wb") as f,
                tqdm(total=total, unit="B", unit_scale=True, desc=filename) as bar,
            ):
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        bar.update(len(chunk))

        except Exception as e:
            print(f"\nFailed to download from {url}: {e}")


class ErsLoginError(Exception):
    pass
