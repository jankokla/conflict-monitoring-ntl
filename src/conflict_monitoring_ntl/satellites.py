import datetime
import os
import re
import time
from abc import ABC, abstractmethod
from pathlib import Path, PosixPath
from typing import Literal
from zoneinfo import ZoneInfo

import geopandas as gpd
import pandas as pd
import rasterio
import requests
import rioxarray
import xarray as xr
from pydantic import ConfigDict, validate_call
from rioxarray.merge import merge_arrays
from shapely import box
from tqdm import tqdm

from conflict_monitoring_ntl.logger import logger


class Satellite(ABC):

    @property
    @abstractmethod
    def DATA_FOLDER(self) -> str:
        pass

    @property
    @abstractmethod
    def FILE_PATTERN(self) -> str:
        pass

    @property
    @abstractmethod
    def DATE_REGEX(self) -> str:
        pass

    @abstractmethod
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable,
    ) -> xr.Dataset | None:
        pass

    def _get_matching_tiles(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
    ) -> dict[datetime.date, list[PosixPath]]:
        """
        Find raster filepaths overlapping features within date range.

        Args:
            gdf: Input features to test for intersection.
            date_range: Date or list of dates for filtering tiles.

        Returns:
            Dictionary mapping dates to lists of matching raster filepaths.

        Raises:
            ValueError: If 'date_range' is empty or invalid.

        Example:
            >>> files = sdgsat.get_matching_tiles(my_gdf, [dt1, dt2])
            >>> print(files)
            {datetime.date(2023, 1, 5): [PosixPath("data/20230105_LH.tif")]}
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
    def _get_patch(gdf: gpd.GeoDataFrame, src: rasterio.DatasetReader) -> xr.DataArray:
        """
        Extract raster patch matching geometry and CRS.

        Args:
            gdf: GeoDataFrame holding clipping geometries.
            src: Open rasterio dataset to sample.

        Returns:
            xarray.DataArray containing area of interest raster patch.
        """
        gdf = gdf.to_crs(src.crs)
        patch = (
            rioxarray.open_rasterio(src, masked=True)
            .rio.clip(gdf.geometry.values, gdf.crs, from_disk=True)
            .load()
        )

        assert isinstance(patch, xr.DataArray)

        return patch


class EnMAP(Satellite):
    """EnMAP satellite raster data loader and band extraction utility."""

    BANDS = [10, 26, 44]

    @property
    def DATA_FOLDER(self) -> str:
        return "enmap"

    @property
    def FILE_PATTERN(self) -> str:
        return "*SPECTRAL_IMAGE_VNIR.TIF"

    @property
    def DATE_REGEX(self) -> str:
        return r"(\d{8})T\d{6}Z"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str = "RGB",
    ) -> xr.Dataset | None:
        """
        Extract raster data for features over a date range.

        Args:
            gdf: Input geospatial features.
            date_range: Single date or list of dates to match tile acquisition.
            variable: Optional band name to extract ("PL", "PH", or "HDR").

        Returns:
            xarray.Dataset containing selected data for each date.

        Raises:
            ValueError: If no matching tiles are found.
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


class SDGSat(Satellite):
    """SDGSat satellite raster data loader and band extraction utility."""

    BAND_MAPPING = {"PL": 1, "PH": 2, "RGB": 3}
    GAIN_BIAS_MAPPING = {"PL": (8.832, 1.67808), "PH": (8.758, 1.83897)}

    @property
    def DATA_FOLDER(self) -> str:
        return "sdgsat"

    @property
    def FILE_PATTERN(self) -> str:
        return "*_LH.tif"

    @property
    def DATE_REGEX(self) -> str:
        return r"\d{8}"

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: Literal["PL", "PH", "HDR"] = "HDR",
    ) -> xr.Dataset | None:
        """
        Extract raster data for features over a date range.

        Args:
            gdf: Input geospatial features.
            date_range: Single date or list of dates to match tile acquisition.
            variable: Optional band name to extract ("PL", "PH", or "HDR").

        Returns:
            xarray.Dataset containing selected data for each date.

        Raises:
            ValueError: If no matching tiles are found.
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
            return None

        combined = (
            xr.concat(data_arrays, dim="time", combine_attrs="drop_conflicts")
            .to_dataset(name=variable, promote_attrs=True)
            .sortby("time")
        )
        combined = combined.rio.reproject("EPSG:4326")

        return combined


class Landsat:

    SERVICE_URL = "https://m2m.cr.usgs.gov/api/api/json/stable"
    DATASET = "landsat_ot_c2_l2"

    def __init__(self) -> None:
        self._api_key = self._prompt_ers_login()
        self._spatial_filter = None
        self._metadata_filter = {
            "filterType": "value",
            "filterId": "61af9273566bb9a8",
            "value": "8",  # landsat-8
        }
        self._cloud_cover_filter = {"min": 0, "max": 20}

    @property
    def DATA_FOLDER(self) -> str:
        return "landsat"

    @property
    def FILE_PATTERN(self) -> str:
        return "*_LH.tif"  # TODO: update

    @property
    def DATE_REGEX(self) -> str:
        return r"\d{8}"  # TODO: update

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: str | None = None,  # TODO: update
        timezone: str = "UTC",
    ) -> xr.Dataset | None:

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
            # Retrieve download URLs
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

            # didn't get all download URLs, retrieve again after 30 seconds
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
        download_payload = {"datasetName": self.DATASET, "entityIds": scene_ids}
        download_options = self._send_request(
            os.path.join(self.SERVICE_URL, "download-options"),
            download_payload,
            self._api_key,
        )

        assert isinstance(download_options, dict)

        # filter out only necessary products
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

        def get_datetime(meta_list, field):
            time_pattern = "%Y-%m-%d %H:%M:%S"
            for d in meta_list:
                if d["fieldName"] == field:
                    # Remove microseconds if present
                    dt_str = d["value"][:19]
                    dt = datetime.datetime.strptime(dt_str, time_pattern)
                    dt_utc = dt.replace(tzinfo=ZoneInfo("UTC"))
                    dt_local = dt_utc.astimezone(ZoneInfo(timezone))
                    return dt_local

        def get_evening_anchor(dt):
            # shifts anchor to 'date of evening' (not after midnight)
            if dt.time() < end_time:
                return (dt - datetime.timedelta(days=1)).date()
            else:
                return dt.date()

        def is_within_evening(dt):
            # time must be between 21:00 - 23:59:59, OR 00:00 - 02:59:59
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
        return self._cloud_cover_filter

    @cloud_cover_filter.setter
    def cloud_cover_filter(self, min: int, max: int) -> None:
        self._cloud_cover_filter = {"min": min, "max": max}

    @property
    def metadata_filter(self) -> dict:
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

            # TODO: untar the file
        except Exception as e:
            print(f"\nFailed to download from {url}: {e}")


class ErsLoginError(Exception):
    pass
