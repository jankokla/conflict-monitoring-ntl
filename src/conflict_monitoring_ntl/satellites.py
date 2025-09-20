import datetime
import re
from pathlib import Path, PosixPath
from typing import Literal, Optional

import geopandas as gpd
import rasterio
import rioxarray
import xarray as xr
from pydantic import ConfigDict, validate_call
from rasterio.windows import from_bounds
from rioxarray.merge import merge_arrays
from shapely import box


class SDGSat:
    """SDGSat satellite raster data loader and band extraction utility."""

    BAND_MAPPING = {"PL": 1, "PH": 2, "HDR": 3}

    @validate_call(config=ConfigDict(arbitrary_types_allowed=True))
    def raster(
        self,
        gdf: gpd.GeoDataFrame,
        date_range: datetime.date | list[datetime.date],
        variable: Optional[Literal["PL", "PH", "HDR"]] = None,
    ) -> xr.Dataset:
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

        if variable is None:
            variable = "HDR"

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

        combined = (
            xr.concat(data_arrays, dim="time", combine_attrs="drop_conflicts")
            .to_dataset(name=variable, promote_attrs=True)
            .sortby("time")
        )
        combined = combined.rio.write_crs("EPSG:4326")

        return combined

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
        # reproject shape geometry to SDGSAT CRS
        patch_bounds = gdf.to_crs(src.crs).total_bounds  # xmin, ymin, xmax, ymax

        # clip SDGSAT patch
        window = from_bounds(*patch_bounds, transform=src.transform)

        # extract offsets and dimensions (cast to int if needed)
        row_off = int(window.row_off)
        col_off = int(window.col_off)
        h = int(window.height)
        w = int(window.width)

        da = rioxarray.open_rasterio(src)

        patch = da.isel(
            y=slice(row_off, row_off + h), x=slice(col_off, col_off + w)
        ).load()

        assert isinstance(patch, xr.DataArray)

        return patch

    @staticmethod
    def _get_matching_tiles(
        gdf: gpd.GeoDataFrame, date_range: datetime.date | list[datetime.date]
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

        start_date = min(date_range)
        end_date = max(date_range)

        data_path = Path(__file__).parent.parent.parent / "data"
        filepaths = list(data_path.rglob("*_LH.tif"))

        matching_files = {}

        for filepath in filepaths:
            try:
                date_str = re.search(r"\d{8}", filepath.stem).group(0)  # type: ignore
            except AttributeError:
                continue

            date = datetime.datetime.strptime(date_str, "%Y%m%d").date()

            if start_date <= date <= end_date:
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
