import re
from datetime import datetime
from pathlib import Path, PosixPath

import rasterio
from geopandas import GeoDataFrame
from shapely import box


class SDGSat:

    def get_matching_tiles(
        self, gdf: GeoDataFrame, date_range: datetime | list[datetime]
    ) -> list[PosixPath]:
        """
        Find all raster files overlapping the given geospatial features and
        within the specified date range.

        Args:
            gdf (GeoDataFrame): Geospatial features to match against raster footprints.
            date_range (datetime or list[datetime]): A single date or a list specifying
                the date interval (inclusive). Raster tiles will be matched to files
                whose acquisition date falls within these bounds.

        Returns:
            list[PosixPath]: List of filepaths to rasters that spatially intersect the
                input features and are within the date range.

        Raises:
            ValueError: If 'date_range' is empty.

        Example:
            >>> matches = sdgsat.get_matching_tiles(my_gdf, [datetime(2023, 1, 1), datetime(2023, 1, 31)]). # noqa: E501
            >>> print(matches)
            [PosixPath('data/20230105_LH.tif'), PosixPath('data/20230129_LH.tif')]
        """
        if isinstance(date_range, datetime):
            date_range = [date_range]

        start_date = min(date_range)
        end_date = max(date_range)

        data_path = Path(__file__).parent.parent.parent / "data"
        filepaths = list(data_path.rglob("*_LH.tif"))

        matching_files = []

        for filepath in filepaths:
            try:
                date_str = re.search(r"\d{8}", filepath.stem).group(0)  # type: ignore
            except AttributeError:
                continue

            date = datetime.strptime(date_str, "%Y%m%d")

            if start_date <= date <= end_date:
                with rasterio.open(filepath) as src:
                    bounds = src.bounds
                    raster_poly = box(
                        bounds.left, bounds.bottom, bounds.right, bounds.top
                    )
                    if gdf.to_crs(src.crs).intersects(raster_poly).any():
                        matching_files.append(filepath)

        return matching_files
