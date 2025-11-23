import datetime
import os
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry
import pygadm
from rasterio.enums import Resampling
from tqdm import tqdm

from conflict_monitoring_ntl.satellites import (
    BlackMarblePy,
    GHSLPopulation,
    GHSLSurface,
    GHSLUrbanization,
)
from conflict_monitoring_ntl.transform import RasterPipeline
from conflict_monitoring_ntl.utils import country_to_continent, get_combined_mask

checkpoint_path = "../results/countries_done.parquet"
removed = ["AUS", "BRA", "CAN", "USA", "RUS", "GRL", "MEX", "CHL", "IDN"]

if os.path.exists(checkpoint_path):
    df = pd.read_parquet(checkpoint_path)
    completed = set(df.gid.to_list())
else:
    df = pd.DataFrame(columns=["country", "gid"])
    completed = set()

completed.update(removed)

error_path = "../results/country_errors.parquet"

if os.path.exists(error_path):
    error_df = pd.read_parquet(error_path)
    error = set(error_df.gid.to_list())
    completed.update(error)
else:
    error_df = pd.DataFrame(columns=["country", "gid", "error"])

rasters = [
    GHSLPopulation(),
    GHSLSurface(),
    BlackMarblePy(frequency="monthly"),
    GHSLUrbanization(),
]
transformations = [
    {"reproject_match": {"resampling": Resampling.sum}},
    {"reproject_match": {"resampling": Resampling.sum}},
    {"reproject_match": {"resampling": Resampling.average}},
    {},
]
date = datetime.date(2020, 1, 1)

with tqdm(pycountry.countries, desc="Retriving Tiles:") as pbar:
    for country_dto in pbar:
        country = country_dto.name
        gid = country_dto.alpha_3

        pbar.set_postfix(country=country)

        if gid in completed:
            continue

        try:
            gdf = pygadm.Items(admin=gid, content_level=1)

            pixels = 0
            conf_mat = np.zeros(4, dtype=np.int64)

            for i in tqdm(range(len(gdf)), desc="Processing Provinces"):
                original_gdf = gdf.iloc[[i]]
                province_gdf = gpd.GeoDataFrame(original_gdf.geometry).set_crs(
                    "EPSG:4326"
                )

                pipeline = RasterPipeline(province_gdf, date, rasters, transformations)
                ds = pipeline.run()

                # make sure we compare non-nan areas
                mask = get_combined_mask(ds)
                ds = ds.where(mask)

                province_df = (
                    ds.to_dataframe().reset_index().drop(columns=["spatial_ref"])
                )
                cols_non_nan = [
                    "ghsl_population",
                    "ghsl_surface",
                    "black_marble_radiance_monthly",
                    "smod_code",
                ]
                province_df = province_df.dropna(subset=cols_non_nan)

                province_df = province_df.astype(
                    {
                        "lat": np.float32,
                        "lon": np.float32,
                        "ghsl_population": np.float32,
                        "ghsl_surface": np.float32,
                        "black_marble_radiance_monthly": np.float32,
                        "smod_code": np.int8,
                    }
                )

                province_df["country_code"] = original_gdf.GID_0.item()
                province_df["admin_1_code"] = (
                    original_gdf.GID_1.item() if "GID_1" in original_gdf else np.nan
                )
                province_df["continent"] = country_to_continent(
                    original_gdf.GID_0.item()
                )
                province_df.to_parquet(
                    path="../results/urbanization",
                    partition_cols=["continent", "country_code"],
                )

            data = [country, gid]
            df = pd.concat(
                [pd.DataFrame([data], columns=df.columns), df], ignore_index=True
            )
            df.to_parquet(checkpoint_path)

            bm_path = Path(os.path.abspath("")).parent / "data" / "black_marble"
            shutil.rmtree(bm_path)
            bm_path.mkdir(exist_ok=True)

        except Exception as e:
            data = [country, gid, str(e)]
            error_df = pd.concat(
                [pd.DataFrame([data], columns=error_df.columns), error_df],
                ignore_index=True,
            )
            error_df.to_parquet(error_path)
