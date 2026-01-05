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
from sklearn.metrics import confusion_matrix
from tqdm import tqdm

from conflict_monitoring_ntl.satellites import BlackMarblePy, GHSLSurface
from conflict_monitoring_ntl.transform import RasterPipeline
from conflict_monitoring_ntl.utils import (
    binarize_xarray,
    get_combined_mask,
    get_non_nan_flat_array,
)

checkpoint_path = "../results/country_confusion_matrix.parquet"
removed = ["AUS", "BRA", "CAN", "USA", "RUS", "GRL", "MEX", "CHL", "IDN"]

if os.path.exists(checkpoint_path):
    df = pd.read_parquet(checkpoint_path)
    completed = set(df.gid.to_list())
else:
    df = pd.DataFrame(columns=["country", "gid", "pixel_count", "TN", "FP", "FN", "TP"])
    completed = set()

completed.update(removed)


error_path = "../results/country_errors.parquet"

if os.path.exists(error_path):
    error_df = pd.read_parquet(error_path)
    error = set(error_df.gid.to_list())
    completed.update(error)
else:
    error_df = pd.DataFrame(columns=["country", "gid", "error"])

rasters = [GHSLSurface(), BlackMarblePy(frequency="monthly")]
transformations = [{"reproject_match": {"resampling": Resampling.sum}}, {}]
date = datetime.date(2020, 1, 1)

with tqdm(pycountry.countries, desc="Calculating confusion matrix:") as pbar:

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

                province_gdf = gpd.GeoDataFrame(gdf.iloc[[i]].geometry).set_crs("EPSG:4326")

                pipeline = RasterPipeline(province_gdf, date, rasters, transformations)
                ds = pipeline.run()

                # make sure we compare non-nan areas
                mask = get_combined_mask(ds)
                ds = ds.where(mask)

                ghsl_surface_binary = binarize_xarray(ds.ghsl_surface, 50000)
                y_true = get_non_nan_flat_array(ghsl_surface_binary)

                bm_binary = binarize_xarray(ds.black_marble_radiance_monthly, 1.0)
                y_pred = get_non_nan_flat_array(bm_binary)

                assert y_pred.shape == y_true.shape

                conf_mat += confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
                pixels += mask.sum().item()

            data = [country, gid, pixels, *conf_mat.tolist()]
            df = pd.concat([pd.DataFrame([data], columns=df.columns), df], ignore_index=True)
            df.to_parquet(checkpoint_path)

            bm_path = Path(os.path.abspath('')).parent / "data" / "black_marble"
            shutil.rmtree(bm_path)  
            bm_path.mkdir(exist_ok=True)

        except Exception as e:

            data = [country, gid, str(e)]
            error_df = pd.concat([pd.DataFrame([data], columns=error_df.columns), error_df], ignore_index=True)
            error_df.to_parquet(error_path)
