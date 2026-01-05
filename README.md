# Conflict Monitoring - Night Time Light

[![Project Status: Research](https://img.shields.io/badge/Status-Research-blue.svg)](https://github.com/epfl-eceo/conflict-monitoring-ntl)

Semester project at the **EPFL ECEO Lab** focusing on the use of Nighttime Light (NTL) satellite data to monitor and characterize human activity in conflict-affected regions.

## Installation

Make sure you have **uv** installed on your system. If not, install it via:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then install project dependencies:

```bash
uv sync
```

## Documentation

Documentation dependencies are managed in the **docs** dependency group. You can serve the MkDocs site locally or check it in [production](https://jankokla.github.io/conflict-monitoring-ntl/):

```bash
uv run --only-group docs mkdocs serve
```

This will start a local development server, usually available at:

```
http://127.0.0.1:8000/
```

## System Architecture

The codebase is organized into a modular pipeline and specialized analysis scripts designed for reproducible geospatial research:

### Core Modules

* **`satellites.py`**: A robust data loader interface (GEE, USGS M2M, and local drivers).
* **`pipeline.py`**: Orchestration engine for spatial alignment and data augmentation.
* **`utils.py`**: Core utilities for GADM admin boundaries and precision/recall metrics.
* **`plotting.py`**: Visualization suite for mapping and statistical exploration.

### Scripts
* **`bm_ghsl_match_global.py`**: Evaluates NTL binarization sensitivity against **GHSL Built-up Surface**.
* **`bm_ghsl_smod_global.py`**: Investigates urbanization baselines using the **GHSL Settlement Model (GHS-SMOD)**.

## Examples

The project is designed to handle complex geospatial workflows through the `RasterPipeline`. Below are examples of how to combine sensors and apply spatial transformations.

```python
import datetime
from rasterio.enums import Resampling
from conflict_monitoring_ntl.satellites import BlackMarblePy, GHSLSurface
from conflict_monitoring_ntl.transform import RasterPipeline

# define sensors
rasters = [GHSLSurface(), BlackMarblePy(frequency="monthly")]

# define transformations:
# 1. reproject GHSL to match Black Marble (base_index=-1 by default)
# 2. no transformation for Black Marble (it serves as the spatial reference)
transformations = [{"reproject_match": {"resampling": Resampling.sum}}, {}]

pipeline = RasterPipeline(province_gdf, datetime.date(2020, 1, 1), rasters, transformations)
ds = pipeline.run()
```