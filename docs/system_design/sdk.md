# SDK-s

## Datasources

### Landsat 8

- **API Access:** [USGS JSON API](https://m2m.cr.usgs.gov/api/docs/json/)
- **Python Tooling:** [`landsatxplore`](https://github.com/yannforget/landsatxplore)
    - Provides a Python and CLI interface to [EarthExplorer](https://earthexplorer.usgs.gov/) for searching and downloading scenes.
    - **Status:** *Not actively maintained* (for alternatives, see project issues or consider direct API calls).

### EnMAP

- **Portal Access:** [EOWEB Geoportal](https://eoweb.dlr.de/egp/main#mainWindowtabExplore)
- **Note:** No direct API or programmatic download available; files must be placed on the local file system before use.

### SDGSat

- Limited SDGSat coverage is available from pre-downloaded files.
- Data location and structure are detailed in the project’s documentation (see *Filip’s notes*).

??? question "`L4A_A` vs `L4A_B`"

    What is the difference between files with same names, but different endings `L4A_A` and `L4A_B`?

### Black Marble

- **Python SDK:** [`blackmarblepy`](https://github.com/worldbank/blackmarblepy)
- **Note:** Some groups (e.g., ETH Zurich) use a [forked version](https://github.com/Starbix/Infrared_Marble) due to concerns about the original SDK’s reliability and support.

---

## Architecture

### Phase 1: SDK Layer

- Define a shared interface for `raster`, following a similar approach to `blackmarblepy`.
- Read files directly using the existing pre-configured filesystem logic.
- Provide a unified raster method that returns an `xarray.Dataset`.

### Phase 2: API

- Use FastAPI to serve local TIFF files through an API.
- Take inspiration from the Sentinel API or Black Marble API design.

![api](../assets/d2/api.d2)

??? info "Color Coding"

    Phase 1 components are shown in blue, existing elements in grey, and API interfaces planned for Phase 2 are highlighted in purple (subject to available time).

---
