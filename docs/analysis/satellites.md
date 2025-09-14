# Satellites

The following satellites are widely used in nighttime light (NTL) analysis and urban monitoring:

- **[Landsat 8](https://landsat.gsfc.nasa.gov/satellites/landsat-8/)**: NASA–USGS satellite for global land monitoring, optimized for daytime land cover and change analysis.
- **[SDGSAT-1](https://www.sdgsat.ac.cn/satellite/intro)**: Chinese Academy of Sciences satellite designed for SDG studies, featuring dedicated nighttime lights imaging.
- **[EnMAP](https://www.enmap.org/)**: German hyperspectral satellite focused on environmental, land, and water quality mapping, mainly under daylight.
- **[Black Marble (VIIRS DNB)](https://blackmarble.gsfc.nasa.gov/)**: NASA’s global nighttime lights data product, providing daily radiance maps.

??? info "Spectral Bands"

    Only spectral bands relevant to nighttime light (NTL) analysis are shown below.

| Satellite    | Spectral Bands            | Imaging Time (local time)   | Coverage Frequency    | Resolution         | Products       |
|--------------|---------------------------|-----------------------------|----------------------|--------------------|-----------------|
| ***Landsat 8***    | Coastal aerosol + RGB     | 10:00 AM (daylight only)    | Every 16 days        | 30 m               |                 |
| ***SDGSAT-1***     | Glimmer (Panchromatic + RGB) | 9:30 PM (nighttime lights) | Every 11 days        | P: 10 m, RGB: 40 m |               |
| ***EnMAP***        | RGB (VNIR bands)          | 11:00 AM (daylight only)    | Every 27 days        | 30 m               |                 |
| ***Black Marble*** | Day/Night Band (DNB)      | 1:30 AM (nightly)           | Daily (global)       | 500 m              |                 |

??? question "Sun-Synchronous Orbit"

    **EnMAP and Landsat 8** are sun-synchronous satellites, acquiring images during consistent daylight hours for each location. As a result, standard data from these satellites do not include nighttime photos.
