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

??? important "Sun-Synchronous Orbit"

    **EnMAP and Landsat 8** are sun-synchronous satellites, acquiring images during consistent daylight hours for each location. As a result, standard data from these satellites do not include nighttime photos. However, one can ask special photos, which might also be during the nighttime.

---

## SDGSAT-1

SDGSAT-1 is equipped with the Glimmer sensor, featuring a panchromatic band at 10-meter resolution.

!!! info "Panchromatic"

    A panchromatic band captures light from the red, green, and blue visible spectrum and merges it into a single channel. This maximizes image brightness and spatial detail, resulting in sharp black-and-white imagery while sacrificing color information for higher ground resolution.

The Glimmer panchromatic imaging mode includes three gain bands [^1]:

- **Panchromatic Low (PL)**:  
  Records lower-intensity light, beneficial for brightly lit scenes to prevent saturation.
- **Panchromatic High (PH)**:  
  Uses higher sensitivity for detecting faint signals, ideal for low-light or nighttime imaging.
- **High Dynamic Range (HDR)**:  
  Merges PL and PH signals to retain details in both bright and dark areas, expanding dynamic range beyond single-gain limitations.

### Technical Details

- **Coordinate Reference System:**  
  EPSG:32631

- **Bit Depth:**  
  12-bit unsigned integer (range: 0–4095)

[^1]: Zhao, Z., Qiu, S., Chen, F., Chen, Y., Qian, Y., Cui, H., Zhang, Y., Khoramshahi, E., & Qiu, Y. (2023). Vessel Detection with SDGSAT-1 Nighttime Light Images. Remote Sensing, 15(17), 4354. [https://doi.org/10.3390/rs15174354](https://doi.org/10.3390/rs15174354)
