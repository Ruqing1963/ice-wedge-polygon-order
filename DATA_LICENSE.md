# Data provenance and licence

This repository redistributes only **derived** per-window summary statistics
(`data/iwp_window_stats.csv`). The raw polygon inventory is **not** redistributed
here; download it from the source below and regenerate with `code/parse_iwp.py`.

## Source
- **Dataset:** Abolt, C.J., Young, M.H. (2019). High-resolution maps of ice wedge
  polygon occurrence and geomorphology [dataset]. PANGAEA,
  https://doi.org/10.1594/PANGAEA.910178
- **Data descriptor:** Abolt, C.J., Young, M.H. (2020). High-resolution mapping of
  spatial heterogeneity in ice wedge polygon geomorphology near Prudhoe Bay, Alaska.
  Scientific Data 7, 87. https://doi.org/10.1038/s41597-020-0423-9
- **File used:** `polygonStatistics.dat` (unified tab-delimited attribute table:
  site, polygon id, centroid UTM-6N x/y (m), area (m^2), relative centre elevation (m)).
  Direct link: https://hs.pangaea.de/Maps/AboltCJ-YoungMH_2019/polygonStatistics.dat
- **Licence:** Creative Commons Attribution 4.0 International (CC-BY-4.0).
  **Required attribution:** "Ice-wedge polygon data from Abolt & Young (2019/2020),
  PANGAEA, CC-BY-4.0."

## Derived data (`data/iwp_window_stats.csv`)
Per-window summary over 90 non-overlapping 1.5 km windows: UTM origin, density CV,
low-centred fraction, median relative centre elevation, N, mean spacing, g(r) peak,
S_min, S(k->0), number-variance slope.
