# Derived data

## `iwp_window_stats.csv`
Per-window summary over 90 non-overlapping 1.5 km windows (output of
`code/analyze_degradation.py`).

| column | units | description |
|---|---|---|
| `utm_e_m`, `utm_n_m` | m | window origin (UTM Zone 6N) |
| `density_CV` | – | coefficient of variation of counts in 200 m cells (large-scale heterogeneity) |
| `low_centered_frac` | – | fraction of polygons with relative centre elevation < 0 (intact) |
| `median_zcen_m` | m | median relative centre elevation |
| `N` | – | polygons in window |
| `spacing_m` | m | mean centre spacing |
| `gr_peak` | – | first peak height of g(r) (CV-contaminated; see paper) |
| `S_min` | – | minimum of the structure factor (CV-robust intrinsic-order metric) |
| `S_k0` | – | structure factor at the lowest wavevector |
| `numvar_slope` | – | log-log slope of number variance vs window radius |

The raw polygon inventory (`iwp.npz`, ~1.16M centroids) is **not** redistributed;
regenerate it with `code/parse_iwp.py` from the PANGAEA download (see `../DATA_LICENSE.md`).
