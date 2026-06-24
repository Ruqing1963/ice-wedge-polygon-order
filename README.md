# Ice-wedge polygon order

Structure-factor analysis of a million-polygon ice-wedge network. Treating the
centres of $1.16\times10^{6}$ lidar-mapped ice-wedge polygons (Prudhoe Bay, Alaska)
as a planar point pattern, we ask how **ordered** the pattern is — and, crucially,
whether the crack network is more ordered than a **random** space-filling tessellation.

This repository contains the code, derived per-window statistics and figures for:

> R. Chen (2026). *Ice-wedge polygon centres are more ordered than a random tessellation:
> structure-factor analysis of a million-polygon network with heterogeneity-limited
> long-range order.* ([`paper/iwp_paper.pdf`](paper/iwp_paper.pdf))

## Key results

| Quantity | Value |
|---|---|
| Polygons analysed (Prudhoe Bay, AK) | 1,162,700 |
| Characteristic spacing (median area 182 m²) | ~15 m |
| Pair correlation | hard-core exclusion + peak at ~15–17 m |
| Structure-factor minimum (uniform windows) | **S_min ≈ 0.09–0.10** at k ≈ 0.16 m⁻¹ |
| vs inhomogeneous-Poisson null (matched density) | S_min far below it in **every** window → genuine |
| vs Poisson–Voronoi null (random tessellation) | S_min below it in uniform windows → **more ordered than random** |
| Long-range (number variance) | σ²(R) ∝ R²·⁷ (> R²) → **not hyperuniform** |
| Thermokarst degradation (low- vs high-centred) | order **persists** (S_min similar, overlapping) |

**Interpretation.** Ice-wedge polygon centres are a locally rigid, intermediate-range-ordered
disordered pattern with a hyperuniform *tendency* (S_min ≈ 0.1, deeper than other geological
spacing patterns). Long-range hyperuniformity is **not** confirmed: real geomorphic
heterogeneity (polygon-size/coverage variation) dominates the largest scales, and the
suppression deepens in more uniform windows — a heterogeneity-limited, not intrinsic, floor.
Thermokarst degradation inverts the microtopographic relief (low→high-centred) but preserves
the planimetric network order.

## Repository layout

```
code/
  stats_common.py        structure factor, g(r), number variance, and nulls
                         (CSR, Poisson-Voronoi centroids, inhomogeneous-Poisson, hexagonal)
  parse_iwp.py           polygonStatistics.dat -> data/iwp.npz
  select_windows.py      tile coverage + uniform-window candidates + overview figure
  analyze_core.py        core result figure (g(r), S(k) vs nulls, number variance)
  analyze_degradation.py 90-window scan + intact-vs-degraded contrast + window-stats CSV
data/
  iwp_window_stats.csv   derived per-window statistics (90 windows; see data/README.md)
figures/   the three paper figures (PDF + PNG)
results/   iwp_window_stats.csv
paper/     manuscript (PDF + LaTeX source + figures)
```

## Reproduce

```bash
pip install -r requirements.txt

# 1. download the raw inventory (CC-BY-4.0, ~54 MB):
#    https://hs.pangaea.de/Maps/AboltCJ-YoungMH_2019/polygonStatistics.dat
cd code
python parse_iwp.py /path/to/polygonStatistics.dat   # -> data/iwp.npz
python select_windows.py                              # -> window candidates + overview
python analyze_core.py                                # -> figures/iwp_core.*
python analyze_degradation.py                         # -> figures/iwp_degradation.* + results CSV
```

## Data source and licence

Ice-wedge polygon inventory: Abolt & Young (2019/2020), PANGAEA
[doi:10.1594/PANGAEA.910178](https://doi.org/10.1594/PANGAEA.910178) and *Scientific Data*
**7**, 87 ([doi:10.1038/s41597-020-0423-9](https://doi.org/10.1038/s41597-020-0423-9)),
**CC-BY-4.0**. This repository redistributes only derived per-window statistics; see
[`DATA_LICENSE.md`](DATA_LICENSE.md). Code is MIT ([`LICENSE`](LICENSE)).

## Companion studies

Two other geological map-pattern spacing analyses with the same tools:
[salt-structure-spacing](https://github.com/Ruqing1963/salt-structure-spacing)
(Zenodo [10.5281/zenodo.20834957](https://doi.org/10.5281/zenodo.20834957); buoyancy-driven,
isotropic hard-core regular) and
[drumlin-anisotropic-spacing](https://github.com/Ruqing1963/drumlin-anisotropic-spacing)
(Zenodo [10.5281/zenodo.20835748](https://doi.org/10.5281/zenodo.20835748); shear-driven,
anisotropic). Among the three, ice-wedge polygons show the deepest intermediate-range
suppression.
