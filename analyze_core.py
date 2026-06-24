"""
analyze_core.py -- core result on the most uniform window:
pair correlation, structure factor vs four nulls, number variance.
Reads ../data/iwp.npz and ../data/window_candidates.npy (run select_windows.py first).
Writes ../figures/iwp_core.{pdf,png}.
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import stats_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(4)

def main():
    d = np.load(os.path.join(HERE, "..", "data", "iwp.npz")); X, Y = d["x"], d["y"]
    C = np.load(os.path.join(HERE, "..", "data", "window_candidates.npy"))
    cand = C[C[:, 4] > 5000]; e0, n0 = cand[np.argmin(cand[:, 2])][:2]  # most uniform large window
    W = 1500.0
    m = (X >= e0) & (X < e0 + W) & (Y >= n0) & (Y < n0 + W)
    P = np.c_[X[m], Y[m]]; P = P - P.min(0); Lx, Ly = P.max(0); N = len(P)
    rho = N / (Lx * Ly); spacing = 1 / np.sqrt(rho)

    kc, S = sc.structure_factor(P, Lx, Ly)
    _, Sc = sc.structure_factor(sc.csr(N, Lx, Ly, rng), Lx, Ly)
    _, Spv = sc.structure_factor(sc.poisson_voronoi_centroids(N, Lx, Ly, rng), Lx, Ly)
    _, Sin = sc.structure_factor(sc.inhomogeneous_poisson(P, Lx, Ly, rng), Lx, Ly)
    rc, g = sc.pair_correlation(P, Lx, Ly, rng=rng)
    Rs = np.array([15, 25, 40, 60, 90, 140, 210.])
    nv = sc.number_variance(P, Lx, Ly, Rs, rng=rng)
    nvc = sc.number_variance(sc.csr(N, Lx, Ly, rng), Lx, Ly, Rs, rng=rng)
    nvp = sc.number_variance(sc.poisson_voronoi_centroids(N, Lx, Ly, rng), Lx, Ly, Rs, rng=rng)
    slope = np.polyfit(np.log(Rs), np.log(nv), 1)[0]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    ax[0].plot(rc, g, "-o", color="#534AB7", ms=3.5); ax[0].axhline(1, ls="--", color="k", lw=.7)
    ax[0].set_xlabel("r (m)"); ax[0].set_ylabel("g(r)")
    ax[0].set_title(f"(a) pair correlation: hard-core + spacing ~{spacing:.0f} m")
    ax[1].plot(kc, S, "-o", color="#534AB7", ms=3.5, label="IWP centroids", zorder=5)
    ax[1].plot(kc, Spv, "-", color="#2a9d8f", lw=1.4, label="Poisson-Voronoi (random tessellation)")
    ax[1].plot(kc, Sin, "-", color="#D9943A", lw=1.4, label="inhomogeneous Poisson (matched density)")
    ax[1].plot(kc, Sc, "-", color="#aaa", lw=1, label="CSR (Poisson)")
    ax[1].axhline(1, ls="--", color="k", lw=.6); ax[1].set_yscale("log")
    ax[1].set_xlabel("k (1/m)"); ax[1].set_ylabel("S(k)"); ax[1].legend(fontsize=7.5, frameon=False)
    ax[1].set_title("(b) structure factor: deep genuine dip, below random tessellation")
    ax[2].loglog(Rs, nv, "-o", color="#534AB7", ms=4, label=f"IWP (slope {slope:.2f})")
    ax[2].loglog(Rs, nvc, "-", color="#aaa", label="CSR (slope 2)")
    ax[2].loglog(Rs, nvp, "-", color="#2a9d8f", label="Poisson-Voronoi")
    ax[2].loglog(Rs, nv[0] * (Rs / Rs[0]), "k:", lw=1, label="hyperuniform (slope 1)")
    ax[2].set_xlabel("window radius R (m)"); ax[2].set_ylabel(r"number variance $\sigma^2(R)$")
    ax[2].legend(fontsize=7.5, frameon=False); ax[2].set_title("(c) number variance: long-range NOT hyperuniform")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "..", "figures", f"iwp_core.{ext}"), dpi=120, bbox_inches="tight")
    print(f"window N={N} spacing={spacing:.1f}m S_min={np.nanmin(S):.3f}@k={kc[np.nanargmin(S)]:.2f} "
          f"PV_min={np.nanmin(Spv):.3f} numvar_slope={slope:.2f}")

if __name__ == "__main__":
    main()
