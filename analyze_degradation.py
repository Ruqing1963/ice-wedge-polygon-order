"""
analyze_degradation.py -- 90-window scan and the intact-vs-degraded contrast.

Scans non-overlapping 1.5 km windows, computes per-window CV / low-centred
fraction / g(r)-peak / S_min / S(k->0) / number-variance slope, writes
../results/iwp_window_stats.csv, and the degradation figure ../figures/iwp_degradation.{pdf,png}.
Key point: S_min is the CV-robust intrinsic-order metric used for the contrast.
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import stats_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(9)

def gpeak(P, spacing):
    t = cKDTree(P); rmax = 2.2 * spacing; edges = np.arange(0, rmax, spacing * 0.15)
    rho = len(P) / ((P[:, 0].max() - P[:, 0].min()) * (P[:, 1].max() - P[:, 1].min()))
    sub = P[rng.choice(len(P), min(600, len(P)), replace=False)]; cnt = np.zeros(len(edges) - 1)
    for p in sub:
        ds = np.sqrt(((P[t.query_ball_point(p, rmax)] - p) ** 2).sum(1)); ds = ds[ds > 0]
        cnt += np.histogram(ds, bins=edges)[0]
    return np.nanmax(cnt / (len(sub) * rho * np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)))

def main(W=1500.0):
    d = np.load(os.path.join(HERE, "..", "data", "iwp.npz")); X, Y, Z = d["x"], d["y"], d["zcen"]
    rows = []
    for e in np.arange(X.min(), X.max() - W, W):
        for n in np.arange(Y.min(), Y.max() - W, W):
            m = (X >= e) & (X < e + W) & (Y >= n) & (Y < n + W); k = int(m.sum())
            if k < 4000:
                continue
            P = np.c_[X[m], Y[m]]; P = P - P.min(0); zc = Z[m]; Lx, Ly = P.max(0)
            H = np.histogram2d(P[:, 0], P[:, 1], bins=[7, 7])[0]; cv = H.std() / H.mean()
            sp = 1 / np.sqrt(k / (Lx * Ly)); kc, S = sc.structure_factor(P, Lx, Ly, mmax=45, nbin=38, kmax=0.5)
            gp = gpeak(P, sp)
            Rs = np.array([20, 35, 60, 100.])
            sl = np.polyfit(np.log(Rs), np.log(sc.number_variance(P, Lx, Ly, Rs, rng=rng, nsamp=150)), 1)[0]
            rows.append((e, n, cv, (zc < 0).mean(), np.median(zc), k, sp, gp, np.nanmin(S), S[0], sl))
    R = np.array(rows)
    import csv
    with open(os.path.join(HERE, "..", "results", "iwp_window_stats.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["utm_e_m", "utm_n_m", "density_CV", "low_centered_frac", "median_zcen_m",
                    "N", "spacing_m", "gr_peak", "S_min", "S_k0", "numvar_slope"])
        for r in R:
            w.writerow([f"{r[0]:.0f}", f"{r[1]:.0f}", f"{r[2]:.3f}", f"{r[3]:.3f}", f"{r[4]:.3f}",
                        int(r[5]), f"{r[6]:.1f}", f"{r[7]:.3f}", f"{r[8]:.3f}", f"{r[9]:.2f}", f"{r[10]:.3f}"])
    cv, low, smin = R[:, 2], R[:, 3], R[:, 8]
    print(f"{len(R)} windows. S_min vs low_frac r={np.corrcoef(low,smin)[0,1]:+.2f}, vs CV r={np.corrcoef(cv,smin)[0,1]:+.2f}")

    # figure: S_min vs CV coloured by low%; matched-pair S(k); S_min boxplot
    band = (cv >= 0.25) & (cv < 0.45); lo = band & (low > 0.5); hi = band & (low < 0.25)
    il = np.where(lo)[0][np.argmin(np.abs(cv[np.where(lo)[0]] - 0.36))]
    ih = np.where(hi)[0][np.argmin(np.abs(cv[np.where(hi)[0]] - 0.36))]
    def win(i):
        e, n = R[i, 0], R[i, 1]; m = (X >= e) & (X < e + W) & (Y >= n) & (Y < n + W)
        P = np.c_[X[m], Y[m]]; return P - P.min(0)
    Pl, Ph = win(il), win(ih); kcl, Sl = sc.structure_factor(Pl, *Pl.max(0)); kch, Sh = sc.structure_factor(Ph, *Ph.max(0))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    s0 = ax[0].scatter(cv, smin, c=low * 100, cmap="coolwarm_r", s=30, edgecolor="k", lw=.3)
    ax[0].set_xlabel("density CV (large-scale heterogeneity)"); ax[0].set_ylabel("S_min")
    ax[0].set_title("(a) S_min is CV-robust; colour = low-centred %"); plt.colorbar(s0, ax=ax[0], label="low-centred %")
    ax[1].plot(kch, Sh, "-o", color="#C0392B", ms=3, label=f"high-centred (degraded, {low[ih]*100:.0f}% low)")
    ax[1].plot(kcl, Sl, "-o", color="#2471A3", ms=3, label=f"low-centred (intact, {low[il]*100:.0f}% low)")
    ax[1].axhline(1, ls="--", color="k", lw=.6); ax[1].set_yscale("log")
    ax[1].set_xlabel("k (1/m)"); ax[1].set_ylabel("S(k)"); ax[1].legend(fontsize=8, frameon=False)
    ax[1].set_title("(b) matched CV: order persists in both")
    bp = ax[2].boxplot([smin[band & (low > 0.4)], smin[band & (low < 0.25)]],
                       tick_labels=["low-centred\n(intact)", "high-centred\n(degraded)"], patch_artist=True)
    for p, c in zip(bp["boxes"], ["#2471A3", "#C0392B"]):
        p.set_facecolor(c); p.set_alpha(.5)
    ax[2].set_ylabel("S_min"); ax[2].set_title("(c) order persists in both\n(S_min difference weak, boxes overlap)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "..", "figures", f"iwp_degradation.{ext}"), dpi=120, bbox_inches="tight")
    print("results/iwp_window_stats.csv and figures/iwp_degradation written")

if __name__ == "__main__":
    main()
