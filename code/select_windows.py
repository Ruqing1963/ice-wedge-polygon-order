"""
select_windows.py -- find gap-free, internally homogeneous analysis windows.

Scans 1.5 km windows over the parsed inventory (../data/iwp.npz), records the
density coefficient of variation (CV, in 200 m cells), low-centred fraction and N,
and writes ../data/window_candidates.npy. Also writes an overview figure.
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

def main(W=1500.0, step=750.0, min_n=4000):
    d = np.load(os.path.join(HERE, "..", "data", "iwp.npz"))
    X, Y, Z = d["x"], d["y"], d["zcen"]
    rows = []
    for e in np.arange(X.min(), X.max() - W, step):
        for n in np.arange(Y.min(), Y.max() - W, step):
            m = (X >= e) & (X < e + W) & (Y >= n) & (Y < n + W); k = int(m.sum())
            if k < min_n:
                continue
            H = np.histogram2d(X[m], Y[m], bins=[7, 7])[0]
            rows.append((e, n, H.std() / H.mean(), (Z[m] < 0).mean(), k))
    R = np.array(rows); np.save(os.path.join(HERE, "..", "data", "window_candidates.npy"), R)
    print(f"{len(R)} candidate windows (N>={min_n}); CV range [{R[:,2].min():.2f},{R[:,2].max():.2f}]")

    # overview figure
    rng = np.random.default_rng(0); idx = rng.choice(len(X), min(120000, len(X)), replace=False)
    fig, ax = plt.subplots(1, 2, figsize=(12, 6.5))
    sc = ax[0].scatter(X[idx] / 1000, Y[idx] / 1000, c=np.clip(Z[idx], -0.3, 0.3), s=1, cmap="RdBu_r", lw=0)
    ax[0].set_aspect("equal"); ax[0].set_xlabel("UTM E (km)"); ax[0].set_ylabel("UTM N (km)")
    ax[0].set_title("(a) IWP centroids (subsample); blue = low-centred")
    plt.colorbar(sc, ax=ax[0], shrink=.6, label="relative centre elevation (m)")
    ax[1].hist(np.clip(Z, -0.4, 0.4), bins=80, color="#534AB7"); ax[1].axvline(0, color="k", ls="--")
    ax[1].set_xlabel("relative centre elevation (m)"); ax[1].set_ylabel("count")
    ax[1].set_title(f"(b) low-centred (<0, {(Z<0).mean()*100:.0f}%) vs high-centred")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "..", "figures", f"iwp_overview.{ext}"), dpi=120, bbox_inches="tight")
    print("overview figure written")

if __name__ == "__main__":
    main()
