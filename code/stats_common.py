"""
stats_common.py -- shared point-pattern statistics for the ice-wedge polygon
order analysis. Pure NumPy/SciPy.
"""
import numpy as np
from scipy.spatial import cKDTree, Voronoi


def structure_factor(pts, Lx, Ly, mmax=70, nbin=46, kmax=0.55):
    """Radially averaged structure factor at the window's allowed wavevectors
    k = 2*pi*(m/Lx, n/Ly) (which suppresses the window-shape artefact)."""
    q = pts - pts.mean(0)
    ms = np.arange(-mmax, mmax + 1)
    KX, KY = np.meshgrid(2 * np.pi * ms / Lx, 2 * np.pi * ms / Ly)
    K = np.c_[KX.ravel(), KY.ravel()]; kk = np.hypot(K[:, 0], K[:, 1])
    sel = (kk > 0) & (kk <= kmax); K = K[sel]; kk = kk[sel]
    S = np.empty(len(K))
    for i in range(0, len(K), 4000):
        S[i:i + 4000] = np.abs(np.exp(-1j * (K[i:i + 4000] @ q.T)).sum(1)) ** 2 / len(q)
    bins = np.linspace(0, kmax, nbin); idx = np.digitize(kk, bins); kc = 0.5 * (bins[1:] + bins[:-1])
    return kc, np.array([S[idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, len(bins))])


def pair_correlation(pts, Lx, Ly, rmax=70, dr=2.0, rng=None, nsub=2000):
    """Isotropic g(r) from a subsample of interior points (robust, fast)."""
    rng = rng or np.random.default_rng(0)
    rho = len(pts) / (Lx * Ly); t = cKDTree(pts)
    edges = np.arange(0, rmax, dr); rc = 0.5 * (edges[1:] + edges[:-1])
    sub = pts[rng.choice(len(pts), min(nsub, len(pts)), replace=False)]
    cnt = np.zeros(len(rc))
    for p in sub:
        ds = np.sqrt(((pts[t.query_ball_point(p, rmax)] - p) ** 2).sum(1)); ds = ds[ds > 0]
        cnt += np.histogram(ds, bins=edges)[0]
    return rc, cnt / (len(sub) * rho * np.pi * (edges[1:] ** 2 - edges[:-1] ** 2))


def number_variance(pts, Lx, Ly, Rs, rng=None, nsamp=600):
    """Variance of point counts in disks of radius R (interior-sampled)."""
    rng = rng or np.random.default_rng(0); t = cKDTree(pts)
    return np.array([np.var([len(t.query_ball_point([a, b], R))
                             for a, b in zip(rng.uniform(R, Lx - R, nsamp), rng.uniform(R, Ly - R, nsamp))])
                     for R in Rs])


def csr(N, Lx, Ly, rng):
    return np.c_[rng.uniform(0, Lx, N), rng.uniform(0, Ly, N)]


def poisson_voronoi_centroids(N, Lx, Ly, rng):
    """Centroids of a Poisson-Voronoi tessellation = a *random* space-filling
    tessellation of matched density (the key null)."""
    s = np.c_[rng.uniform(0, Lx, N), rng.uniform(0, Ly, N)]; vor = Voronoi(s); c = []
    for reg in (vor.regions[i] for i in vor.point_region):
        if -1 in reg or len(reg) < 3:
            continue
        c.append(vor.vertices[reg].mean(0))
    c = np.array(c); return c[(c[:, 0] > 0) & (c[:, 0] < Lx) & (c[:, 1] > 0) & (c[:, 1] < Ly)]


def inhomogeneous_poisson(pts, Lx, Ly, rng, bw=120.0):
    """Density-matched inhomogeneous-Poisson null: kernel-resample observed
    positions (reproduces large-scale density variation, destroys short-range order)."""
    q = pts[rng.integers(0, len(pts), len(pts))] + rng.normal(0, bw, (len(pts), 2))
    return q[(q[:, 0] > 0) & (q[:, 0] < Lx) & (q[:, 1] > 0) & (q[:, 1] < Ly)]


def hex_lattice(N, Lx, Ly):
    a = np.sqrt(2 / (np.sqrt(3) * (N / (Lx * Ly)))); pts = []; j = 0; yv = 0.0
    while yv < Ly:
        xv = (a / 2 if j % 2 else 0.0)
        while xv < Lx:
            pts.append((xv, yv)); xv += a
        yv += a * np.sqrt(3) / 2; j += 1
    return np.array(pts)
