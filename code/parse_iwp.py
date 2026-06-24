"""
parse_iwp.py -- parse the Abolt & Young polygonStatistics.dat into a compact
NumPy archive used by the analysis scripts.

Download the file first (CC-BY-4.0; ~54 MB):
  https://hs.pangaea.de/Maps/AboltCJ-YoungMH_2019/polygonStatistics.dat

Usage:
  python parse_iwp.py /path/to/polygonStatistics.dat
Writes ../data/iwp.npz with arrays x, y (m, UTM-6N), area (m^2), zcen (m).
Columns in the .dat: site, pid, cenx (m), ceny (m), area (m^2), zcen (m).
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

def main(path):
    x = []; y = []; ar = []; z = []
    with open(path, encoding="latin1") as f:
        next(f)
        for ln in f:
            p = ln.rstrip().split("\t")
            if len(p) < 6:
                continue
            x.append(float(p[2])); y.append(float(p[3])); ar.append(float(p[4])); z.append(float(p[5]))
    out = os.path.join(HERE, "..", "data", "iwp.npz")
    np.savez(out, x=np.array(x), y=np.array(y), area=np.array(ar), zcen=np.array(z))
    print(f"wrote {out}: {len(x)} polygons; median area {np.median(ar):.0f} m^2 "
          f"(spacing ~{2*np.sqrt(np.median(ar)/np.pi):.0f} m); "
          f"low-centred {np.mean(np.array(z)<0)*100:.0f}%")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    main(sys.argv[1])
