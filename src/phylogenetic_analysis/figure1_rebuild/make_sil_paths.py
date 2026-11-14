#!/usr/bin/env python3
"""Convert the potrace SVG silhouettes in ./sil/*.svg into normalised vector paths (sil_paths.npy)
used by build_figure.py.  Each entry: family -> (vertices, path codes, width/height aspect); height = 1."""
import glob, os, re
import numpy as np
from svgpath2mpl import parse_path

out = {}
for f in sorted(glob.glob("sil/*.svg")):
    s = open(f).read()
    d = " ".join(re.findall(r'<path[^>]*\sd="([^"]+)"', s))
    p = parse_path(d)
    V = p.vertices.copy() * np.array([0.1, -0.1])      # potrace: scale(0.1,-0.1)
    V -= V.min(0)
    h = V[:, 1].max()
    V /= h
    V[:, 1] = 1 - V[:, 1]                             # y up (matplotlib axes)
    out[os.path.basename(f)[:-4]] = (V, p.codes.copy(), float(V[:, 0].max()))
np.save("sil_paths.npy", out, allow_pickle=True)
print(len(out), "silhouettes ->", "sil_paths.npy")
