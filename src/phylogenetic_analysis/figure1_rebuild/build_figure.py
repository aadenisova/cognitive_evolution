#!/usr/bin/env python3
"""Figure 1 (tree + phylogenetic signal + foreground/background), rebuilt.
Inputs (all in ./):  tree_edges.csv tree_tips.csv tree_depth.csv  (exported from R, see 04_export_tree.R)
                     data/tree_vgp_0_vgp.nwk  signal_50.csv  boxplot_data.csv  wilcox.csv  sil_paths.npy
"""
import itertools, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.patches import PathPatch, Polygon, Rectangle, Wedge
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
import colorsys

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.5, "axes.linewidth": 0.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "mathtext.default": "regular",
})
GREEN = "#35B779"; GREY = "#BFBFBF"; INK = "#2B2B2B"; VGP_BLUE = "#1E9BF0"

# ------------------------------------------------------------------ data
E = pd.read_csv("tree_edges.csv"); T = pd.read_csv("tree_tips.csv"); DEP = pd.read_csv("tree_depth.csv")
T["Order"] = T.Order.str.upper(); T["Family"] = T.Family.str.upper()
ntip = len(T); depth = DEP.depth.values; D = depth.max()
children = {}
for p, c in zip(E.parent, E.child): children.setdefault(p, []).append(c)

# VGP tips
vgp_tips = set(re.findall(r"[(,]([A-Za-z_]+):", open("data/tree_vgp_0_vgp.nwk").read()))
is_vgp = T.Species.isin(vgp_tips).values
print("VGP tips in tree:", is_vgp.sum())

# ------------------------------------------------------------------ layout (angles)
ANG0, ANG1 = np.deg2rad(2.0), np.deg2rad(178.0)          # tip 0 at the right, last tip at the left
tip_ang = ANG0 + (np.arange(ntip) + 0.5) / ntip * (ANG1 - ANG0)
tip_ang = ANG0 + (np.arange(ntip)) / (ntip - 1) * (ANG1 - ANG0)
node_ang = {i + 1: tip_ang[i] for i in range(ntip)}
def get_ang(n):
    if n in node_ang: return node_ang[n]
    a = [get_ang(c) for c in children[n]]
    node_ang[n] = 0.5 * (min(a) + max(a)); return node_ang[n]
sys.setrecursionlimit(10000)
root = ntip + 1; get_ang(root)

# ------------------------------------------------------------------ geometry in inches
R = 3.85                      # tree radius (tips)
rad = lambda node: R * depth[node - 1] / D
CMAP = mcolors.LinearSegmentedColormap.from_list("innov", ["#5C6E91", "#5F93A8", "#8FBFA0", "#E9D758", "#F29E1F", "#D7301F", "#8E0F2D"], N=512)
vmin, vmax = -0.62, 2.7
norm = mcolors.PowerNorm(gamma=1.15, vmin=vmin, vmax=vmax, clip=True)

fig_w, fig_h = 12.0, 10.65
fig = plt.figure(figsize=(fig_w, fig_h))
A_BOT = 3.7 / fig_h
axA = fig.add_axes([0, A_BOT, 1, 1 - A_BOT])
axA.set_xlim(-fig_w / 2, fig_w / 2); axA.set_ylim(-1.05, (1 - A_BOT) * fig_h - 1.05); axA.set_aspect("equal"); axA.axis("off")

# --- time grid: faint concentric arcs every 20 Ma (behind the tree)
th = np.linspace(0, np.pi, 300)


# --- branches (colour = ancestral state reconstruction, interpolated along radial segments)
segs, cols, sval = [], [], []
for p, c, pv, cv in zip(E.parent, E.child, E.pv, E.cv):
    a = node_ang[c]; r0, r1 = rad(p), rad(c)
    k = max(1, int(np.ceil((r1 - r0) / 0.06)))
    rr = np.linspace(r0, r1, k + 1); vv = np.linspace(pv, cv, k + 1)
    for i in range(k):
        segs.append([(rr[i] * np.cos(a), rr[i] * np.sin(a)), (rr[i + 1] * np.cos(a), rr[i + 1] * np.sin(a))])
        cols.append(CMAP(norm(0.5 * (vv[i] + vv[i + 1])))); sval.append(0.5 * (vv[i] + vv[i + 1]))
nodeval = {}
for p, pv in zip(E.parent, E.pv): nodeval[p] = pv
for n, ch in children.items():
    a = [node_ang[c] for c in ch]; lo, hi = min(a), max(a)
    tt = np.linspace(lo, hi, max(2, int(np.ceil((hi - lo) / np.deg2rad(0.8))) + 1)); r = rad(n)
    pts = np.c_[r * np.cos(tt), r * np.sin(tt)]
    for i in range(len(pts) - 1):
        segs.append([tuple(pts[i]), tuple(pts[i + 1])]); cols.append(CMAP(norm(nodeval[n]))); sval.append(nodeval[n])
order_ = np.argsort(sval, kind="stable")   # draw high-innovation branches last so they stay visible
lw_ = [0.5 + 0.35 * float(norm(sval[k]) > 0.55) for k in order_]
axA.add_collection(LineCollection([segs[k] for k in order_], colors=[cols[k] for k in order_], linewidths=lw_, capstyle="round", zorder=2))

# ------------------------------------------------------------------ rings
ord_runs, fam_runs = [], []
for key, grp in itertools.groupby(enumerate(T.Order), key=lambda t: t[1]):
    idx = [g[0] for g in grp]; ord_runs.append((key, idx[0], idx[-1]))
for key, grp in itertools.groupby(enumerate(T.Family), key=lambda t: t[1]):
    idx = [g[0] for g in grp]; fam_runs.append((key, idx[0], idx[-1]))
half = (ANG1 - ANG0) / (ntip - 1) / 2

def ring(a0, a1, r0, r1, fc, z=3):
    tt = np.linspace(a0, a1, max(3, int(np.ceil((a1 - a0) / np.deg2rad(0.5))) + 1))
    xy = np.r_[np.c_[r1 * np.cos(tt), r1 * np.sin(tt)], np.c_[r0 * np.cos(tt[::-1]), r0 * np.sin(tt[::-1])]]
    axA.add_patch(Polygon(xy, closed=True, fc=fc, ec=fc, lw=0.15, zorder=z))

# Order colours: muted hues, neighbours kept apart by stepping through the hue wheel
norder = len(ord_runs); step = 11
def order_color(i, v=0.74, s=0.42):
    h = ((i * step) % norder) / norder
    return colorsys.hsv_to_rgb(h, s, v)
OCOL = {o: order_color(i) for i, (o, _, _) in enumerate(ord_runs)}
OTXT = {o: order_color(i, v=0.52, s=0.62) for i, (o, _, _) in enumerate(ord_runs)}

R_VGP = R + 0.075; R_F0, R_F1 = R + 0.13, R + 0.165; R_O0, R_O1 = R + 0.19, R + 0.28
for i, (f, i0, i1) in enumerate(fam_runs):
    ring(tip_ang[i0] - half, tip_ang[i1] + half, R_F0, R_F1, "#C9C9C9" if i % 2 == 0 else "#E6E6E6")
for o, i0, i1 in ord_runs:
    ring(tip_ang[i0] - half, tip_ang[i1] + half, R_O0, R_O1, OCOL[o])

# VGP markers: small filled circles just beyond the tips
for i in np.where(is_vgp)[0]:
    axA.plot(R_VGP * np.cos(tip_ang[i]), R_VGP * np.sin(tip_ang[i]), "o", ms=2.6, mfc=VGP_BLUE, mec="white", mew=0.25, zorder=6)

# ------------------------------------------------------------------ label placement helpers
def pav_spread(desired, min_gap, lo=None, hi=None):
    """Positions closest (least squares) to `desired` with x[i+1]-x[i] >= min_gap[i]  (pool-adjacent-violators)."""
    n = len(desired); off = np.r_[0, np.cumsum(min_gap)]
    y = np.asarray(desired) - off
    blocks = [[y[i], 1, i, i] for i in range(n)]
    out = []
    for b in blocks:
        out.append(b)
        while len(out) > 1 and out[-2][0] > out[-1][0]:
            b2 = out.pop(); b1 = out.pop()
            w = b1[1] + b2[1]; out.append([(b1[0] * b1[1] + b2[0] * b2[1]) / w, w, b1[2], b2[3]])
    x = np.zeros(n)
    for m, w, i0, i1 in out: x[i0:i1 + 1] = m
    return x + off

def put_radial_text(r, ang, s, color, fs, style="normal", weight="normal"):
    deg = np.rad2deg(ang)
    if deg <= 90: rot, ha = deg, "left"
    else: rot, ha = deg + 180, "right"
    axA.text(r * np.cos(ang), r * np.sin(ang), s, rotation=rot, rotation_mode="anchor", ha=ha, va="center",
             fontsize=fs, color=color, style=style, fontweight=weight, zorder=8)

# ------------------------------------------------------------------ silhouettes: exactly one per Order, label beside it
from PIL import Image
from scipy import ndimage
SIL = np.load("sil_paths.npy", allow_pickle=True).item()
ORDER_SIL = {   # Order -> family silhouette (largest family in the data that has a silhouette)
 "ACCIPITRIFORMES": "ACCIPITRIDAE", "CATHARTIFORMES": "CATHARTIDAE", "CORACIIFORMES": "MEROPIDAE", "PICIFORMES": "PICIDAE",
 "BUCEROTIFORMES": "BUCEROTIDAE", "TROGONIFORMES": "TROGONIDAE", "PSITTACIFORMES": "PSITTACIDAE", "PASSERIFORMES": "CORVIDAE",
 "FALCONIFORMES": "FALCONIDAE", "CAPRIMULGIFORMES": "swift.png", "PROCELLARIIFORMES": "PROCELLARIIDAE", "SPHENISCIFORMES": "SPHENISCIDAE",
 "PELECANIFORMES": "ARDEIDAE", "SULIFORMES": "SULIDAE", "CICONIIFORMES": "CICONIIDAE", "GAVIIFORMES": "GAVIIDAE",
 "PHAETHONTIFORMES": "PHAETHONTIDAE", "CHARADRIIFORMES": "SCOLOPACIDAE", "GRUIFORMES": "RALLIDAE", "OPISTHOCOMIFORMES": "OPISTHOCOMIDAE",
 "COLUMBIFORMES": "COLUMBIDAE", "PTEROCLIFORMES": "PTEROCLIDAE", "CUCULIFORMES": "CUCULIDAE", "OTIDIFORMES": "OTIDIDAE",
 "PODICIPEDIFORMES": "PODICIPEDIDAE", "PHOENICOPTERIFORMES": "PHOENICOPTERIDAE", "GALLIFORMES": "PHASIANIDAE", "ANSERIFORMES": "ANATIDAE",
 "STRUTHIONIFORMES": "APTERYGIDAE"}
assert set(ORDER_SIL) == {o for o, _, _ in ord_runs}
SIL_H, SIL_WMAX = 0.215, 0.33
EXTRA_SHRINK = {"PODICIPEDIFORMES": 1.5, "CATHARTIFORMES": 1.5, "STRUTHIONIFORMES": 2.0}
SHRINK = {"PODICIPEDIFORMES", "GALLIFORMES", "PICIFORMES", "COLUMBIFORMES", "PTEROCLIFORMES", "OTIDIFORMES", "CHARADRIIFORMES"}
def sil_aspect(key):
    if key.endswith(".png"):
        im = Image.open("sil/" + key); return im.size[0] / im.size[1]
    return SIL[key][2]
osz = {}
for o, _, _ in ord_runs:
    w = sil_aspect(ORDER_SIL[o]); sc = min(SIL_H, SIL_WMAX / w); osz[o] = (sc * w, sc)
    f_ = (1.5 if o in SHRINK else 1.0) * EXTRA_SHRINK.get(o, 1.0)
    osz[o] = (sc * w / f_, sc / f_)
TIER_R = [R + 0.48, R + 0.80]; R_LAB = R + 1.06; FS_O = 6.6
o_names = [o for o, _, _ in ord_runs]
o_mid = np.array([0.5 * (tip_ang[i0] + tip_ang[i1]) for _, i0, i1 in ord_runs])
n_o = len(o_names); tier_of = [i % 2 for i in range(n_o)]
text_gap = (FS_O / 72 * 1.4) / R_LAB
def req(i, j):   # angular clearance needed between same-tier silhouettes i, j
    r = TIER_R[tier_of[i]]; return ((osz[o_names[i]][0] + osz[o_names[j]][0]) / 2 + 0.06) / r
gaps = []
for i in range(n_o - 1):
    cands = [text_gap]
    if i - 1 >= 0: cands.append(0.5 * req(i - 1, i + 1))
    if i + 2 < n_o: cands.append(0.5 * req(i, i + 2))
    gaps.append(max(cands))
o_pos = pav_spread(o_mid, np.array(gaps))
swift_rgba = None
for i, o in enumerate(o_names):
    sw, sh = osz[o]; r = TIER_R[tier_of[i]]; a = o_pos[i]; am = o_mid[i]
    col = OTXT[o]; cx, cy = r * np.cos(a), r * np.sin(a)
    # leader: ring -> silhouette -> (label sits just beyond on the same ray)
    r0 = R_O1 + 0.02; rk = R + 0.34; r1 = r - sh / 2 - 0.03
    axA.plot([r0 * np.cos(am), r1 * np.cos(a)], [r0 * np.sin(am), r1 * np.sin(a)],
             color=col, lw=0.4, alpha=0.9, zorder=4, solid_joinstyle="round")
    key = ORDER_SIL[o]
    if key.endswith(".png"):
        im = np.asarray(Image.open("sil/" + key).convert("RGBA")).astype(float)
        al = im[..., 3] / 255.0 if im[..., 3].min() < 255 else 1 - im[..., :3].mean(-1) / 255.0
        halo = ndimage.binary_dilation(al > 0.3, iterations=6).astype(float)
        rgb = np.array(mcolors.to_rgb(col))
        for layer, c_, z in ((halo, np.ones(3), 7.0), (al, rgb, 7.5)):
            arr = np.zeros(al.shape + (4,)); arr[..., :3] = c_; arr[..., 3] = layer
            axA.imshow(arr, extent=[cx - sw / 2, cx + sw / 2, cy - sh / 2, cy + sh / 2], zorder=z, interpolation="bilinear")
    else:
        V, C, w = SIL[key]
        Vs = V.copy(); Vs[:, 0] = (Vs[:, 0] - w / 2) * sh + cx; Vs[:, 1] = (Vs[:, 1] - 0.5) * sh + cy
        axA.add_patch(PathPatch(Path(Vs, C), fc=col, ec="white", lw=0.9, zorder=7, joinstyle="round"))
        axA.add_patch(PathPatch(Path(Vs, C), fc=col, ec="none", zorder=7.5))
    put_radial_text(R_LAB, a, o.capitalize(), col, FS_O)

# ------------------------------------------------------------------ time axis (as in Stiller et al. 2024, Fig. 1b):
# ticks every 20 Myr, dashed arcs at the K-Pg (66.0 Ma) and Pg-Neogene (23.0 Ma) boundaries, labelled with arrows below the axis
Y_AX = -0.16
axA.plot([-R, 0], [Y_AX, Y_AX], color=INK, lw=0.5, zorder=5)   # time axis on the LEFT, below the longest (palaeognath) branch
for Tm in (0, 20, 40, 60, 80, 100):
    x = -R * (1 - Tm / D)
    axA.plot([x, x], [Y_AX, Y_AX + 0.05], color=INK, lw=0.5, zorder=5)
    axA.text(x, Y_AX - 0.07, str(Tm), ha="center", va="top", fontsize=6.3, color=INK)
th_k = np.linspace(0, np.pi, 300)
for Tb, lab in ((66.04, "K–Pg boundary"), (23.03, "Pg–Neogene boundary")):
    rb = R * (1 - Tb / D); xb = -rb
    axA.plot(rb * np.cos(th_k), rb * np.sin(th_k), color="#222222", lw=0.6, ls=(0, (4, 3)), zorder=1)
    axA.annotate("", xy=(xb, Y_AX - 0.22), xytext=(xb, Y_AX - 0.38), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.6, mutation_scale=5), zorder=6)
    axA.text(xb, Y_AX - 0.41, lab, ha="center", va="top", fontsize=6.0, color=INK)
axA.text(-R * 0.5, Y_AX - 0.62, "Time (Ma)", ha="center", va="top", fontsize=6.5, color=INK)

# ------------------------------------------------------------------ colour bar + legend (left, below the baseline)
cb_w, cb_h = 1.45, 0.075; cb_x0 = R - cb_w; cb_y0 = Y_AX - 0.12
grad = np.linspace(0, 1, 256)[None, :]
axA.imshow(CMAP(norm(np.linspace(vmin, vmax, 256)))[None, :, :], extent=[cb_x0, cb_x0 + cb_w, cb_y0, cb_y0 + cb_h], aspect="auto", zorder=5)
axA.set_xlim(-fig_w / 2, fig_w / 2); axA.set_ylim(-1.05, (1 - A_BOT) * fig_h - 1.05)
for v in (0, 1, 2):
    x = cb_x0 + (v - vmin) / (vmax - vmin) * cb_w
    axA.plot([x, x], [cb_y0, cb_y0 - 0.03], color=INK, lw=0.5, zorder=5)
    axA.text(x, cb_y0 - 0.06, str(v), ha="center", va="top", fontsize=6, color=INK)
axA.text(cb_x0 + cb_w / 2, cb_y0 - 0.26, "Innovation rate (ancestral-state reconstruction)", ha="center", va="top", fontsize=6.5, color=INK)
lx = cb_x0 - 0.62
axA.plot([lx], [cb_y0 + cb_h / 2], "o", ms=3.6, mfc=VGP_BLUE, mec="white", mew=0.3, zorder=6)
_logo = np.asarray(Image.open("data/VGP_tag.png").convert("RGBA"))
_h = 0.128; _w = _h * _logo.shape[1] / _logo.shape[0]
axA.imshow(_logo, extent=[lx + 0.07, lx + 0.07 + _w, cb_y0 + cb_h / 2 - _h / 2, cb_y0 + cb_h / 2 + _h / 2], zorder=6, interpolation="bilinear")
axA.set_xlim(-fig_w / 2, fig_w / 2); axA.set_ylim(-1.05, (1 - A_BOT) * fig_h - 1.05)
# ================================================================== panels B and C
sig = pd.read_csv("signal_50.csv")
name_map = {"Technical Innovation": "Technical innovation", "Food type innovation": "Food innovation", "Total innovation": "Total innovation",
            "Brain size": "Brain size", "Relative brain size": "Relative brain size", "Mass": "Body mass", "Diet Breadth": "Diet breadth",
            "Habitat Breadth": "Habitat breadth", "Generation Length": "Generation length", "Range Size": "Range size"}
sig["label"] = sig.Trait.map(name_map)
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]; return np.minimum(q, 1)
sig["Kq"] = bh(sig.K_p.values); sig["Lq"] = bh(sig.lambda_p.values)
sig = sig.sort_values("K", ascending=False).reset_index(drop=True)
INNOV = {"Technical innovation", "Food innovation", "Total innovation"}
bar_c = [GREEN if l in INNOV else "#7A8B99" for l in sig.label]
x = np.arange(len(sig))

left, wB = 0.065, 0.50
inch = lambda x0, y0, w, h: [x0 / fig_w, y0 / fig_h, w / fig_w, h / fig_h]
axK = fig.add_axes(inch(0.95, 2.15, 5.65, 1.25)); axL = fig.add_axes(inch(0.95, 0.82, 5.65, 1.25))
def star(q): return "*" if q < 0.05 else "ns"
# K on a log scale, bars drawn from the Brownian-motion expectation K = 1
axK.set_yscale("log"); axK.set_ylim(0.025, 80)
for xi, k, c in zip(x, sig.K, bar_c):
    axK.bar(xi, k - 1, bottom=1, width=0.68, color=c, lw=0, zorder=3)
axK.axhline(1, color=INK, lw=0.6, ls=(0, (3, 2)), zorder=4)
for xi, k, q in zip(x, sig.K, sig.Kq):
    if k >= 1: axK.text(xi, k * 1.12, star(q), ha="center", va="bottom", fontsize=13 if q < .05 else 9.5, color=INK if q < .05 else "#777")
    else: axK.text(xi, k / 1.12, star(q), ha="center", va="top", fontsize=13 if q < .05 else 9.5, color=INK if q < .05 else "#777")
axK.set_yticks([0.03, 0.1, 0.3, 1, 3, 10, 30]); axK.set_yticklabels(["0.03", "0.1", "0.3", "1", "3", "10", "30"])
axK.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axK.set_ylabel("Blomberg's K", fontsize=8.5); axK.set_xticks(x); axK.set_xticklabels([])
axK.text(len(sig) - 0.45, 1.13, "Brownian motion", fontsize=6.5, color="#555", ha="right", va="bottom", style="italic")
axL.bar(x, sig["lambda"], width=0.68, color=bar_c, lw=0, zorder=3)
axL.axhline(1, color=INK, lw=0.6, ls=(0, (3, 2)), zorder=4)
for xi, l, q in zip(x, sig["lambda"], sig.Lq): axL.text(xi, l + 0.012, star(q), ha="center", va="bottom", fontsize=13 if q < .05 else 9.5, color=INK if q < .05 else "#777")
axL.set_ylim(0, 1.22); axL.set_yticks([0, 0.25, 0.5, 0.75, 1]); axL.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
axL.set_ylabel("Pagel's λ", fontsize=8.5)
axL.set_xticks(x); axL.set_xticklabels(sig.label, rotation=38, ha="right", rotation_mode="anchor", fontsize=7.8)
for ax in (axK, axL):
    ax.set_xlim(-0.6, len(sig) - 0.4); ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="x", length=0 if ax is axK else 2.2)
    ax.grid(axis="y", color="#EEEEEE", lw=0.4, zorder=0); ax.set_axisbelow(True)
axK.spines["bottom"].set_visible(False)
# legend for bar colours
for i, (lab, c) in enumerate([("Innovation traits", GREEN), ("Morphology / life history / ecology", "#7A8B99")]):
    axK.add_patch(Rectangle((len(sig) - 6.2 if i == 0 else len(sig) - 6.2, 0), 0, 0))
hand = [Rectangle((0, 0), 1, 1, fc=GREEN), Rectangle((0, 0), 1, 1, fc="#7A8B99")]
axK.legend(hand, ["Innovation traits", "Other traits"], loc="upper right", bbox_to_anchor=(0.995, 0.86), fontsize=7.2, frameon=False, handlelength=0.9, handleheight=0.9, borderaxespad=0)

# ---- C: foreground vs background
bx = pd.read_csv("boxplot_data.csv"); w = pd.read_csv("wilcox.csv").iloc[0]
axC = fig.add_axes(inch(8.15, 0.82, 3.65, 2.58))
rng = np.random.default_rng(3)
groups = [("Background", bx[~bx.fg].log10_range.values, GREY), ("Foreground", bx[bx.fg].log10_range.values, GREEN)]
from scipy.stats import gaussian_kde
for i, (nm, v, c) in enumerate(groups):
    q1, med, q3 = np.percentile(v, [25, 50, 75]); iqr = q3 - q1
    lo, hi = v[v >= q1 - 1.5 * iqr].min(), v[v <= q3 + 1.5 * iqr].max()
    ys = np.linspace(v.min(), v.max(), 400); dens = gaussian_kde(v, bw_method=0.18)(ys); hw = 0.44 * dens / dens.max()
    dark = "#6B6B6B" if i == 0 else "#1E8A58"
    axC.fill_betweenx(ys, i - hw, i + hw, facecolor=c, alpha=0.30, edgecolor="none", zorder=2)       # translucent violin
    axC.plot(np.r_[i - hw, (i + hw)[::-1], (i - hw)[0]], np.r_[ys, ys[::-1], ys[0]], color=dark, lw=0.8, zorder=4, solid_joinstyle="round")
    axC.plot([i, i], [lo, hi], color=INK, lw=0.7, zorder=5)                                    # whisker range
    axC.add_patch(Rectangle((i - 0.035, q1), 0.07, iqr, fc=INK, ec=INK, lw=0.4, zorder=6))     # interquartile range
    axC.plot([i], [med], "o", ms=3.8, mfc="white", mec=INK, mew=0.6, zorder=7)                 # median
axC.set_xlim(-0.6, 1.6); axC.set_ylim(4.0, 8.55)
axC.set_xticks([0, 1]); axC.set_xticklabels([f"Background\n(n = {len(groups[0][1])})", f"Foreground: rapid-innovation lineages\n(n = {len(groups[1][1])})"], fontsize=7.6)
axC.set_ylabel(r"log$_{10}$(range size)", fontsize=8.5)
axC.spines[["top", "right"]].set_visible(False); axC.grid(axis="y", color="#EEEEEE", lw=0.4, zorder=0); axC.set_axisbelow(True)
p = w.p; ex = int(np.floor(np.log10(p))); mant = p / 10 ** ex
axC.plot([0, 0, 1, 1], [8.15, 8.25, 8.25, 8.15], color=INK, lw=0.6)
axC.text(0.5, 8.3, f"Wilcoxon rank-sum, p = {mant:.1f} × 10$^{{{ex}}}$", ha="center", va="bottom", fontsize=7.4, color=INK)

# ---- panel letters
for lab, (xx, yy) in {"A": (0.012, 0.992), "B": (0.012, 3.68 / fig_h), "C": (0.60, 3.68 / fig_h)}.items():
    fig.text(xx, yy, lab, fontsize=12, fontweight="bold", ha="left", va="top")

for ext in ("pdf", "png"):
    fig.savefig(f"final_figure_1.3.{ext}", dpi=300)
print("saved")
