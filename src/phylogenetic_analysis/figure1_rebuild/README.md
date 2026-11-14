# Figure 1 rebuild (tree + phylogenetic signal + foreground/background)

Python/matplotlib version of `final_figure_1.2`: circular tree with order silhouettes, VGP marker, time scale (A),
Blomberg's K / Pagel's lambda bars (B), foreground vs background range size violins (C).

## Inputs (put in `data/`, not tracked here)
`tree_clootl_overlap_25.nwk`, `inno_db_25_clootl_overlap.csv`, `tree_vgp_0_vgp.nwk`, `VGP_tag.png`

Tracked here: `signal_50.csv` (K / lambda, ResEff50, from Table 1), `species_50.csv` (830 species of the ResEff50 set).

## Run order
```
python3 download_silhouettes.py data/inno_db_25_clootl_overlap.csv -o sil   # PhyloPic SVGs + attribution.csv (needs internet)
python3 make_sil_paths.py            # sil/*.svg -> sil_paths.npy
Rscript 03_foreground.R              # boxplot_data.csv, wilcox.csv, foreground_tips.csv
Rscript 04_export_tree.R             # tree_edges.csv, tree_tips.csv, tree_depth.csv
python3 build_figure.py              # final_figure_1.3.pdf / .png
```
`sil/swift.png` (Caprimulgiformes) is a custom silhouette; the other credits are in `silhouette_credits.txt`.

Requires: Python (numpy, pandas, scipy, matplotlib, Pillow, svgpath2mpl); R (ape, phytools).
Time scale follows Stiller et al. 2024, Nature 629, 851-860 (Fig. 1b).
