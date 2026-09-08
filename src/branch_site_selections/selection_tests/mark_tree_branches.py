from ete3 import Tree
import pandas as pd
import os
from sys import argv

PATH_TO_DATA = ""
PREFIX= argv[1]

# newick_path = f"{PATH_TO_DATA}/cooney_tree.tree"
# tsv_path = f"{PATH_TO_DATA}/selection_TOGA.tsv"
newick_file = argv[2]
gene = argv[3]

df = pd.read_csv(f"{PATH_TO_DATA}{PREFIX}absrel_e_summary_per_gene_filter.csv", sep=",")

branches = df[df["gene"] == gene]["sig_branches"].values[0]
print(branches)

t = Tree(newick_file, format=1)

for leaf in t.iter_leaves():
    if "{TEST}" in leaf.name:
        leaf.name = leaf.name.replace("{TEST}", "")

    if leaf.name in branches:
        leaf.name = leaf.name+"{TEST}"
print(t)

t.write(format=1, outfile=f"{PATH_TO_DATA}meme_trees/{gene}.newick")