from ete3 import Tree
import pandas as pd
from sys import argv


PATH_TO_PROJECT="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

# tree_path = f"{PATH_TO_PROJECT}/initial_data/trees_labeled/pruned_tree_TOGA_inno_marked.newick"
foreground = pd.read_csv(f"{PATH_TO_PROJECT}/initial_data/foreground_species.tsv", sep = "\s+")

# foreground = foreground[foreground["Inno"] == "high"]["Scientific_Name"].tolist()

tree_path = f"{PATH_TO_PROJECT}/initial_data/trees_labeled/pruned_tree_TOGA_small_selection_marked.newick"

# print(foreground)

background = argv[1].split(",")
gene_name = argv[2]
phenotype = argv[3]

foreground = foreground[foreground["Inno"] == phenotype]["Scientific_Name"].tolist()

all_leaves = []
for i in background:
    if i in foreground:
        i = i+"{TEST}"
    all_leaves.append(i)

# foreground_test = [i+"{TEST}" for i in foreground]

# all_leaves = foreground_test + background
print("Leaves to preserve:", all_leaves)
print("Total leaves to preserve:", len(all_leaves))

t = Tree(tree_path)
t.prune(
    all_leaves, 
    preserve_branch_length=True
)
t.write(outfile=f"{PATH_TO_PROJECT}/tmp/trees/pruned_tree_{gene_name}.newick")

# python /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/branch_site_selections/selection_tests/filter_tree.py \
#     'Falco_peregrinus,Egretta_garzetta,Agelaius_phoeniceus,Sturnus_vulgaris,Phalacrocorax_carbo,Sterna_hirundo,Parus_major,Haliaeetus_leucocephalus,Falco_tinnunculus,Passer_domesticus' \
#     'Rhinopomastus_cyanomelas,Upupa_epops,Baryphthengus_martii,Mesembrinibis_cayennensis,Psittacula_krameri,Guaruba_guarouba' \
#     'gene'