from ete3 import Tree
import pandas as pd
from sys import argv


PATH_TO_PROJECT="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"
foreground = pd.read_csv(f"{PATH_TO_PROJECT}/initial_data/foreground_species_with_filtration.tsv", sep = "\t")

tree_path = argv[1]
gene_name = argv[2]
phenotype = argv[3]
tree_dir = argv[4]

species_to_mark = foreground[foreground["inno_selection"] == phenotype]["sci_name_2025"].tolist()
t = Tree(tree_path, format = 1)

is_foreground = False
for leaf in t.get_leaves():
    if leaf.name in species_to_mark:
        is_foreground = True
        leaf.name = leaf.name+"{TEST}"

t.write(outfile=f"{tree_dir}/{gene_name}.newick")

# print(f"n = {len(t.get_leaves())}")
# if is_foreground:
#     print(f"Tree for {gene_name} has foreground species marked with {{TEST}}")
# python /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/branch_site_selections/selection_tests/filter_tree.py \
#     'Falco_peregrinus,Egretta_garzetta,Agelaius_phoeniceus,Sturnus_vulgaris,Phalacrocorax_carbo,Sterna_hirundo,Parus_major,Haliaeetus_leucocephalus,Falco_tinnunculus,Passer_domesticus' \
#     'Rhinopomastus_cyanomelas,Upupa_epops,Baryphthengus_martii,Mesembrinibis_cayennensis,Psittacula_krameri,Guaruba_guarouba' \
#     'gene'