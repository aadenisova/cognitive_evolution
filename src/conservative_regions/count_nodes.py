from ete3 import Tree
import os

# dir_path = "rer_trees_birds_no_na_by_final_cnees_long/"

# for file in os.listdir(dir_path)[:10]:
    
#     t = Tree(open(f"{dir_path}/{file}", "r").read(), format=1)   
#     print(t.count_nodes())


path_to_file = "/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/conserved_regions/Sackton_et_al/rer_trees_birds_no_na_by_final_cnees_long/mCE1909781.tree"

t = Tree(open(f"{path_to_file}", "r").read().split("\t")[1], format=1)   
print(len(t.get_leaf_names()))
