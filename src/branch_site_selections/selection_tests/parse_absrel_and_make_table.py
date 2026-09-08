import json
import glob
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests

prefix = "trimmed_prank_"
files = glob.glob(f"{prefix}absrel_small_results_inno/*.json")

rows = []

for f in files:
    print(f)
    gene = os.path.basename(f).replace(".json", "")

    if os.path.getsize(f) == 0:
        continue

    with open(f) as fh:
        try: 
            data = json.load(fh)
        except json.JSONDecodeError:
            print(f"Ошибка декодирования JSON в файле {f}, пропуск.")
            continue

    tested = data["tested"]["0"]
    test_branches = [b for b, t in tested.items() if t == "test"]

    branch_data = data["branch attributes"]["0"]

    for branch in test_branches:
        pval = branch_data[branch]["Corrected P-value"]
        rows.append({
            "gene": gene,
            "branch": branch,
            "raw_p": pval
        })

df = pd.DataFrame(rows)

# ---- FDR внутри каждой ВЕТВИ (вида) по всем генам ----
df["fdr"] = (
    df.groupby("branch")["raw_p"]
      .transform(lambda p: multipletests(p, method="fdr_bh")[1])
)

# Теперь фильтр: у каких генов ВСЕ тестируемые ветви значимы
significant_genes = (
    df.groupby("gene")
      .apply(lambda x: (x["fdr"] < 0.05).all())
)

sig_genes = significant_genes[significant_genes].index.tolist()

# Summary
summary = (
    df.assign(sig=df["fdr"] < 0.05)
      .groupby("gene")
      .agg(
          n_branches=("branch", "count"),
          n_sig=("sig", "sum"),
          sig_branches=("branch", lambda x: list(x[df.loc[x.index, "fdr"] < 0.05])),
          not_sig_branches=("branch", lambda x: list(x[df.loc[x.index, "fdr"] >= 0.05])),
      )
)

summary["fraction_sig"] = summary["n_sig"] / summary["n_branches"]
summary.sort_values(by="fraction_sig", ascending=False).to_csv(f"{prefix}absrel_e_summary_per_gene_filter.csv")

df_filtered = df[df["gene"].isin(sig_genes)]

print("Гены, где все ветви под FDR < 0.05:")
print(sig_genes)

df.to_csv(f"{prefix}absrel_full_table_filter.csv", index=False)
#df_filtered.to_csv("absrel_genes_all_branches_selected_filter.csv", index=False)

print("Готово.")