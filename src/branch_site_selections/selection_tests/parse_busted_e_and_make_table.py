import json
import glob
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests

# files = glob.glob("absrel_results_correct_inno_10/*.json")
prefix = "trimmed_prank_"
files = glob.glob(f"{prefix}busted_e_small_results_inno/*.json")

rows = []

for f in files:
    gene = os.path.basename(f).replace(".json", "")


    if os.path.getsize(f) == 0:
        continue

    with open(f) as fh:
        data = json.load(fh)

    LRT_e = data["test results"]["LRT"]
    p_val_e = data["test results"]["p-value"]

    # fits = data['fits']["Unconstrained model"]['Rate Distributions']

    # rates_background = fits['Background']
    # rates_test = fits['Test']
    # rates_syn = fits["Synonymous site-to-site rates"]

    
    # error_sink_syn = sum(
    #     v["proportion"]
    #     for v in rates_syn.values()
    #     if v["rate"] > 5   # или >10, зависит от твоей калибровки
    # )
    
    # possel_test_prop = 0
    # possel_test_maxomega = 0
    # for v in rates_test.values():
    #     if v["omega"] > 1:
    #         possel_test_prop += v["proportion"]
    #         possel_test_maxomega = max(possel_test_maxomega, v["omega"])

    # omega0_test_prop = sum(
    #     v["proportion"]
    #     for v in rates_test.values()
    #     if v["omega"] == 0
    # )

    # omega_0_proportion = 0
    # for cat, values in rates_syn.items():
    #     omega = values.get('rate', 0)
    #     if omega == 0:  
    #         omega_0_proportion += values.get('proportion', 0)

    # error_sink_proportion_test = 0
    # for cat, values in rates_test.items():
    #     omega = values.get('omega', 0)
    #     if omega >= 100:  # error-sink
    #         error_sink_proportion_test += values.get('proportion', 0)

    rows.append({
            "gene": gene,
            # "p-value": p_val,
            # "LRT": LRT,
            "p-value_e": p_val_e,
            "LRT_e": LRT_e,
        })

df = pd.DataFrame(rows)
 
# ---- FDR внутри каждой ВЕТВИ (вида) по всем генам ----
# df["fdr"] = multipletests(df["p-value"], method="fdr_bh")[1]
df["fdr_e"] = multipletests(df["p-value_e"], method="fdr_bh")[1]
# df["fdr"] = df["p-value"].apply((lambda p: multipletests(p, method="fdr_bh")[1]))

df.sort_values(by = ["fdr_e"], ascending=True).to_csv(f"{prefix}busted_e_small_table.csv", index=False)
# df_filtered.to_csv("busted_e_genes_all_branches_selected.csv", index=False)

df[
    (df['fdr_e'] < 0.05) 
].to_csv(f"{prefix}busted_e_small_significant_genes.csv", index=False)