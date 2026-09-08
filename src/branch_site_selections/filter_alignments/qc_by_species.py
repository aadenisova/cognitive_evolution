from pathlib import Path
from scipy import stats
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import pandas as pd

PATH_TO_PROJECT="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

FASTA_DIR = Path(f"{PATH_TO_PROJECT}/tmp/cds_TOGA_clean_fasta")
OUT = f"{PATH_TO_PROJECT}/tmp/species_qc_summary.tsv"

species_stats = defaultdict(list)

df_alignment = pd.read_csv(f"{PATH_TO_PROJECT}/tmp/alignment_qc_summary_filtered.tsv", sep = "\t")

# n = 0
for gene in df_alignment["gene"]:
    fasta = f"{PATH_TO_PROJECT}/tmp/cds_TOGA_clean_fasta/{gene}.fa"

    records = list(SeqIO.parse(fasta, "fasta"))
    aln_len = len(records[0].seq)

    for r in records:
        species = r.id #.split()[0].replace("vs_", "")
        gap_frac = r.seq.count("-") / aln_len
        species_stats[species].append(gap_frac)
    # n += 1
    # if n > 10:
    #     break

summary = []

for sp, gaps in species_stats.items():
    summary.append({
        "species": sp,
        "n_genes": len(gaps),
        "mean_gap_frac": np.mean(gaps),
        "median_gap_frac": np.median(gaps),
        "max_gap_frac": np.max(gaps),
    })

df = pd.DataFrame(summary)

foreground_df = pd.read_csv(f"{PATH_TO_PROJECT}/initial_data/foreground_species.tsv", sep = "\s+")

foreground_df = foreground_df[foreground_df["Inno"] == "high"]
print(f'Number of branches tested: {len(foreground_df["Scientific_Name"].unique())}')


foreground = set(foreground_df["Scientific_Name"].tolist())
df["class"] = df["species"].apply(
    lambda x: "foreground" if x in foreground else "background"
)

df.to_csv(OUT, sep="\t", index=False)

print("Saved species QC:", OUT)
print(df.groupby("class")[["mean_gap_frac", "median_gap_frac"]].mean())

res = stats.mannwhitneyu(
    df[df["class"] == "foreground"]["mean_gap_frac"],
    df[df["class"] == "background"]["mean_gap_frac"], 
    alternative="two-sided",
    method="exact",
)

print(res)