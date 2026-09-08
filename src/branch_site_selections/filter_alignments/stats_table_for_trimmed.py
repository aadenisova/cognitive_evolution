from pathlib import Path
import subprocess
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import numpy as np

SMALL_DATA = True

PATH_TO_PROJECT="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

CDS_DIR = Path(f"{PATH_TO_PROJECT}/tmp/trimmed_small_prank_fa")
# OUT_FASTA_DIR = Path(f"{PATH_TO_PROJECT}/tmp/cds_small_TOGA_clean_fasta")
# OUT_PHYLIP_DIR = Path(f"{PATH_TO_PROJECT}/tmp/cds_small_TOGA_clean_phylip")

SPECIES_TABLE = f"{PATH_TO_PROJECT}/initial_data/toga_with_innovation_scores_in_cooney.csv"

SMALL_SELECTION = f"{PATH_TO_PROJECT}/initial_data/small_selection_toga.tsv"
SUMMARY_OUT = f"{PATH_TO_PROJECT}/tmp/trimmed_alignment_qc_summary.tsv"

species_df = pd.read_csv(SPECIES_TABLE, sep=",")

if SMALL_DATA:
    small_df = pd.read_csv(SMALL_SELECTION, sep="\t")
    species_df = species_df[species_df["Scientific_Name"].apply(lambda x: x in small_df["Scientific_Name"].values)]

# allowed_species = set(species_df["Assembly name", "Scientific_Name"])
species_df.index = species_df["Assembly name"]
allowed_species = species_df["Scientific_Name"].to_dict()

foreground = set(pd.read_csv(f"{PATH_TO_PROJECT}/initial_data/foreground_species.tsv", sep = "\s+")["Scientific_Name"].tolist())

def is_cds_valid(seq: str):
    """Check CDS length and stop codons"""
    seq = seq.replace("-", "").upper()
    length_ok = (len(seq) % 3 == 0)
    stop_codons = ["TAA", "TAG", "TGA"]
    stops = sum(seq[i:i+3] in stop_codons for i in range(0, len(seq)-2, 3))
    return length_ok, stops

def has_error(seq):
    return "!" in seq

summary = []

foreground_sp_in_alignment = set()

# n = 0
for aln_file in CDS_DIR.glob("*.fa*"):
    # if n > 10:
    #     break
    # n+=1
    gene = aln_file.stem
    records = list(SeqIO.parse(aln_file, "fasta"))

    species_map = defaultdict(list)

    for rec in records:
        sp = rec.id
        species_map[sp].append(rec)
        if sp in foreground:
            foreground_sp_in_alignment.add(sp)

    kept_records = []
    duplicate_species = []

    for sp, recs in species_map.items():
        if len(recs) == 1:
            best = recs[0]

        else:
            duplicate_species.append(sp)
            best = min(recs, key=lambda r: r.seq.count("-"))
        
        if not has_error(best):
            best.id = sp
            best.description = ""
            kept_records.append(best)
            

    if len(kept_records) < 4:
        summary.append({
            "gene": gene,
            "status": "discarded_too_few_species",
        })
        continue

    aln_len = len(kept_records[0].seq)

    # matrix of alignment
    aln_array = np.array([list(str(r.seq)) for r in kept_records])

    # gap columns
    gap_cols = np.all(aln_array == "-", axis=0)
    gap_cols_frac = gap_cols.mean()

    # columns with >50% gaps
    gap50_cols = (np.mean(aln_array == "-", axis=0) > 0.5).mean()

    # per-sequence stats
    gap_fracs = []
    n_fracs = []
    cds_ok_count = 0
    stop_codons_total = 0


    for r in kept_records:
        seq = str(r.seq)
        gap_fracs.append(seq.count("-") / aln_len)
        n_fracs.append(seq.upper().count("N") / aln_len)
        length_ok, stops = is_cds_valid(seq)
        cds_ok_count += int(length_ok)
        stop_codons_total += stops

    cds_ok_frac = cds_ok_count / len(kept_records)

    summary.append({
        "gene": gene,
        "status": "kept",
        "n_species": len(kept_records),
        "n_duplicates": len(duplicate_species),
        "duplicates": ",".join(duplicate_species),
        "alignment_length": aln_len,
        "gap_columns_frac": gap_cols_frac,
        "gap50_columns_frac": gap50_cols,
        "mean_gap_frac": np.mean(gap_fracs),
        "max_gap_frac": np.max(gap_fracs),
        "mean_N_frac": np.mean(n_fracs),
        "cds_ok_frac": cds_ok_frac,
        "stop_codons_total": stop_codons_total,
        "number_of_foreground_species": len(foreground_sp_in_alignment),
        "foreground_species": ','.join(list(foreground_sp_in_alignment)),
        "all_species": ",".join(list(species_map.keys())),
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(SUMMARY_OUT, sep="\t", index=False)

print("Done")
# print(summary_df.value_counts())
