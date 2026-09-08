#!/usr/bin/env python3
"""
06_generate_excel_v3.py
Generate final Excel with 12 sheets for v3.1:
- GO split by ontology (BP/CC/MF) per organism
- Gene symbols in Robust_CNEEs and GO leadingEdge
- DAR enrichment without 'all' direction (no NaN)
- All sheets sorted by significance
"""
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
from openpyxl.utils.dataframe import dataframe_to_rows
import subprocess

BASE = "/mnt/results/rerconverge_cont_v2"
OUTPUT_WORKSPACE = "/workspace/rerconverge_log_noexon_summary_v3.xlsx"
OUTPUT_FINAL = "/mnt/results/rerconverge_log_noexon_summary_v3.xlsx"

GREEN_FILL = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
YELLOW_FILL = PatternFill(start_color="FFEB9C", end_color="FFEB9C", fill_type="solid")
HEADER_FILL = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
HEADER_FONT = Font(color="FFFFFF", bold=True, size=11)
THIN_BORDER = Border(
    left=Side(style="thin"), right=Side(style="thin"),
    top=Side(style="thin"), bottom=Side(style="thin")
)

def write_df_to_sheet(ws, df, start_row=1):
    for r_idx, row in enumerate(dataframe_to_rows(df, index=False, header=True), start_row):
        for c_idx, val in enumerate(row, 1):
            cell = ws.cell(row=r_idx, column=c_idx, value=val)
            if r_idx == start_row:
                cell.fill = HEADER_FILL
                cell.font = HEADER_FONT
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
                cell.border = THIN_BORDER
    for i, col in enumerate(df.columns, 1):
        max_len = max(
            len(str(col)),
            df[col].astype(str).str.len().max() if len(df) > 0 else 0
        )
        if i <= 26:
            col_letter = chr(64 + i)
        else:
            col_letter = chr(64 + (i - 1) // 26) + chr(65 + (i - 1) % 26)
        ws.column_dimensions[col_letter].width = min(max_len + 2, 35)
    ws.freeze_panes = "A2"
    return start_row + len(df) + 1

def highlight_padj(ws, df, padj_col_name, fill_color=GREEN_FILL, threshold=0.05):
    """Highlight rows where padj < threshold."""
    if padj_col_name not in df.columns:
        return
    padj_col = list(df.columns).index(padj_col_name) + 1
    n_cols = len(df.columns)
    for r_idx in range(2, len(df) + 2):
        val = ws.cell(row=r_idx, column=padj_col).value
        if val is not None and not pd.isna(val) and val < threshold:
            for c in range(1, n_cols + 1):
                ws.cell(row=r_idx, column=c).fill = fill_color

def main():
    wb = Workbook()

    # === Sheet 1: Robust_CNEEs ===
    print("Building Robust_CNEEs sheet...")
    robust = pd.read_csv(f"{BASE}/robust_cnees_annotated_v3.csv")
    # Add geneSymbol from comprehensive gene_symbol_map (covers ALL geneIds incl. Distal Intergenic)
    sym_map = pd.read_csv(f"{BASE}/gene_symbol_map.csv")
    id_to_sym = dict(zip(sym_map["geneId"].astype(str), sym_map["symbol"]))
    robust["geneSymbol"] = robust["geneId"].apply(
        lambda x: id_to_sym.get(str(int(x)) if pd.notna(x) else "", "") if pd.notna(x) else ""
    )
    robust = robust.sort_values("permpval")
    cols = ["CNEE", "Rho", "permpval", "P", "observed_stat", "sign_stability",
            "direction", "attributed", "category", "geneId", "geneSymbol",
            "distanceToTSS", "annotation"]
    robust = robust[[c for c in cols if c in robust.columns]]
    ws1 = wb.active
    ws1.title = "Robust_CNEEs"
    write_df_to_sheet(ws1, robust)

    # === Sheet 2: Gene_Stats ===
    print("Building Gene_Stats sheet...")
    gene = pd.read_csv(f"{BASE}/gene_level_stats_v3.csv")
    # Update geneSymbol from comprehensive map (ensures all symbols are current)
    sym_map = pd.read_csv(f"{BASE}/gene_symbol_map.csv")
    id_to_sym_gs = dict(zip(sym_map["geneId"].astype(str), sym_map["symbol"]))
    gene["geneSymbol"] = gene["geneId"].apply(
        lambda x: id_to_sym_gs.get(str(int(x)) if pd.notna(x) else "", "") if pd.notna(x) else ""
    )
    # Sort by most significant binom p-value (min of pos/neg)
    gene["_min_p"] = gene[["binom_p_pos", "binom_p_neg"]].min(axis=1)
    gene = gene.sort_values("_min_p")
    gene = gene.drop(columns=["_min_p"])
    ws2 = wb.create_sheet("Gene_Stats")
    write_df_to_sheet(ws2, gene)

    # === Sheets 3-8: GO by organism x ontology ===
    go_sheets = [
        ("GO_Chicken_BP", "chicken_BP"),
        ("GO_Chicken_CC", "chicken_CC"),
        ("GO_Chicken_MF", "chicken_MF"),
        ("GO_Human_BP",   "human_BP"),
        ("GO_Human_CC",   "human_CC"),
        ("GO_Human_MF",   "human_MF"),
    ]
    for sheet_name, file_key in go_sheets:
        print(f"Building {sheet_name} sheet...")
        f = f"{BASE}/phase_d_gsea/gsea_v3_{file_key}.csv"
        df = pd.read_csv(f)
        df = df.sort_values("padj")
        ws = wb.create_sheet(sheet_name)
        write_df_to_sheet(ws, df)
        highlight_padj(ws, df, "padj")

    # === Sheet 9: DAR_Enrichment ===
    print("Building DAR_Enrichment sheet...")
    dar = pd.read_csv(f"{BASE}/dar_enrichment_23annotated_v3.csv")
    # Sort: negative first (by padj), then positive (by padj)
    dar["_dir_order"] = dar["direction"].map({"negative": 0, "positive": 1})
    dar = dar.sort_values(["_dir_order", "padj_wilcoxon_23"])
    dar = dar.drop(columns=["_dir_order"])
    ws9 = wb.create_sheet("DAR_Enrichment")
    write_df_to_sheet(ws9, dar)
    # Highlight: green = FDR<0.05 + effect>=0.015, yellow = FDR<0.05 + effect<0.015
    padj_col = list(dar.columns).index("padj_wilcoxon_23") + 1
    effect_col = list(dar.columns).index("effect_size") + 1
    n_cols = len(dar.columns)
    for r_idx in range(2, len(dar) + 2):
        fdr_val = ws9.cell(row=r_idx, column=padj_col).value
        eff_val = ws9.cell(row=r_idx, column=effect_col).value
        if fdr_val is not None and not pd.isna(fdr_val) and fdr_val < 0.05:
            if eff_val is not None and not pd.isna(eff_val) and eff_val >= 0.015:
                fill = GREEN_FILL
            else:
                fill = YELLOW_FILL
            for c in range(1, n_cols + 1):
                ws9.cell(row=r_idx, column=c).fill = fill

    # === Sheet 10: Element_Enrichment ===
    print("Building Element_Enrichment sheet...")
    ws10 = wb.create_sheet("Element_Enrichment")
    annot = pd.read_csv(f"{BASE}/phase_d_gsea/cnee_annotation_cont_noexon.csv",
                        usecols=["name", "annotation"])
    annot["category"] = annot["annotation"].str.split("(").str[0].str.strip()
    genome_counts = annot["category"].value_counts()
    n_genome = len(annot)
    robust_annot = pd.read_csv(f"{BASE}/robust_cnees_annotated_v3.csv")
    freq_rows = []
    for cat in ["Distal Intergenic", "Intron", "Promoter", "Downstream", "3' UTR", "5' UTR"]:
        g_count = int(genome_counts.get(cat, 0))
        g_prop = g_count / n_genome
        for direction in ["all", "positive", "negative"]:
            if direction == "all":
                sub = robust_annot
            else:
                sub = robust_annot[robust_annot["direction"] == direction]
            r_count = int(sub["category"].eq(cat).sum())
            r_prop = r_count / len(sub) if len(sub) > 0 else 0
            ratio = r_prop / g_prop if g_prop > 0 else np.nan
            freq_rows.append({
                "category": cat, "direction": direction,
                "robust_count": r_count, "robust_prop": f"{r_prop:.4f}",
                "genome_count": g_count, "genome_prop": f"{g_prop:.4f}",
                "ratio": f"{ratio:.2f}"
            })
    freq_df = pd.DataFrame(freq_rows)
    next_row = write_df_to_sheet(ws10, freq_df)
    next_row += 1
    fisher = pd.read_csv(f"{BASE}/phase_e_elements/element_fisher_summary_v2.csv")
    fisher["_dir_order"] = fisher["direction"].map({"negative": 0, "positive": 1, "all": 2})
    fisher = fisher.sort_values(["_dir_order", "padj"])
    fisher = fisher.drop(columns=["_dir_order"])
    for r_idx, row in enumerate(dataframe_to_rows(fisher, index=False, header=True), next_row):
        for c_idx, val in enumerate(row, 1):
            cell = ws10.cell(row=r_idx, column=c_idx, value=val)
            if r_idx == next_row:
                cell.fill = HEADER_FILL
                cell.font = HEADER_FONT
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
                cell.border = THIN_BORDER
    ws10.freeze_panes = "A2"

    # === Sheet 11: LeaveOut_Stability ===
    print("Building LeaveOut_Stability sheet...")
    ws11 = wb.create_sheet("LeaveOut_Stability")
    leaveout = pd.read_csv(f"{BASE}/phase_c_leave_out/leave_out_species_results.csv")
    summary_rows = [
        {"metric": "Spearman rho (min)", "value": f"{leaveout['spearman_rho'].min():.4f}"},
        {"metric": "Spearman rho (median)", "value": f"{leaveout['spearman_rho'].median():.4f}"},
        {"metric": "Spearman rho (max)", "value": f"{leaveout['spearman_rho'].max():.4f}"},
        {"metric": "Top100 overlap (min)", "value": f"{leaveout['top100_overlap'].min():.3f}"},
        {"metric": "Top100 overlap (median)", "value": f"{leaveout['top100_overlap'].median():.3f}"},
        {"metric": "Top1000 overlap (min)", "value": f"{leaveout['top1000_overlap'].min():.3f}"},
        {"metric": "Top1000 overlap (median)", "value": f"{leaveout['top1000_overlap'].median():.3f}"},
        {"metric": "Sign stability (all 7,350 robust CNEEs)", "value": "100% (all sign_stability=1.0)"},
        {"metric": "N leave-out iterations", "value": "92 (one per species)"},
    ]
    summary_df = pd.DataFrame(summary_rows)
    next_row = write_df_to_sheet(ws11, summary_df)
    next_row += 1
    leaveout = leaveout.sort_values("spearman_rho")
    for r_idx, row in enumerate(dataframe_to_rows(leaveout, index=False, header=True), next_row):
        for c_idx, val in enumerate(row, 1):
            cell = ws11.cell(row=r_idx, column=c_idx, value=val)
            if r_idx == next_row:
                cell.fill = HEADER_FILL
                cell.font = HEADER_FONT
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
                cell.border = THIN_BORDER

    # === Sheet 12: Methods ===
    print("Building Methods sheet...")
    ws12 = wb.create_sheet("Methods")
    methods_text = [
        ("RERconverge: непрерывный анализ (log, noexon, permutation FDR) — v3.2", True),
        ("", False),
        ("Данные: 272,811 CNEE, 92 вида птиц, фенотип TECHINNO2025_ResEff", False),
        ("Трансформация: log (лучшая калибровка, lambda GC=19.59, выбрана в Phase A)", False),
        ("Удаление экзонов: 84 CNEE удалены, осталось 272,730", False),
        ("Пермутации: 2,000, permpval для всех CNEE", False),
        ("FDR: BH на permpval (0 значимых), global FDR (0 значимых)", False),
        ("Robust CNEE: permpval<0.05 + sign stability>=80% -> 7,350 (все 100% stable)", False),
        ("  Positive Rho: 3,804 / Negative Rho: 3,546", False),
        ("Leave-out: 92 итерации, Spearman 0.968-0.9999", False),
        ("", False),
        ("Атрибуция CNEE к генам (v3):", True),
        ("  Только ген + близкая регуляторная область:", False),
        ("    Intron + 5'UTR + 3'UTR + Promoter (<=3kb от TSS) + Downstream (<=300bp)", False),
        ("    Distal Intergenic ИСКЛЮЧЕН из ген-уровневой статистики", False),
        ("  Атрибутировано: 43,572 CNEE -> 3,878 генов", False),
        ("  Robust атрибутировано: 1,141 CNEE -> 677 генов", False),
        ("  Background rate: 2.62%", False),
        ("  Max k_pos=5, Max k_neg=15, Median n_total=3", False),
        ("", False),
        ("Маппинг символов генов (v3.2 — исправлено):", True),
        ("  org.Gg.eg.db для ВСЕХ geneId (4,443), включая Distal Intergenic", False),
        ("  Ранее: только 3,878 атрибутированных генов -> пропуски (напр. 428134/TRPC4AP)", False),
        ("  Теперь: 4,440 из 4,443 замапплено (99.9%), 3 unmapped (770433, 102465680, 102466994)", False),
        ("  geneSymbol добавлен в Robust_CNEEs и обновлён в Gene_Stats", False),
        ("", False),
        ("observed_stat в Robust_CNEEs:", True),
        ("  t-статистика Pearson корреляции между RER и фенотипом (из RERconverge)", False),
        ("  Коррелирует с Rho и P: чем больше |t|, тем сильнее связь", False),
        ("", False),
        ("GSEA (v3.1 — онтологии разделены):", True),
        ("  GO термы разделены по онтологиям: BP, CC, MF — отдельные таблицы", False),
        ("  Методы: count-based + binomial, fgsea scoreType='pos', BH FDR", False),
        ("  Организмы: chicken (org.Gg.eg.db) + human orthologs (homologene)", False),
        ("  leadingEdge содержит как Entrez ID, так и символы генов", False),
        ("  Chicken: BP=440, CC=212, MF=300 термов, 0 значимых", False),
        ("  Human: BP=3132, CC=808, MF=816 термов, 1 значимый", False),
        ("    GO:0001669 (acrosomal vesicle, CC) — NES=1.235, padj=0.047", False),
        ("  Меньше мощность чем v2 (3,878 генов vs 5,554, 438 non-zero vs 1,775)", False),
        ("  fgsea стохастичен: p-values могут слегка меняться между запусками", False),
        ("", False),
        ("Element enrichment: Fisher exact, 6 категорий x 3 направления", False),
        ("  0 значимых (FDR<0.05). Intron neg: OR=1.13, p=0.009, FDR=0.054", False),
        ("", False),
        ("DAR enrichment (v3.1 — без omnibus 'all'):", True),
        ("  Перекрытие: >50% длины CNEE должно быть внутри DAR", False),
        ("  Направления: только negative и positive Rho (omnibus 'all' убран)", False),
        ("  Метрика: -log10(permpval), Wilcoxon rank-sum + analytical permutation z-score", False),
        ("  23 annotated cell types, BH FDR на 23 тестов per direction", False),
        ("  Результат: 7 значимых (negative Rho, все green effect>=0.015)", False),
        ("    OL, AST, EXC_GLU-8, OPC, INH_VIP, IOL, EXC_GLU-5", False),
        ("  Positive: 0 значимых", False),
        ("  Цвет: зеленый = FDR<0.05 + effect>=0.015, желтый = FDR<0.05 + effect<0.015", False),
        ("", False),
        ("Итог v3.2: 0 геномно-значимых CNEE, 1 GO-терм (human CC), 7 значимых типов клеток", False),
        ("Сигнал коллективный, атрибуция строгая (только ген+промоутер)", False),
        ("Все таблицы отсортированы по значимости (наиболее значимые находки вверху)", False),
        ("", False),
        ("Фигуры (v3.2):", True),
        ("  A: Филогенетическое дерево (rectangular) + innovation values", False),
        ("  B1: Top 10 intron CNEEs — Rho, distance category, p-value", False),
        ("  B2: Promoter <=2kb vs >2kb — |Rho| comparison (Wilcoxon)", False),
        ("  D: DAR enrichment — effect sizes for 23 cell types", False),
        ("  S: Element enrichment distribution (supplement)", False),
        ("  Combined figure assembled with cowplot", False),
        ("  Дерево экспортировано: species_tree_rooted.nwk (Latin names, rooted)", False),
        ("  Innovation table: innovation_table.csv (92 species)", False),
    ]
    for r_idx, (text, is_header) in enumerate(methods_text, 1):
        cell = ws12.cell(row=r_idx, column=1, value=text)
        if is_header:
            cell.font = Font(bold=True, size=14)
        cell.alignment = Alignment(horizontal="left", vertical="top")
    ws12.column_dimensions["A"].width = 95

    # Save
    print(f"\nSaving to {OUTPUT_WORKSPACE}...")
    wb.save(OUTPUT_WORKSPACE)
    subprocess.run(["cp", OUTPUT_WORKSPACE, OUTPUT_FINAL], check=True)
    print(f"Copied to {OUTPUT_FINAL}")

    import os
    size = os.path.getsize(OUTPUT_FINAL)
    print(f"File size: {size:,} bytes")

    from openpyxl import load_workbook
    wb2 = load_workbook(OUTPUT_FINAL, read_only=True)
    print(f"Sheets: {wb2.sheetnames}")
    for name in wb2.sheetnames:
        ws = wb2[name]
        print(f"  {name}: {ws.max_row} rows x {ws.max_column} cols")
    wb2.close()

if __name__ == "__main__":
    main()
