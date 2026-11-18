"""Collect the RMSD tables of all pairs into one overview table.

One row per pair: how good the top-ranked pose is, how good the best of the
top N poses is, and how many of them are at least CAPRI 'acceptable'.
"""

import sys

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

GOOD = ("acceptable", "medium", "high")

rows = []
for pair, csv in zip(snakemake.params.pairs, snakemake.input):
    df = pd.read_csv(csv)
    top1 = df.iloc[0]
    best = df.loc[df["DockQ"].idxmax()]            # best pose among the top N
    rows.append({
        "Pair": pair,
        "Models": len(df),
        "Top1_L_RMSD": top1["L_RMSD"],
        "Top1_DockQ": top1["DockQ"],
        "Top1_CAPRI": top1["CAPRI"],
        "Best_rank": int(best["Rank"]),
        "Best_L_RMSD": best["L_RMSD"],
        "Best_DockQ": best["DockQ"],
        "Best_CAPRI": best["CAPRI"],
        "Acceptable_or_better": int(df["CAPRI"].isin(GOOD).sum()),
    })

summary = pd.DataFrame(rows)
summary.to_csv(snakemake.output[0], index=False)
print(summary.to_string(index=False))
