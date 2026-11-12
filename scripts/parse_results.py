import re
import sys

import matplotlib
matplotlib.use("Agg")          # draw plots without a screen
import matplotlib.pyplot as plt
import pandas as pd

# Send errors and messages to the Snakemake log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

rank_file = snakemake.input.rank
csv_file = snakemake.output.csv
plot_file = snakemake.output.plot
pair = snakemake.wildcards.pair
top_n = snakemake.params.top

# --- Read the LightDock ranking ---
rows = []
with open(rank_file) as f:
    columns = next(f).split()                       # first line = column names
    for line in f:
        if not line.strip():
            continue                                # skip empty lines
        coords = re.search(r"\(.*?\)", line).group(0)
        fields = line.replace(coords, "COORDS").split()
        fields[2] = coords                          # put the coordinates back
        rows.append(fields)

df = pd.DataFrame(rows, columns=columns)

# Text -> numbers
for col in ["Swarm", "Glowworm", "Clashes"]:
    df[col] = df[col].astype(int)
for col in ["Luciferin", "Scoring"]:
    df[col] = df[col].astype(float)

# Rank (file is already sorted, best first) and full model path
df.insert(0, "Rank", range(1, len(df) + 1))
df["Model"] = "swarm_" + df["Swarm"].astype(str) + "/" + df["PDB"]

# --- Save the table ---
df[["Rank", "Swarm", "Glowworm", "Model", "Clashes", "Luciferin", "Scoring"]].to_csv(csv_file, index=False)
print(f"{pair}: {len(df)} poses, best score {df['Scoring'].iloc[0]:.3f}")

# --- Plot the 20 best scores ---
top = df.head(20)
colors = ["tab:green" if r <= top_n else "tab:gray" for r in top["Rank"]]
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(top["Rank"], top["Scoring"], color=colors)
ax.set_xlabel("Rank")
ax.set_ylabel("LightDock score (higher = better)")
ax.set_title(f"{pair}: 20 best docking poses (green = top {top_n} models)")
ax.set_xticks(top["Rank"])
fig.tight_layout()
fig.savefig(plot_file, dpi=150)