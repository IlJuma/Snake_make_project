"""
Parse the LightDock ranking of one protein pair.
Making the outcome from the parsing readable by converting it into a .csv file
This script is also shwing a visualization of the best scores and theire height

Called by the Snakemake rule `parse_results` (one run per pair).

Input:
    rank_by_scoring.list - all clustered docking poses, sorted by score (best first)
Output:
    scores.csv     - table of all poses: rank, swarm, glowworm, model file, clashes, score
    top_scores.png - bar chart of the 20 best scores, top N models highlighted
"""

import re
import sys

import matplotlib
matplotlib.use("Agg")          # draw plots without a screen
import matplotlib.pyplot as plt
import pandas as pd

# Send errors and messages to the Snakemake log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

#geting the paths and file names from the snakemakefile to retrife the right files
#the snakemake object already holds this information (is created for each run 
#by the snakemakefile)

rank_file = snakemake.input.rank
csv_file = snakemake.output.csv
plot_file = snakemake.output.plot
pair = snakemake.wildcards.pair
#number from config file how many top models should be used
top_n = snakemake.params.top

# --- Read the LightDock ranking ---
rows = []
with open(rank_file) as f:
    columns = next(f).split()                       # first line = column names
    for line in f:
        if not line.strip():
            continue                                # skip empty lines
        #Finding the coordinates by the searching for values in brackets
        coords = re.search(r"\(.*?\)", line).group(0)
        #split the coordinates in seperate columns
        fields = line.replace(coords, "COORDS").split()
        fields[2] = coords                          # put the coordinates back
        rows.append(fields)

df = pd.DataFrame(rows, columns=columns)

# Text -> numbers
#Values from the file are strings, here they are converted back to numbers to use
#them later for plotting and sorting
for col in ["Swarm", "Glowworm", "Clashes"]:
    df[col] = df[col].astype(int)
for col in ["Luciferin", "Scoring"]:
    df[col] = df[col].astype(float)

# Rank (file is already sorted, best first) and full model path
df.insert(0, "Rank", range(1, len(df) + 1))
#Every swarm has files needed in the folder, therfore it is needed to identify the 
#model file
df["Model"] = "swarm_" + df["Swarm"].astype(str) + "/" + df["PDB"]

# --- Save the table ---
#Only keeping the important columns (the coordinates for example are not needed)
df[["Rank", "Swarm", "Glowworm", "Model", "Clashes", "Luciferin", "Scoring"]].to_csv(csv_file, index=False)
#Writes a summary in the log file
print(f"{pair}: {len(df)} poses, best score {df['Scoring'].iloc[0]:.3f}")

# --- Plot the 20 best scores ---
#The top scores (defined in config) are green, the rest is grey in total best 20
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