"""More sampling -> higher scores. Two panels."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

root = Path("analysis_all_poses")
LEVELS = {5: root / "sampling_5swarm", 50: root / "sampling_50swarm", 100: root}
PAIRS = ["1efn", "2g45", "2goo", "3fap"]


def load(p, n):
    f = LEVELS[n] / f"{p}_all_poses_rmsd.csv"
    return pd.read_csv(f) if f.exists() else None


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))

# ---- Panel A: score distribution for 1efn at three sampling depths ----
levels = [5, 50, 100]
data, labels, bests = [], [], []
for n in levels:
    d = load("1efn", n)
    data.append(d["Score"].values)
    labels.append(f"{n} swarms\n({len(d)} poses)")
    bests.append(d["Score"].max())

YLO, YHI = -48, max(bests) + 9

bp = ax1.boxplot(data, patch_artist=True, widths=0.55, showfliers=False,
                 whis=(5, 95), medianprops=dict(color="#1a1a1a", linewidth=2))
for patch in bp["boxes"]:
    patch.set_facecolor("#cfe0f2")
    patch.set_edgecolor("#4b7bab")

ax1.plot(range(1, len(levels) + 1), bests, "o-", color="tab:red",
         markersize=9, linewidth=2, zorder=5, label="best score found")
for i, b in enumerate(bests, start=1):
    ax1.annotate(f"{b:.1f}", (i, b), textcoords="offset points",
                 xytext=(16, -3), ha="left", fontsize=11,
                 color="tab:red", fontweight="bold")

for i, (n, arr) in enumerate(zip(levels, data), start=1):
    med = np.median(arr)
    ax1.annotate(f"median {med:.1f}", (i, med), textcoords="offset points",
                 xytext=(0, -20), ha="center", fontsize=9.5, color="#1a1a1a")
    if arr.min() < YLO:
        ax1.annotate(f"tail to {arr.min():.0f}", (i, YLO), textcoords="offset points",
                     xytext=(0, 10), ha="center", fontsize=9, color="0.4",
                     arrowprops=dict(arrowstyle="-|>", color="0.5", lw=1))

ax1.set_ylim(YLO, YHI)
ax1.set_xticklabels(labels, fontsize=10.5)
ax1.set_ylabel("LightDock score      higher = better", fontsize=11)
ax1.set_title("1efn: the score distribution tightens and rises\n"
              "box = middle 50% of poses, axis cut at -48", fontsize=12.5)
ax1.legend(fontsize=10, loc="upper left")
ax1.grid(axis="y", alpha=0.25)

# ---- Panel B: best score, 5 vs 100 swarms, all four pairs ----
colors = {"1efn": "tab:blue", "2g45": "tab:orange",
          "2goo": "tab:green", "3fap": "tab:purple"}
for p in PAIRS:
    lo = load(p, 5)["Score"].max()
    hi = load(p, 100)["Score"].max()
    ax2.plot([0, 1], [lo, hi], "o-", color=colors[p], linewidth=2.4,
             markersize=9, zorder=4)
    ax2.annotate(f"{p}   +{hi - lo:.1f}", (1, hi), textcoords="offset points",
                 xytext=(12, 0), va="center", fontsize=11,
                 color=colors[p], fontweight="bold")
    ax2.annotate(f"{lo:.1f}", (0, lo), textcoords="offset points",
                 xytext=(-12, 0), va="center", ha="right",
                 fontsize=10, color=colors[p])

ax2.set_xlim(-0.35, 1.52)
ax2.set_xticks([0, 1])
ax2.set_xticklabels(["5 swarms", "100 swarms"], fontsize=11.5)
ax2.set_ylabel("best score found", fontsize=11)
ax2.set_title("Best score improves for every pair\n"
              "20x more sampling, higher optimum each time", fontsize=12.5)
ax2.grid(axis="y", alpha=0.25)

fig.suptitle("More sampling finds higher-scoring poses", fontsize=15)
fig.subplots_adjust(top=0.80, bottom=0.12, left=0.07, right=0.97, wspace=0.26)
out = Path("analysis_all_poses") / "sampling_vs_score.png"
fig.savefig(out, dpi=200)
print("written:", out.resolve())
