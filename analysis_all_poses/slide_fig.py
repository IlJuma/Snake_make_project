"""Presentation figure: score vs L-RMSD for two contrasting pairs."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

base = Path(__file__).parent
s100 = base / "s100"

PANELS = [
    ("1efn", "HIV-1 Nef / Fyn SH3", "best outcome"),
    ("2goo", "BMP-2 / ActRII", "worst ranking"),
]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))

for ax, (p, name, tag) in zip(axes, PANELS):
    df = pd.read_csv(s100 / f"{p}_all_poses_rmsd.csv")
    top5 = df.nsmallest(5, "Rank")
    good = df[df["CAPRI"] != "incorrect"].sort_values("L_RMSD")
    bad = df[df["CAPRI"] == "incorrect"]
    r = df["Score"].corr(df["L_RMSD"])

    ax.scatter(bad["L_RMSD"], bad["Score"], s=20, color="0.72",
               linewidths=0, zorder=3, label=f"incorrect (n={len(bad)})")
    ax.scatter(good["L_RMSD"], good["Score"], s=95, color="tab:orange",
               edgecolors="#7a3d00", linewidths=1.0, zorder=6,
               label=f"acceptable (n={len(good)})")
    ax.scatter(top5["L_RMSD"], top5["Score"], facecolors="none",
               edgecolors="tab:red", s=150, linewidths=2, zorder=7,
               label="top 5 (delivered)")

    ax.axvline(10, color="k", linestyle="--", linewidth=1.1, zorder=1)

    # never clip an acceptable pose off the axis
    ylo = min(df["Score"].quantile(0.05), good["Score"].min() - 7)
    yhi = df["Score"].max() + 3
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(left=0)
    ax.text(10.8, yhi - 0.03 * (yhi - ylo), "acceptable threshold",
            fontsize=9.5, color="0.25", va="top")

    mid = (ylo + yhi) / 2
    for row in good.itertuples():
        dy = 26 if row.Score < mid else -30
        ax.annotate(f"rank {int(row.Rank)} of {len(df)}",
                    xy=(row.L_RMSD, row.Score),
                    xytext=(46, dy), textcoords="offset points",
                    fontsize=10.5, color="#7a3d00", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="#fdf0e2",
                              edgecolor="#d6914a", linewidth=0.8),
                    arrowprops=dict(arrowstyle="->", color="#7a3d00", lw=1.4))

    ax.set_title(f"{p}  —  {name}   ({tag})\n"
                 f"{len(df)} poses    r(score, L-RMSD) = {r:+.3f}",
                 fontsize=12.5)
    ax.set_xlabel("L-RMSD to native complex (Å)      lower = better", fontsize=11)
    ax.set_ylabel("LightDock score      higher = better", fontsize=11)
    ax.tick_params(labelsize=10)
    ax.legend(fontsize=9.5, loc="lower right", framealpha=0.95)

fig.suptitle("The scoring function cannot find its own best poses",
             fontsize=15)
fig.subplots_adjust(top=0.81, bottom=0.13, left=0.07, right=0.98, wspace=0.22)
out = base / "slide_score_vs_rmsd.png"
fig.savefig(out, dpi=200)
print("written:", out)
