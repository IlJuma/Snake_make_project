"""Master summary: ceiling vs achieved, for all four pairs.

Combines
  - positive_control_results.txt : the best result achievable with these inputs
  - *_all_poses_rmsd.csv         : what the docking actually produced
"""
import re
from pathlib import Path
import pandas as pd

base = Path(__file__).parent
s100 = base / "s100"
ctrl_file = base / "positive_control_results.txt"

NAMES = {"1efn": "HIV-1 Nef / Fyn SH3",
         "2g45": "USP5 ZnF / ubiquitin",
         "2goo": "BMP-2 / ActRII",
         "3fap": "FKBP12 / FRB"}

# --- parse the positive control ---
ceiling = {}
for line in ctrl_file.read_text().splitlines():
    m = re.match(r"\s*(\w+)\s+([\d.]+) A\s+([\d.]+) A\s+([\d.]+)\s+([\d.]+)\s+"
                 r"([\d.]+)\s+([\d.]+)\s+(\w+)", line)
    if m:
        p, br, bl, lr, ir, fn, dq, cls = m.groups()
        ceiling[p] = dict(d_rec=float(br), d_lig=float(bl), I_RMSD=float(ir),
                          fnat=float(fn), DockQ=float(dq), CAPRI=cls)

rows = []
for p, name in NAMES.items():
    df = pd.read_csv(s100 / f"{p}_all_poses_rmsd.csv")
    top5 = df.nsmallest(5, "Rank")
    good = df[df["CAPRI"] != "incorrect"]
    best = df.loc[df["DockQ"].idxmax()]
    c = ceiling[p]

    rows.append({
        "pair": p,
        "complex": name,
        "conf. change rec (A)": c["d_rec"],
        "conf. change lig (A)": c["d_lig"],
        "CEILING class": c["CAPRI"],
        "CEILING DockQ": round(c["DockQ"], 3),
        "poses": len(df),
        "acceptable found": len(good),
        "best DockQ found": round(best["DockQ"], 3),
        "best L-RMSD (A)": round(df["L_RMSD"].min(), 2),
        "rank of best pose": int(best["Rank"]),
        "in top 5?": "yes" if best["Rank"] <= 5 else "no",
        "top5 best L-RMSD (A)": round(top5["L_RMSD"].min(), 2),
        "top5 score spread": round(top5["Score"].max() - top5["Score"].min(), 2),
        "r(score, L-RMSD)": round(df["Score"].corr(df["L_RMSD"]), 3),
    })

out = pd.DataFrame(rows)
out.to_csv(base / "master_summary.csv", index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
print(out.to_string(index=False))

print("\n" + "=" * 78)
print("READING: 'CEILING' is the best result achievable with these unbound inputs,")
print("measured by superposing each partner onto its native chain.")
print()
for r in rows:
    if r["CEILING class"] == "incorrect":
        print(f"{r['pair']}: ceiling is 'incorrect' -- the receptor moves "
              f"{r['conf. change rec (A)']} A on binding, so rigid-body docking")
        print(f"      cannot reproduce the interface (fnat = 0 even for a perfect pose).")
    else:
        print(f"{r['pair']}: ceiling '{r['CEILING class']}' (DockQ {r['CEILING DockQ']}) "
              f"-> winnable. Found {r['acceptable found']} acceptable pose(s); "
              f"best ranked {r['rank of best pose']} of {r['poses']}.")
