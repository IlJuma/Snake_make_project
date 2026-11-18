"""Compare docked poses with the experimental (native) complex.

For every pose of a pair this computes the standard CAPRI docking metrics:
  - L-RMSD : superpose the receptor of the pose on the native receptor and
             measure how far the ligand backbone is from the native ligand
             ("is the ligand in the right place?")
  - I-RMSD : backbone RMSD of the native interface residues (within 10 A of
             the partner) after superposing only those residues
  - fnat   : fraction of native residue-residue contacts (heavy atoms within
             5 A) that the pose reproduces
  - DockQ  : one number between 0 (wrong) and 1 (perfect) combining the three
  - CAPRI  : quality class: incorrect / acceptable / medium / high

The docking inputs are unbound structures, so residue numbers and lengths do
not match the native complex. Residues are therefore paired by a sequence
alignment of each chain, not by their number.
"""

import sys

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

BACKBONE = ("N", "CA", "C", "O")
CONTACT_CUTOFF = 5.0      # A, heavy atoms, for fnat
INTERFACE_CUTOFF = 10.0   # A, heavy atoms, for I-RMSD

_parser = PDBParser(QUIET=True)


# ---------------------------------------------------------------- structure
def load_chain(pdb_file, chain_id):
    """Return the amino-acid residues of one chain (no waters/ligands)."""
    chain = _parser.get_structure("s", pdb_file)[0][chain_id]
    return [res for res in chain if res.id[0] == " " and "CA" in res]


def heavy_atoms(residue):
    return [a for a in residue if a.element not in ("H", "D")]


def map_residues(model_res, native_res):
    """Pair model and native residues with a global sequence alignment.

    Returns a list of (model_residue, native_residue). Gaps at the ends are
    free, because the unbound structures often miss or add terminal residues.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    aligner.end_gap_score = 0

    seq_m = "".join(seq1(r.get_resname()) for r in model_res)
    seq_n = "".join(seq1(r.get_resname()) for r in native_res)
    aln = aligner.align(seq_m, seq_n)[0]

    pairs = []
    for (m0, m1), (n0, n1) in zip(*aln.aligned):
        pairs.extend(zip(model_res[m0:m1], native_res[n0:n1]))
    return pairs


def backbone_coords(pairs):
    """Backbone coordinates of paired residues (only atoms present in both)."""
    model, native = [], []
    for rm, rn in pairs:
        for name in BACKBONE:
            if name in rm and name in rn:
                model.append(rm[name].coord)
                native.append(rn[name].coord)
    return np.array(model, dtype=float), np.array(native, dtype=float)


# ------------------------------------------------------------------- maths
def kabsch(P, Q):
    """Rotation R and translation t that best map P onto Q (P @ R.T + t ~ Q)."""
    p0, q0 = P.mean(axis=0), Q.mean(axis=0)
    H = (P - p0).T @ (Q - q0)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))          # avoid a mirror image
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, q0 - p0 @ R.T


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(axis=1).mean()))


def superposed_rmsd(P, Q):
    R, t = kabsch(P, Q)
    return rmsd(P @ R.T + t, Q)


def contacts(rec_residues, lig_residues, cutoff):
    """Set of (i, j) index pairs whose heavy atoms are closer than cutoff."""
    def stack(residues):
        coords, owner = [], []
        for i, res in enumerate(residues):
            for atom in heavy_atoms(res):
                coords.append(atom.coord)
                owner.append(i)
        return np.array(coords, dtype=float), np.array(owner)

    rc, ri = stack(rec_residues)
    lc, li = stack(lig_residues)
    dist = np.linalg.norm(rc[:, None, :] - lc[None, :, :], axis=2)
    close_r, close_l = np.nonzero(dist < cutoff)
    return set(zip(ri[close_r].tolist(), li[close_l].tolist()))


def capri_class(fnat, lrms, irms):
    if fnat >= 0.5 and (lrms <= 1.0 or irms <= 1.0):
        return "high"
    if fnat >= 0.3 and (lrms <= 5.0 or irms <= 2.0):
        return "medium"
    if fnat >= 0.1 and (lrms <= 10.0 or irms <= 4.0):
        return "acceptable"
    return "incorrect"


def dockq(fnat, lrms, irms):
    return (fnat + 1 / (1 + (irms / 1.5) ** 2) + 1 / (1 + (lrms / 8.5) ** 2)) / 3


# ---------------------------------------------------------------- one pose
def evaluate_pose(model_pdb, native_pdb, rec_chain, lig_chain,
                  model_rec_chain="R", model_lig_chain="L"):
    """All metrics for one pose (chains R/L) against the native complex."""
    rec_pairs = map_residues(load_chain(model_pdb, model_rec_chain), load_chain(native_pdb, rec_chain))
    lig_pairs = map_residues(load_chain(model_pdb, model_lig_chain), load_chain(native_pdb, lig_chain))

    # L-RMSD: fit on the receptor, measure the ligand without refitting
    rec_m, rec_n = backbone_coords(rec_pairs)
    lig_m, lig_n = backbone_coords(lig_pairs)
    R, t = kabsch(rec_m, rec_n)
    lrms = rmsd(lig_m @ R.T + t, lig_n)

    # fnat: native contacts (on paired residues) that the pose also has
    native_contacts = contacts([n for _, n in rec_pairs], [n for _, n in lig_pairs], CONTACT_CUTOFF)
    model_contacts = contacts([m for m, _ in rec_pairs], [m for m, _ in lig_pairs], CONTACT_CUTOFF)
    fnat = len(native_contacts & model_contacts) / len(native_contacts)

    # I-RMSD: backbone of the native interface residues, both partners together
    interface = contacts([n for _, n in rec_pairs], [n for _, n in lig_pairs], INTERFACE_CUTOFF)
    rec_if = sorted({i for i, _ in interface})
    lig_if = sorted({j for _, j in interface})
    if_m1, if_n1 = backbone_coords([rec_pairs[i] for i in rec_if])
    if_m2, if_n2 = backbone_coords([lig_pairs[j] for j in lig_if])
    irms = superposed_rmsd(np.vstack([if_m1, if_m2]), np.vstack([if_n1, if_n2]))

    return {
        "L_RMSD": round(lrms, 3),
        "I_RMSD": round(irms, 3),
        "fnat": round(fnat, 3),
        "DockQ": round(dockq(fnat, lrms, irms), 3),
        "CAPRI": capri_class(fnat, lrms, irms),
        "Mapped_rec": len(rec_pairs),
        "Mapped_lig": len(lig_pairs),
        "Native_contacts": len(native_contacts),
    }


# --------------------------------------------------------------- snakemake
def main(smk):
    import matplotlib
    matplotlib.use("Agg")          # draw plots without a screen
    import matplotlib.pyplot as plt
    import pandas as pd

    sys.stderr = sys.stdout = open(smk.log[0], "w")
    pair = smk.wildcards.pair
    native = smk.params.native

    scores = pd.read_csv(smk.input.scores).set_index("Rank")

    rows = []
    # models are given in rank order: top_1.pdb, top_2.pdb, ...
    for rank, model in enumerate(smk.input.models, start=1):
        m = evaluate_pose(model, native["pdb"], native["receptor_chain"], native["ligand_chain"])
        rows.append({"Rank": rank, "Model": model, "Score": scores.loc[rank, "Scoring"], **m})
        print(f"{pair} top_{rank}: L-RMSD {m['L_RMSD']:6.2f} A  I-RMSD {m['I_RMSD']:6.2f} A  "
              f"fnat {m['fnat']:.2f}  DockQ {m['DockQ']:.2f}  {m['CAPRI']}")

    df = pd.DataFrame(rows)
    df.to_csv(smk.output.csv, index=False)
    print(f"{pair}: {df['Mapped_rec'].iloc[0]} receptor and {df['Mapped_lig'].iloc[0]} ligand "
          f"residues paired with the native, {df['Native_contacts'].iloc[0]} native contacts")

    # --- Plot: score vs L-RMSD, coloured by CAPRI class ---
    colors = {"high": "tab:green", "medium": "tab:blue", "acceptable": "tab:orange", "incorrect": "tab:gray"}
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for cls, color in colors.items():
        sub = df[df["CAPRI"] == cls]
        ax.scatter(sub["L_RMSD"], sub["Score"], color=color, label=cls, s=60, zorder=3)
    for _, r in df.iterrows():
        ax.annotate(f"top_{r['Rank']}", (r["L_RMSD"], r["Score"]),
                    textcoords="offset points", xytext=(5, 4), fontsize=8)
    for x, name in [(1, "high"), (5, "medium"), (10, "acceptable")]:
        ax.axvline(x, color="lightgray", linestyle="--", zorder=1)
        ax.text(x, ax.get_ylim()[1], f" {name}", fontsize=8, color="gray", va="top")
    ax.set_xlabel("L-RMSD to native complex (A, lower = better)")
    ax.set_ylabel("LightDock score (higher = better)")
    ax.set_title(f"{pair}: docked poses vs native structure")
    ax.legend(title="CAPRI class", fontsize=8)
    fig.tight_layout()
    fig.savefig(smk.output.plot, dpi=150)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821  (injected by Snakemake)
