"""POSITIVE CONTROL for the RMSD reference structures.

Question: are the native_*.pdb files a valid reference for our docked models?

Test: build the *ideal* answer. Take the same preprocessed unbound receptor and
ligand that LightDock docks, and superpose each one separately onto its chain in
the native complex. The result is the best pose the docking could ever produce
with these inputs. Feed it through the exact same evaluate_pose() the workflow
uses.

  - If it comes back "high" with L-RMSD near 0, the reference is sound and the
    bad results are genuine docking failure.
  - If it comes back poor, the reference (or the chain assignment) is wrong and
    every RMSD number in the project is suspect.
"""
import sys
sys.path.insert(0, "scripts")

import numpy as np
import rmsd as R

PAIRS = {
    "1efn": ("B", "A"),
    "2g45": ("A", "B"),
    "2goo": ("A", "C"),
    "3fap": ("A", "B"),
}


def superpose_onto(mobile_pdb, native_pdb, native_chain):
    """Superpose a single-chain PDB onto a chain of the native, via sequence mapping.

    Returns (all_atom_lines, rotation, translation, backbone_rmsd, n_mapped).
    """
    mob_res = R.load_chain(mobile_pdb, "A")
    nat_res = R.load_chain(native_pdb, native_chain)
    pairs = R.map_residues(mob_res, nat_res)
    P, Q = R.backbone_coords(pairs)          # P = mobile, Q = native
    rot, tran = R.kabsch(P, Q)
    bb = R.rmsd(P @ rot.T + tran, Q)
    return rot, tran, bb, len(pairs)


def transform_lines(pdb_file, rot, tran, chain):
    out = []
    for line in open(pdb_file):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        x, y, z = xyz @ rot.T + tran
        out.append(f"{line[:21]}{chain}{line[22:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
    return out


print(f"{'pair':6} {'bb-fit rec':>11} {'bb-fit lig':>11} "
      f"{'L-RMSD':>8} {'I-RMSD':>8} {'fnat':>6} {'DockQ':>7}  CAPRI")
print("-" * 76)

for pair, (rc, lc) in PAIRS.items():
    native = f"data/raw/{pair}_native.pdb"
    rec = f"results/{pair}/prep/receptor_clean.pdb"
    lig = f"results/{pair}/prep/ligand_clean.pdb"

    rot_r, tran_r, bb_r, nr = superpose_onto(rec, native, rc)
    rot_l, tran_l, bb_l, nl = superpose_onto(lig, native, lc)

    ideal = f"/tmp/{pair}_ideal.pdb"
    with open(ideal, "w") as fh:
        fh.writelines(transform_lines(rec, rot_r, tran_r, "R"))
        fh.write("TER\n")
        fh.writelines(transform_lines(lig, rot_l, tran_l, "L"))
        fh.write("TER\nEND\n")

    m = R.evaluate_pose(ideal, native, rc, lc)
    print(f"{pair:6} {bb_r:8.2f} A  {bb_l:8.2f} A  "
          f"{m['L_RMSD']:7.2f}  {m['I_RMSD']:7.2f}  {m['fnat']:6.2f}  "
          f"{m['DockQ']:7.3f}  {m['CAPRI']}")

print("\nbb-fit = backbone RMSD of the unbound chain fitted onto its native chain")
print("         (this is the unbound->bound conformational change, our hard floor)")
