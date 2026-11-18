"""Sanity checks for scripts/rmsd.py.

Run from the project root with the rmsd environment:
    python -m pytest tests/
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import rmsd  # noqa: E402

NATIVE = ROOT / "data/raw/1efn_native.pdb"
REC, LIG = "B", "A"  # Nef = receptor, SH3 = ligand


def write_model(path, rec_shift=0, lig_shift=0, lig_move=(0.0, 0.0, 0.0), drop_first=0):
    """Write chains B/A of the native as a fake docking pose (chains R/L).

    rec_shift / lig_shift renumber the residues, lig_move translates the
    ligand, drop_first removes the first residues of each chain: the script
    must cope with all of these, like it has to with real poses.
    """
    out = []
    seen = {"R": set(), "L": set()}
    for line in open(NATIVE):
        if not line.startswith("ATOM") or line[21] not in (REC, LIG):
            continue
        chain = "R" if line[21] == REC else "L"
        resnum = int(line[22:26])
        seen[chain].add(resnum)
        if len(seen[chain]) <= drop_first:
            continue
        shift = rec_shift if chain == "R" else lig_shift
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        if chain == "L":
            x, y, z = x + lig_move[0], y + lig_move[1], z + lig_move[2]
        line = (line[:21] + chain + f"{resnum + shift:4d}" + line[26:30]
                + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
        out.append(line)
    path.write_text("".join(out) + "END\n")
    return path


def evaluate(model):
    return rmsd.evaluate_pose(str(model), str(NATIVE), REC, LIG)


def test_native_against_itself_is_perfect(tmp_path):
    m = evaluate(write_model(tmp_path / "self.pdb"))
    assert m["L_RMSD"] == pytest.approx(0, abs=1e-3)
    assert m["I_RMSD"] == pytest.approx(0, abs=1e-3)
    assert m["fnat"] == pytest.approx(1)
    assert m["DockQ"] == pytest.approx(1)
    assert m["CAPRI"] == "high"


def test_renumbered_and_truncated_model_still_matches(tmp_path):
    # different numbering and missing residues, like the unbound inputs
    m = evaluate(write_model(tmp_path / "renum.pdb", rec_shift=100, lig_shift=-50, drop_first=3))
    assert m["L_RMSD"] == pytest.approx(0, abs=1e-3)
    assert m["fnat"] > 0.9
    assert m["CAPRI"] == "high"


def test_ligand_moved_10A_gives_lrmsd_10(tmp_path):
    m = evaluate(write_model(tmp_path / "moved.pdb", lig_move=(10.0, 0.0, 0.0)))
    assert m["L_RMSD"] == pytest.approx(10, abs=1e-2)
    assert m["fnat"] < 0.5
    assert m["CAPRI"] in ("incorrect", "acceptable")


def test_kabsch_recovers_rotation():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(20, 3))
    theta = 0.7
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    Q = P @ R.T + np.array([1.0, -2.0, 3.0])
    assert rmsd.superposed_rmsd(P, Q) == pytest.approx(0, abs=1e-8)


@pytest.mark.parametrize("fnat,lrms,irms,expected", [
    (0.6, 0.8, 0.9, "high"),
    (0.4, 4.0, 3.0, "medium"),
    (0.2, 9.0, 5.0, "acceptable"),
    (0.05, 2.0, 1.0, "incorrect"),   # fnat too low, whatever the RMSD
    (0.5, 15.0, 6.0, "incorrect"),
])
def test_capri_classes(fnat, lrms, irms, expected):
    assert rmsd.capri_class(fnat, lrms, irms) == expected
