"""Rewrite a LightDock pose so receptor and ligand are two separate chains.

lgd_top.py writes the whole complex as a single chain A with overlapping
residue numbers, so viewers (ChimeraX, PyMOL) see one molecule. The receptor
atoms always come first in the file, so counting the atoms of the receptor
structure LightDock actually used (lightdock_receptor.pdb, heavy atoms only)
gives the split point. Receptor becomes chain R, ligand chain L, separated by
a TER record.
"""

import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

model_file = snakemake.input.model
receptor_file = snakemake.input.receptor
out_file = snakemake.output[0]

ATOM_RECORDS = ("ATOM", "HETATM")

# number of atoms in the receptor = position of the receptor/ligand border
n_receptor = sum(1 for line in open(receptor_file) if line.startswith(ATOM_RECORDS))

atom_count = 0
wrote_ter = False

with open(model_file) as fh, open(out_file, "w") as out:
    for line in fh:
        if line.startswith(ATOM_RECORDS):
            atom_count += 1
            # column 22 (index 21) holds the chain ID
            if atom_count <= n_receptor:
                line = line[:21] + "R" + line[22:]
            else:
                # close the receptor chain before the first ligand atom
                if not wrote_ter:
                    out.write("TER\n")
                    wrote_ter = True
                line = line[:21] + "L" + line[22:]
        elif line.startswith(("TER", "END")):
            # drop the original records, we write our own
            continue
        out.write(line)
    out.write("TER\nEND\n")

n_ligand = atom_count - n_receptor
print(f"{snakemake.wildcards.pair} {out_file}: chain R = {n_receptor} atoms, chain L = {n_ligand} atoms")

if n_ligand <= 0:
    raise ValueError(
        f"no ligand atoms left after the split: model has {atom_count} atoms, "
        f"receptor has {n_receptor}"
    )
