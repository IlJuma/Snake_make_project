"""
Render an image of the best docked pose of one protein pair with PyMOL.

Called by the Snakemake rule `visualize` (one run per pair).

Input:
    top_1.pdb              - best-ranked docking model (receptor + docked ligand)
    lightdock_receptor.pdb - receptor as prepared by LightDock (used to count its atoms)
Output:
    top_1.png              - image: receptor blue, ligand orange, interface residues as sticks
"""

import sys

from pymol import cmd

#send print() statements and error messages to the lop file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

#Get the paths from the snakemake file to the needed files
model_file = snakemake.input.model
png_file = snakemake.output[0]
pair = snakemake.wildcards.pair

# --- Label receptor and ligand as separate chains ---
# LightDock writes both proteins as chain A, with overlapping residue numbers.
# Problem is PyMol would merge residues from both proteins when they are overlapping 
# --> which is not correct
# The receptor comes first in the file, so it count its atoms
# and rename the chains: R = receptor, L = ligand. 
# to be able to seperate both for creating the picture


#counts the atoms in the pdb file from the receptor from the receptor
n_receptor = sum(1 for line in open(receptor_file) if line.startswith(("ATOM", "HETATM")))

#Renames the chain from the first n-atoms which belong to the receptor and the rest is Ligand
lines = []
atom_count = 0
for line in open(model_file):
    if line.startswith(("ATOM", "HETATM")):
        atom_count += 1
        chain = "R" if atom_count <= n_receptor else "L"
        line = line[:21] + chain + line[22:]          # column 22 = chain ID
    lines.append(line)

cmd.read_pdbstr("".join(lines), "complex")
# The pose already has the receptor as chain R and the ligand as chain L
# (written by rule split_chains).
cmd.load(model_file, "complex")

# --- Style ---
cmd.dss("complex")                                  # calculate helices and sheets
cmd.hide("everything")
cmd.show("cartoon")                 #style of the protein to visualize
cmd.color("skyblue", "chain R")     #color of the receptor
cmd.color("orange", "chain L")      #color of the ligand

# Residues within 4 Å of the other protein = the interface are shown as sticks
cmd.select("interface", "byres ((chain L within 4 of chain R) or (chain R within 4 of chain L))")
cmd.show("sticks", "interface and not hydro and not name N+C+O")
print(f"{pair}: {cmd.count_atoms('interface and name CA')} interface residues")

# --- Render ---
cmd.orient("complex")
cmd.bg_color("white")
cmd.png(png_file, width=1200, height=900, dpi=150, ray=1)