import sys

from pymol import cmd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

model_file = snakemake.input.model
png_file = snakemake.output[0]
pair = snakemake.wildcards.pair

# The pose already has the receptor as chain R and the ligand as chain L
# (written by rule split_chains).
cmd.load(model_file, "complex")

# --- Style ---
cmd.dss("complex")                                  # calculate helices and sheets
cmd.hide("everything")
cmd.show("cartoon")
cmd.color("skyblue", "chain R")
cmd.color("orange", "chain L")

# Residues within 4 Å of the other protein = the interface
cmd.select("interface", "byres ((chain L within 4 of chain R) or (chain R within 4 of chain L))")
cmd.show("sticks", "interface and not hydro and not name N+C+O")
print(f"{pair}: {cmd.count_atoms('interface and name CA')} interface residues")

# --- Render ---
cmd.orient("complex")
cmd.bg_color("white")
cmd.png(png_file, width=1200, height=900, dpi=150, ray=1)