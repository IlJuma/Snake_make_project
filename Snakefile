#Reading the config file with the configs needed for snakemake like the pairs of ligand and
#receptor, the parameters for the docking. all is stored in a nested dict.
#The config file makes it easy to configure the run at a central point.
configfile: "config/config.yaml"

#makes a list of all protein pairs which are tested for binding
PAIRS = list(config["pairs"].keys())

#Produces the swarms needed for docking number can be defined in config file
SWARMS = range(config["docking"]["swarms"])

#Writes numbers from the config file into a variable for easier acces during later workflow
STEPS = config["docking"]["steps"]
GLOWWORMS = config["docking"]["glowworms"]
TOP = config["docking"]["top_models"]

#Is the target rule -> with only input, at the end of the workflow these files need to exist
#with expand the wildcards are inserted and creating a file/path for each value in the wildcard
rule all:
    input:
        expand("results/{pair}/docking/top_{n}.pdb", pair=PAIRS, n=range(1, TOP + 1)),
        #same poses but with receptor and ligand as two chains (ready for viewers/RMSD)
        expand("results/{pair}/complex/top_{n}.pdb", pair=PAIRS, n=range(1, TOP + 1)),
        expand("results/{pair}/analysis/scores.csv", pair=PAIRS),
        expand("results/{pair}/analysis/top_scores.png", pair=PAIRS),
        expand("results/{pair}/analysis/top_1.png", pair=PAIRS)

#Rule which is responsible for the preprocessing the protein structures: removes water and non
#non-protein molecules, add missing residues, add missing heavy atoms and add hydrogens to match pH7
rule preprocess:
    input:
        #creates name for the file which need to be used (pdb file) needs a function because
        #file names don't follow a pattern which can be used with wildcards wc.pair defines the 
        #protein and wc.part if the receptor or the ligand is called
        lambda wc: config["pairs"][wc.pair][wc.part]
    output:
        #name of the outputfile which is generated with the wildcards for all proteins and for 
        #part defined in input from rule lightdock_setup
        "results/{pair}/prep/{part}_clean.pdb"
    log:
        #creates log file for each time rule is runned (for each wildcard)
        #everything the command is printing goes into log file, needed for error handling
        "logs/{pair}/preprocess_{part}.log"
    benchmark:
        #creates benchmark file for each time rule is runned
        #documents runtime and resources needed for optimization and identifing slow steps
        "benchmarks/{pair}/preprocess_{part}.tsv"
    conda:
        #acitvates the prepared enviroment needed for this job
        "envs/prep.yaml"
    shell:
        #Snakemake replaceec the placeholders and runs the command written in ""
        "pdbfixer {input} --output={output} --keep-heterogens=none "
        "--add-residues --add-atoms=all --ph=7.0 > {log} 2>&1"
    

rule lightdock_setup:
    #This rule prepares the search by checking the surface from receptor and placing swarms
    #with glowworms: copies of the ligand, each in a random position and orientation.
    input:
        #needs both files (receptor and ligand for each pair)
        rec="results/{pair}/prep/receptor_clean.pdb",
        lig="results/{pair}/prep/ligand_clean.pdb"
    output:
        #Lightdock produced multiple files (setup.json, swarm/0....) only one needs to be declared 
        #for snakemake as signe that the step is finished
        "results/{pair}/docking/setup.json"
    params:
        #extracts the parameter for docking from the config file
        swarms=config["docking"]["swarms"],
        glowworms=config["docking"]["glowworms"],
        outdir="results/{pair}/docking"
    log:
        "logs/{pair}/lightdock_setup.log"
    benchmark:
        "benchmarks/{pair}/lightdock_setup.tsv"
    conda:
        "envs/lightdock.yaml"
    shell:
        """
        #starts from an empty folder (removes old stuff from prev runs)
        rm -rf {params.outdir} && mkdir -p {params.outdir}
        #copy the clean structures into the folder
        cp {input.rec} {params.outdir}/receptor.pdb
        cp {input.lig} {params.outdir}/ligand.pdb
        #runs the setup in the folder  (Lightdock writes into the current folder)
        (cd {params.outdir} && lightdock3_setup.py receptor.pdb ligand.pdb \
            -s {params.swarms} -g {params.glowworms} --noxt --noh --now) > {log} 2>&1
        """

rule lightdock_run:
    #This rule runs the actual docking (calculates a score for each glowworm how well the ligand 
    #fits, moves each glowworm a bit toward a brighter neighbor in the swarm)
    input:
        #one run per set-up defined before
        "results/{pair}/docking/setup.json"
    output:
        #one output file per swarm the {{pair}} keeps this wc for multiple wc {s/steps}
        expand("results/{{pair}}/docking/swarm_{s}/gso_{steps}.out", s=SWARMS, steps=STEPS)
    params:
        steps=STEPS,
        outdir="results/{pair}/docking"
    #this says it can use up to 4 cores
    threads: 4
    log:
        "logs/{pair}/lightdock_run.log"
    benchmark:
        "benchmarks/{pair}/lightdock_run.tsv"
    conda:
        "envs/lightdock.yaml"
    shell:
        "(cd {params.outdir} && lightdock3.py setup.json {params.steps} -c {threads}) > {log} 2>&1"

rule generate_models:
    #This rule builds a real 3D structure for each glowworm: it takes the receptor, 
    #moves and rotates the ligand according to the numbers, and writes both into a PDB file.
    input:
        #creates the file name by combining strings with +
        gso="results/{pair}/docking/swarm_{s}/gso_" + str(STEPS) + ".out"
    output:
        #one file per glowworm
        expand("results/{{pair}}/docking/swarm_{{s}}/lightdock_{i}.pdb", i=range(GLOWWORMS))
    params:
        swarmdir="results/{pair}/docking/swarm_{s}",
        steps=STEPS,
        glowworms=GLOWWORMS
    log:
        "logs/{pair}/generate_models_swarm_{s}.log"
    benchmark:
        "benchmarks/{pair}/generate_models_swarm_{s}.tsv"
    conda:
        "envs/lightdock.yaml"
    shell:
        "(cd {params.swarmdir} && lgd_generate_conformations.py ../receptor.pdb ../ligand.pdb "
        "gso_{params.steps}.out {params.glowworms}) > {log} 2>&1"

rule cluster:
    #Glowworms in a swarm move toward each other during the docking, so many of them end up at 
    #almost the same pose. Ranking all 50 would fill the top of the list with near-duplicates. 
    #Clustering groups similar poses and keeps only the best one per group (the representative).
    input:
        gso="results/{pair}/docking/swarm_{s}/gso_" + str(STEPS) + ".out",
        models="results/{pair}/docking/swarm_{s}/lightdock_0.pdb"
    output:
        "results/{pair}/docking/swarm_{s}/cluster.repr"
    params:
        swarmdir="results/{pair}/docking/swarm_{s}",
        steps=STEPS
    log:
        "logs/{pair}/cluster_swarm_{s}.log"
    benchmark:
        "benchmarks/{pair}/cluster_swarm_{s}.tsv"
    conda:
        "envs/lightdock.yaml"
    shell:
        "(cd {params.swarmdir} && lgd_cluster_bsas.py gso_{params.steps}.out) > {log} 2>&1"


rule rank:
#lgd_rank.py (comes from LightDock) goes through all swarm folders, collects the cluster representatives from each 
#cluster.repr, looks up their full data in the gso_10.out files, and sorts them all together 
#by score, best first.
    input:
        expand("results/{{pair}}/docking/swarm_{s}/cluster.repr", s=SWARMS)
    output:
        "results/{pair}/docking/rank_by_scoring.list"
    params:
        outdir="results/{pair}/docking",
        swarms=len(SWARMS),
        steps=STEPS
    log:
        "logs/{pair}/rank.log"
    benchmark:
        "benchmarks/{pair}/rank.tsv"
    conda:
        "envs/lightdock.yaml"
    shell:
        "(cd {params.outdir} && lgd_rank.py {params.swarms} {params.steps}) > {log} 2>&1"

rule top_models:
#The ranking tells you which poses are best, but the matching model files 
#are scattered over the swarm folders (swarm_0/lightdock_17.pdb, swarm_3/lightdock_39.pdb, …). 
#This rule creates the N best poses as PDB files with simple names in the docking folder: 
#top_1.pdb is the best pose, top_2.pdb the second best, and so on.
    input:
        "results/{pair}/docking/rank_by_scoring.list"
    output:
        expand("results/{{pair}}/docking/top_{n}.pdb", n=range(1, TOP + 1))
    params:
        outdir="results/{pair}/docking",
        top=TOP
    log:
        "logs/{pair}/top_models.log"
    benchmark:
        "benchmarks/{pair}/top_models.tsv"
    conda:
        "envs/lightdock.yaml"
    shell:
        "(cd {params.outdir} && lgd_top.py receptor.pdb ligand.pdb "
        "rank_by_scoring.list {params.top}) > {log} 2>&1"

rule split_chains:
# lgd_top.py writes receptor and ligand as one single chain A, so viewers like
# ChimeraX or PyMOL load the pose as one molecule. This rule writes the same pose
# with the receptor as chain R and the ligand as chain L (separated by TER), which
# is what downstream steps (visualisation, RMSD, interface analysis) need.
    input:
        model="results/{pair}/docking/top_{n}.pdb",
        #the receptor structure LightDock really used, defines where the ligand starts
        receptor="results/{pair}/docking/lightdock_receptor.pdb"
    output:
        "results/{pair}/complex/top_{n}.pdb"
    log:
        "logs/{pair}/split_chains_top_{n}.log"
    benchmark:
        "benchmarks/{pair}/split_chains_top_{n}.tsv"
    conda:
        "envs/pymol.yaml"
    script:
        "scripts/split_chains.py"

rule parse_results:
# Parse the LightDock ranking (rank_by_scoring.list) of each pair:
# - save all ranked poses with rank, swarm, glowworm, model file and score as a CSV table
# - plot the 20 best scores as a bar chart (top N models highlighted in green)
    input:
        rank="results/{pair}/docking/rank_by_scoring.list"
    output:
        csv="results/{pair}/analysis/scores.csv",
        plot="results/{pair}/analysis/top_scores.png"
    params:
        top=TOP
    log:
        "logs/{pair}/parse_results.log"
    benchmark:
        "benchmarks/{pair}/parse_results.tsv"
    conda:
        "envs/pymol.yaml"
    script:
        "scripts/parse_results.py"

rule visualize:
# Render an image of the best docked pose (top_1.pdb) of each pair with PyMOL:
# - show both proteins as cartoon (receptor blue, ligand orange)
# - show the interface residues (within 4 Å of the other protein) as sticks
    input:
        #two-chain version written by rule split_chains
        model="results/{pair}/complex/top_1.pdb"
    output:
        "results/{pair}/analysis/top_1.png"
    log:
        "logs/{pair}/visualize.log"
    benchmark:
        "benchmarks/{pair}/visualize.tsv"
    conda:
        "envs/pymol.yaml"
    script:
        "scripts/visualize.py"