configfile: "config/config.yaml"

PAIRS = list(config["pairs"].keys())


SWARMS = range(config["docking"]["swarms"])
STEPS = config["docking"]["steps"]

GLOWWORMS = config["docking"]["glowworms"]

TOP = config["docking"]["top_models"]

rule all:
    input:
        expand("results/{pair}/docking/top_{n}.pdb", pair=PAIRS, n=range(1, TOP + 1)),
        expand("results/{pair}/analysis/scores.csv", pair=PAIRS),
        expand("results/{pair}/analysis/top_scores.png", pair=PAIRS),
        expand("results/{pair}/analysis/top_1.png", pair=PAIRS)


rule preprocess:
    input:
        lambda wc: config["pairs"][wc.pair][wc.part]
    output:
        "results/{pair}/prep/{part}_clean.pdb"
    log:
        "logs/{pair}/preprocess_{part}.log"
    benchmark:
        "benchmarks/{pair}/preprocess_{part}.tsv"
    conda:
        "envs/prep.yaml"
    shell:
        "pdbfixer {input} --output={output} --keep-heterogens=none "
        "--add-residues --add-atoms=all --ph=7.0 > {log} 2>&1"
    

rule lightdock_setup:
    input:
        rec="results/{pair}/prep/receptor_clean.pdb",
        lig="results/{pair}/prep/ligand_clean.pdb"
    output:
        "results/{pair}/docking/setup.json"
    params:
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
        rm -rf {params.outdir} && mkdir -p {params.outdir}
        cp {input.rec} {params.outdir}/receptor.pdb
        cp {input.lig} {params.outdir}/ligand.pdb
        (cd {params.outdir} && lightdock3_setup.py receptor.pdb ligand.pdb \
            -s {params.swarms} -g {params.glowworms} --noxt --noh --now) > {log} 2>&1
        """

rule lightdock_run:
    input:
        "results/{pair}/docking/setup.json"
    output:
        expand("results/{{pair}}/docking/swarm_{s}/gso_{steps}.out", s=SWARMS, steps=STEPS)
    params:
        steps=STEPS,
        outdir="results/{pair}/docking"
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
    input:
        gso="results/{pair}/docking/swarm_{s}/gso_" + str(STEPS) + ".out"
    output:
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

rule parse_results:
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
    input:
        model="results/{pair}/docking/top_1.pdb"
    output:
        "results/{pair}/analysis/top_1.png"
    params:
        receptor="results/{pair}/docking/lightdock_receptor.pdb"
    log:
        "logs/{pair}/visualize.log"
    benchmark:
        "benchmarks/{pair}/visualize.tsv"
    conda:
        "envs/pymol.yaml"
    script:
        "scripts/visualize.py"