rule quality_summary:
    input:
        caln50  = lambda wildcards: expand("{resultdir}/stats/caln50/{sample}-{assembler}.caln", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards)),
        quast = lambda wildcards: expand("{resultdir}/stats/metaquast/{sample}-{assembler}/report.html", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards))
    output:
        touch(f"{config['tmpdir']}/quality_summary.done")

#### CalN50 ####################################################################

localrules: download_caln50_script

rule download_caln50_script:
    output:
        "{tmpdir}/scripts/caln50.js"
    message: "Download the calN50 script from GitHub"
    params:
        url = "https://raw.githubusercontent.com/lh3/calN50/master/calN50.js"
    shell:
        "wget -O {output} {params.url}"

rule caln50:
    input:
        prg = lambda wildcards: f"{config['tmpdir']}/scripts/caln50.js",
        fasta = "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
    output:
        "{resultdir}/stats/caln50/{sample}-{assembler}.caln"
    message: "Calculate N50 related statistics of the contigs: {wildcards.sample}"
    container: "docker://quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
    resources:
        mem_mb = 4000,
        runtime = 120,
        slurm_partition = "short",
        cores = 1
    shell:
        """
        k8 {input.prg} {input.fasta} > {output}
        """

################################################################################

#### MetaQUAST #################################################################

rule metaQUAST:
    input:
        "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
    output:
        "{resultdir}/stats/metaquast/{sample}-{assembler}/report.html"
    message: "Run metaQUAST on the contigs: {wildcards.sample}"
    container: "docker://quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"
    resources:
        mem_mb = lambda wildcards, attempt: 24000 + attempt * 24000,
        runtime = 1400,
        slurm_partition = "standard",
        cores = 8,
        metaquast = 1
    params:
        outdir = "{resultdir}/stats//metaquast/{sample}-{assembler}",
        reffa = lambda wildcards: f"-r {config['reference_database']}" if config['reference_database'] != "" else "", 
    threads: 8
    shell:
        """
        metaquast.py \
            -o {params.outdir} \
            {params.reffa} \
            --threads {threads} \
            --ambiguity-usage all \
            --no-icarus \
            --no-read-stats \
            {input} && touch {output}
        """
################################################################################
