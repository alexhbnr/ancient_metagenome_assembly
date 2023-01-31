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
    conda: "../envs/ENVS_minimap2.yaml"
    resources:
        mem = 4,
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
    conda: "../envs/ENVS_quast.yaml"
    resources:
        mem = lambda wildcards, attempt: 24 + attempt * 24,
        cores = 8,
        metaquast = 1
    params:
        outdir = "{resultdir}/stats//metaquast/{sample}-{assembler}"
    threads: 8
    shell:
        """
        metaquast.py \
            -o {params.outdir} \
            --threads {threads} \
            --ambiguity-usage all \
            --no-icarus \
            --no-read-stats \
            {input} && touch {output}
        """
################################################################################
