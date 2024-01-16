import yaml

import pandas as pd


#### Auxilliary functions ######################################################

def path_to_r(sample, tmpdir, suffix, prefix):
    path = sampletsv.at[sample, suffix]
    if type(path) is float or path == "" or path == "NA":
        return ""
    else:
        if config['readcorrection']:
            return f"-{prefix} {tmpdir}/error_correction/{sample}-wreadcorr_{suffix[-1]}.fastq.gz"
        else:
            return f"-{prefix} {path}"


def path_to_corrected_r(wildcards, suffix):
    with open(checkpoints.error_correction.get(**wildcards).output[0], "rt") as yamlfile:
        corr = yaml.safe_load(yamlfile)
    if suffix == 1:
        return corr[0]['left reads'][0]
    elif suffix == 2:
        return corr[0]['right reads'][0]
    elif suffix == 0:
        return " ".join(sorted(corr[0]['single reads']))

################################################################################

localrules: assembly_workflow

checkpoint assembly_workflow:
    input:
        f"{config['tmpdir']}/sampletsv.validated",
        expand("{tmpdir}/assembly/{sample}_{assembler}.done", tmpdir=[config['tmpdir']], sample=SAMPLES, assembler=[config['assembler']])
    output:
        f"{config['tmpdir']}/successful_samples.txt"
    message: "Evaluate which samples were successfully assembled with at least one contig"
    run:
        with open(output[0], "wt") as outfile:
            for s in SAMPLES:
                with gzip.open(f"{config['resultdir']}/assembly/{s}-{config['assembler']}.fa.gz", "rt") as fastafile:
                    line = fastafile.readline()
                    if line.startswith(">"):
                        outfile.write(s + "\n")

#### Error correction with SPAdes-hammer ######################################

checkpoint error_correction:
    input:
        "{tmpdir}/sampletsv.validated"
    output:
        "{tmpdir}/error_correction/spadeshammer_{sample}/corrected/corrected.yaml"
    message: "Correct reads of sample {wildcards.sample} using SPAdes hammer"
    conda: "../envs/ENVS_metaSPAdes.yaml"
    group: "spadeshammer"
    resources:
        mem = config['assembly_mem'],
        cores = 18
    params: 
        pe1 = lambda wildcards: sampletsv.at[wildcards.sample, 'R1'],
        pe2 = lambda wildcards: sampletsv.at[wildcards.sample, 'R2'],
        pe0 = lambda wildcards: f"-s {sampletsv.at[wildcards.sample, 'R0']}" if sampletsv.at[wildcards.sample, 'R0'] != "NA" else "",
        fileprefix = lambda wildcards: os.path.basename(sampletsv.at[wildcards.sample, 'R1']).split("_")[0],
        filesuffix = lambda wildcards: sampletsv.at[wildcards.sample, 'R1'].split(".")[-2],
        memory = int(config['assembly_mem'] * 0.9),
        outdir = "{tmpdir}/error_correction/spadeshammer_{sample}"
    threads: 18
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.outdir}/tmp
        metaspades.py -o {params.outdir} \
            -1 {params.pe1} \
            -2 {params.pe2} \
            {params.pe0} \
            --tmp-dir {params.outdir}/tmp \
            --threads {threads} \
            --memory {params.memory} \
            --only-error-correction
        """

rule rename_spadeshammer_fastqs:
    input:
        "{tmpdir}/error_correction/spadeshammer_{sample}/corrected/corrected.yaml"
    output:
        pe1 = "{tmpdir}/error_correction/{sample}-wreadcorr_1.fastq.gz",
        pe2 = "{tmpdir}/error_correction/{sample}-wreadcorr_2.fastq.gz"
    message: "Rename the corrected FastQs: {wildcards.sample}"
    group: "spadeshammer"
    params:
        singleend_data = lambda wildcards: sampletsv.at[wildcards.sample, 'R0'] != "NA",
        pe1 = lambda wildcards: path_to_corrected_r(wildcards, 1),
        pe2 = lambda wildcards: path_to_corrected_r(wildcards, 2),
        pe0 = lambda wildcards: path_to_corrected_r(wildcards, 0),
        output_pe0 = "{tmpdir}/error_correction/{sample}-wreadcorr_0.fastq.gz"
    shell:
        """
        mv {params.pe1} {output.pe1}
        mv {params.pe2} {output.pe2}
        if [[ "{params.singleend_data}" = "True" ]]; then
            cat {params.pe0} > {params.output_pe0}
        fi
        """

##### Assembly ##################################################################

if config['assembler'] == "megahit":

    rule megahit:
        input:
            sample_validation = "{tmpdir}/sampletsv.validated",
            pe1 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_1.fastq.gz" if config['readcorrection'] else [],
            pe2 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_2.fastq.gz" if config['readcorrection'] else []
        output:
            "{tmpdir}/assembly/megahit/{sample}/final.contigs.fa"
        message: "Assemble the metagenomic data using MEGAHIT: {wildcards.sample}"
        conda: "../envs/ENVS_MEGAHIT.yaml"
        resources:
            mem = lambda wildcards: int(config['assembly_mem'] / 2),
            cores = 24,
            assembly = 1
        params: 
            pe1 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "R1", '1'),
            pe2 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "R2", '2'),
            pe0 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "R0", 'r'),
            memory = f"{int(config['assembly_mem'] / 2)}00000000",
            prefix = "{tmpdir}/assembly/megahit/{sample}",
            tmpdir = "tmp",
            minlength = lambda wildcards: config['min_contiglength']
        threads: 24
        shell:
            """
            if [[ -d {params.prefix} ]]; then
                rm -r {params.prefix}
            fi
            megahit \
                -t {threads} \
                {params.pe1} \
                {params.pe2} \
                {params.pe0} \
                --tmp-dir {params.tmpdir} \
                --min-contig-len {params.minlength} \
                --out-dir {params.prefix}
            """

    checkpoint cleanup_megahit:
        input:
            "{tmpdir}/assembly/megahit/{sample}/final.contigs.fa"
        output:
            "{tmpdir}/assembly/{sample}_megahit.done"
        message: "Clean up assembly folder of MEGAHIT: {wildcards.sample}"
        conda: "../envs/ENVS_unixEssentials.yaml"
        resources:
            mem = 4,
            cores = 4
        params:
            dir = "{tmpdir}/assembly/megahit/{sample}",
            tar = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-megahit.tar.gz",
            fasta = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-megahit.fa.gz",
            logfile = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-megahit.log"
        threads: 4
        shell:
            """
            mkdir -p $(dirname {params.tar})
            tar -czvf {params.tar} {params.dir}/intermediate_contigs/
            pigz -p {threads} -c {input} > {params.fasta}  
            cp {params.dir}/log {params.logfile}
            touch {output}
            """

elif config['assembler'] == "metaspades":

    rule metaspades:
        input:
            sample_validation = "{tmpdir}/sampletsv.validated",
            pe1 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_1.fastq.gz" if config['readcorrection'] else [],
            pe2 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_2.fastq.gz" if config['readcorrection'] else []
        output:
            "{tmpdir}/assembly/metaspades/{sample}/contigs.fasta"
        message: "Assemble the metagenomic data with metaSPAdes: {wildcards.sample}"
        conda: "../envs/ENVS_metaSPAdes.yaml"
        resources:
            mem = config['assembly_mem'],
            cores = 24,
            assembly = 1
        params: 
            pe1 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "R1", '1'),
            pe2 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "R2", '2'),
            pe0 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "R0", 's'),
            outdir = "{tmpdir}/assembly/metaspades/{sample}",
            memory = config['assembly_mem'],
            kmers = config['metaspades_kmers']
        threads: 24
        shell:
            """
            metaspades.py -o {params.outdir} \
                {params.pe1} \
                {params.pe2} \
                {params.pe0} \
                -k {params.kmers} \
                --threads {threads} \
                --tmp-dir {params.outdir}/tmp \
                --memory {params.memory} \
                --only-assembler
            """

    checkpoint cleanup_metaspades:
        input:
            "{tmpdir}/assembly/metaspades/{sample}/contigs.fasta"
        output:
            "{tmpdir}/assembly/{sample}_metaspades.done"
        message: "Clean up assembly folder of metaSPAdes: {wildcards.sample}"
        conda: "../envs/ENVS_unixEssentials.yaml"
        resources:
            mem = 4,
            cores = 4
        params:
            dir = "{tmpdir}/assembly/metaspades/{sample}",
            scaffolds = "{tmpdir}/assembly/metaspades/{sample}/scaffolds.fasta", 
            tar = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-metaspades.tar.gz",
            fasta = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-metaspades.fa.gz",
            logfile = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-metaspades.log"
        threads: 4
        shell:
            """
            mkdir -p $(dirname {params.tar})
            tar -czvf {params.tar} {params.dir}/{{assembly_graph.fastg,before_rr.fasta,contigs.fasta,contigs.paths,scaffolds.fasta,scaffolds.paths}}
            pigz -p {threads} -c {params.scaffolds} > {params.fasta}  
            cp {params.dir}/spades.log {params.logfile} 
            touch {output}
            """
