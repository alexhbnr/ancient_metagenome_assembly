import numpy as np
import pandas as pd

#### Auxilliary functions ######################################################

def path_to_r0(sample, tmpdir):
    # TODO: Test on data that also has single-end data
    if np.isnan(sampletsv.at[sample, 'R0']):  # no single-end data
        return ""
    else:
        if config['readcorrection'] == "true":
            return f"{tmpdir}/error_correction/{sample}-wreadcorr_0.fastq.gz"
        else:
            return sampletsv.at[sample, 'R0']

################################################################################

rule assembly_workflow:
    input:
        f"{config['tmpdir']}/sampletsv.validated",
        expand("{tmpdir}/assembly/{sample}_{assembler}.done", tmpdir=[config['tmpdir']], sample=SAMPLES, assembler=[config['assembler']])
    output:
        touch(f"{config['tmpdir']}/assembly.done")

#### Error correction with SPAdes-hammer ######################################

rule error_correction:
    input:
        "{tmpdir}/sampletsv.validated"
    output:
        pe1 = "{tmpdir}/error_correction/{sample}-wreadcorr_1.fastq.gz",
        pe2 = "{tmpdir}/error_correction/{sample}-wreadcorr_2.fastq.gz",
    message: "Correct reads of sample {wildcards.sample} using SPAdes hammer"
    conda: "../envs/ENVS_metaSPAdes.yaml"
    resources:
        mem = 80,
        cores = 18
    params: 
        pe1 = lambda wildcards: sampletsv.at[wildcards.sample, 'R1'],
        pe2 = lambda wildcards: sampletsv.at[wildcards.sample, 'R2'],
        pe0 = lambda wildcards: f"-s {sampletsv.at[wildcards.sample, 'R0']}" if not np.isnan(sampletsv.at['SPM001', 'R0']) else "",
        filesuffix = lambda wildcards: sampletsv.at[wildcards.sample, 'R1'].split(".")[-2],
        output_pe0 = "{tmpdir}/error_correction/{sample}-wreadcorr_0.fastq.gz",
        singleend_data = lambda wildcards: ~np.isnan(sampletsv.at['SPM001', 'R0']),
        memory = 72,
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
        mv {params.outdir}/corrected/{wildcards.sample}_1.{params.filesuffix}.00.0_0.cor.fastq.gz {output.pe1}
        mv {params.outdir}/corrected/{wildcards.sample}_2.{params.filesuffix}.00.0_0.cor.fastq.gz {output.pe2}
        if [[ "{params.singleend_data}" = "True" ]]; then
            mv {params.outdir}/corrected/{wildcards.sample}__unpaired.00.0_0.cor.fastq.gz {params.output_pe0}
        fi
        rm -r {params.outdir}
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
            cores = 24
        params: 
            pe1 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_1.fastq.gz" if config['readcorrection'] else sampletsv.at[wildcards.sample, 'R1'],
            pe2 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_2.fastq.gz" if config['readcorrection'] else sampletsv.at[wildcards.sample, 'R2'],
            pe0 = lambda wildcards: f"-r {path_to_r0(wildcards.sample, wildcards.tmpdir)}" if path_to_r0(wildcards.sample, wildcards.tmpdir) != "" else "",
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
                -1 {params.pe1} \
                -2 {params.pe2} \
                {params.pe0} \
                --tmp-dir {params.tmpdir} \
                --min-contig-len {params.minlength} \
                --out-dir {params.prefix}
            """

    rule cleanup_megahit:
        input:
            "{tmpdir}/assembly/megahit/{sample}/final.contigs.fa"
        output:
            touch("{tmpdir}/assembly/{sample}_megahit.done"),
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
        params: 
            pe1 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_1.fastq.gz" if config['readcorrection'] else sampletsv.at[wildcards.sample, 'R1'],
            pe2 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_2.fastq.gz" if config['readcorrection'] else sampletsv.at[wildcards.sample, 'R2'],
            pe0 = lambda wildcards: f"-r {path_to_r0(wildcards.sample, wildcards.tmpdir)}" if path_to_r0(wildcards.sample, wildcards.tmpdir) != "" else "",
            outdir = "{tmpdir}/assembly/metaspades/{sample}",
            memory = config['assembly_mem'],
            kmers = config['metaspades_kmers']
        threads: 24
        shell:
            """
            metaspades.py -o {params.outdir} \
                -1 {params.pe1} \
                -2 {params.pe2} \
                {params.pe0} \
                -k {params.kmers} \
                --threads {threads} \
                --tmp-dir {params.outdir}/tmp \
                --memory {params.memory} \
                --only-assembler
            """

    rule cleanup_metaspades:
        input:
            "{tmpdir}/assembly/metaspades/{sample}/contigs.fasta"
        output:
            touch("{tmpdir}/assembly/{sample}_metaspades.done")
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
            """
