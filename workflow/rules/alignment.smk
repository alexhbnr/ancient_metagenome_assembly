import numpy as np
import pandas as pd

MINLENGTH = config['min_contiglength']

rule alignment_workflow:
    input: 
        expand("{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam.bai", tmpdir=[config['tmpdir']], assembler=[config['assembler']], sample=SAMPLES)
    output:
        touch(f"{config['tmpdir']}/alignment.done")

if config['assembler'] == "metaspades":

    rule filter_by_contig_length:
        input:
            "{tmpdir}/assembly/{sample}_metaspades.done"
        output:
            temp("{tmpdir}/alignment/{assembler}/{sample}.raw.fasta")
        message: "Filter contigs for a minimal length of >= {MINLENGTH} bp: {wildcards.sample}"
        conda: "../envs/ENVS_bioawk.yaml"
        resources:
            mem = 4,
            cores = 1
        params:
            fasta = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-metaspades.fa.gz"
        threads: 1
        shell:
            """
            bioawk -c fastx '{{if (length($seq) >= {MINLENGTH}){{print ">"$name; print $seq}}}}' {params.fasta} > {output}
            """

elif config['assembler'] == "megahit":

    rule link_contigs:
        input:
            "{tmpdir}/assembly/{sample}_megahit.done"
        output:
            temp("{tmpdir}/alignment/{assembler}/{sample}.raw.fasta")
        message: "Link FastA file with contigs from MEGAHIT: {wildcards.sample}"
        conda: "../envs/ENVS_unixEssentials.yaml"
        resources:
            mem = 4,
            cores = 1
        params:
            fasta = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-megahit.fa.gz"
        threads: 1
        shell:
            """
            pigz -d -c {params.fasta} > {output}
            """

rule build_BowTie2_index:
    input:
        "{tmpdir}/alignment/{assembler}/{sample}.raw.fasta"
    output:
        index1 = temp("{tmpdir}/alignment/{assembler}/{sample}.1.bt2"),
        index2 = temp("{tmpdir}/alignment/{assembler}/{sample}.2.bt2"),
        index3 = temp("{tmpdir}/alignment/{assembler}/{sample}.3.bt2"),
        index4 = temp("{tmpdir}/alignment/{assembler}/{sample}.4.bt2"),
        rev_index1 = temp("{tmpdir}/alignment/{assembler}/{sample}.rev.1.bt2"),
        rev_index2 = temp("{tmpdir}/alignment/{assembler}/{sample}.rev.2.bt2")
    message: "Index the contigs for alignment using BowTie2: {wildcards.sample}"
    conda: "../envs/ENVS_bowtie2.yaml"
    resources:
        mem = 8,
        cores = 4
    params:
        index = "{tmpdir}/alignment/{assembler}/{sample}"
    threads: 4
    shell:
        """
        bowtie2-build -f \
            {input} \
            {params.index}
        """

rule BowTie2_alignment:
    input:
        index1 = "{tmpdir}/alignment/{assembler}/{sample}.1.bt2",
        index2 = "{tmpdir}/alignment/{assembler}/{sample}.2.bt2",
        index3 = "{tmpdir}/alignment/{assembler}/{sample}.3.bt2",
        index4 = "{tmpdir}/alignment/{assembler}/{sample}.4.bt2",
        rev_index1 = "{tmpdir}/alignment/{assembler}/{sample}.rev.1.bt2",
        rev_index2 = "{tmpdir}/alignment/{assembler}/{sample}.rev.2.bt2"
    output:
        pipe("{tmpdir}/alignment/{assembler}/{sample}.sorted.raw.sam")
    message: "Align reads back to the uncorrected contigs using BowTie2's very-sensitive setting: {wildcards.sample}"
    conda: "../envs/ENVS_bowtie2.yaml"
    resources:
        mem = 16,
        cores = 16
    params:
        index = "{tmpdir}/alignment/{assembler}/{sample}",
        n_mismatches = lambda wildcards: config['bowtie2_seed_nmismatches'],
        pe1 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_1.fastq.gz" if config['readcorrection'] else sampletsv.at[wildcards.sample, 'R1'],
        pe2 = lambda wildcards: f"{wildcards.tmpdir}/error_correction/{wildcards.sample}-wreadcorr_2.fastq.gz" if config['readcorrection'] else sampletsv.at[wildcards.sample, 'R2'],
        pe0 = lambda wildcards: f"-U {path_to_r0(wildcards.sample, wildcards.tmpdir)}" if path_to_r0(wildcards.sample, wildcards.tmpdir) != "" else "",
    threads: 16
    shell:
        """
        bowtie2 -p {threads} --very-sensitive -N {params.n_mismatches} -x {params.index} \
                -1 {params.pe1} -2 {params.pe2} {params.pe0} -S {output}
        """

rule samtools_sort:
    input:
        bam = "{tmpdir}/alignment/{assembler}/{sample}.sorted.raw.sam",
        fasta = "{tmpdir}/alignment/{assembler}/{sample}.raw.fasta"
    output:
        bam = "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam",
        bai = "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam.bai"
    message: "Sort the sequencing data: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 2
    threads: 2
    shell:
        """
        samtools view -Sb {input.bam} | \
        samtools calmd -u /dev/stdin {input.fasta} | \
        samtools sort -l 4 -o {output.bam} -
        samtools index {output.bam}
        """

rule samtools_depth:
    input:
        "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam"
    output:
        temp("{tmpdir}/alignment/{assembler}/{sample}.samtools_depth")
    message: "Determine the depth along the contigs with samtools: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 1
    shell:
        """
        samtools depth -a {input} > {output}
        """
