rule samtools_sort_by_name:
    input:
        bam = "{tmpdir}/alignment/{assembler}/{sample}.sorted.renamed.bam",
        bai = "{tmpdir}/alignment/{assembler}/{sample}.sorted.renamed.bam.bai"
    output:
        pipe("{tmpdir}/alignment/{assembler}/{sample}.nsorted.bam")
    message: "Sort BAM file by name: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    group: "samtools_fixmate"
    resources:
        mem = 16,
        cores = 4
    threads: 4
    shell:
        """
        samtools sort -n -l 0 -@ {threads} -o {output} {input.bam}
        """

rule samtools_fixmate:
    input:
        "{tmpdir}/alignment/{assembler}/{sample}.nsorted.bam"
    output:
        pipe("{tmpdir}/alignment/{assembler}/{sample}.fixmate.bam")
    message: "Apply samtools fixmate: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    group: "samtools_fixmate"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        samtools fixmate -mc -@ {threads} {input} {output}
        """

rule samtools_sort_by_coord:
    input:
        "{tmpdir}/alignment/{assembler}/{sample}.fixmate.bam"
    output:
        temp("{tmpdir}/alignment/{assembler}/{sample}.fixmate.sorted.bam")
    message: "Sort BAM file back to coordinates: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    group: "samtools_fixmate"
    resources:
        mem = 16,
        cores = 4
    threads: 4
    shell:
        """
        samtools sort -l 0 -@ {threads} -o {output} {input}
        """

rule samtools_markup:
    input:
        bam = lambda wildcards: f"{config['tmpdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.fixmate.sorted.bam"
    output:
        "{resultdir}/alignment/{assembler}/{sample}.sorted.dedup.bam"
    message: "Mark duplicate reads: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    resources:
        mem = 16,
        cores = 4
    threads: 4
    log: "{resultdir}/stats/markdup/{sample}-{assembler}_samtoolsmarkdup.log"
    shell:
        """
        samtools markdup -r -s -@ {threads} {input.bam} {output} 2> {log}
        """

rule samtools_flagstat:
    input:
        "{resultdir}/alignment/{assembler}/{sample}.sorted.dedup.bam"
    output:
        "{resultdir}/stats/flagstat/{sample}-{assembler}.flagstat"
    message: "Samtools flagstat: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule samtools_index:
    input:
        bam = "{resultdir}/alignment/{assembler}/{sample}.sorted.dedup.bam",
        flagstat = "{resultdir}/stats/flagstat/{sample}-{assembler}.flagstat"
    output:
        "{resultdir}/alignment/{assembler}/{sample}.sorted.dedup.bam.bai"
    message: "Index deduplicated BAM file: {wildcards.sample}"
    conda: "../envs/ENVS_samtools.yaml"
    resources:
        mem = 4,
        cores = 1
    shell:
        """
        samtools index {input.bam}
        """
