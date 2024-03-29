import pandas as pd

def path_to_nonudg_r0(sample, tmpdir):
    if sampletsv.at[sample, 'nonUDG_R0'] == "NA":  # no single-end data
        return ""
    else:
        return sampletsv.at[sample, 'nonUDG_R0']

if config['pydamage']:

    rule build_BowTie2_index_pydamage:
        input:
            lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
        output:
            index1 = temp("{tmpdir}/pydamage/{assembler}/{sample}.1.bt2"),
            index2 = temp("{tmpdir}/pydamage/{assembler}/{sample}.2.bt2"),
            index3 = temp("{tmpdir}/pydamage/{assembler}/{sample}.3.bt2"),
            index4 = temp("{tmpdir}/pydamage/{assembler}/{sample}.4.bt2"),
            rev_index1 = temp("{tmpdir}/pydamage/{assembler}/{sample}.rev.1.bt2"),
            rev_index2 = temp("{tmpdir}/pydamage/{assembler}/{sample}.rev.2.bt2")
        message: "Index the corrected contigs for alignment using BowTie2: {wildcards.sample}"
        conda: "../envs/ENVS_bowtie2.yaml"
        resources:
            mem = 8,
            cores = 4
        params:
            index = "{tmpdir}/pydamage/{assembler}/{sample}"
        threads: 4
        shell:
            """
            bowtie2-build -f \
                {input} \
                {params.index}
            """

    rule BowTie2_alignment_pydamage:
        input:
            index1 = "{tmpdir}/pydamage/{assembler}/{sample}.1.bt2",
            index2 = "{tmpdir}/pydamage/{assembler}/{sample}.2.bt2",
            index3 = "{tmpdir}/pydamage/{assembler}/{sample}.3.bt2",
            index4 = "{tmpdir}/pydamage/{assembler}/{sample}.4.bt2",
            rev_index1 = "{tmpdir}/pydamage/{assembler}/{sample}.rev.1.bt2",
            rev_index2 = "{tmpdir}/pydamage/{assembler}/{sample}.rev.2.bt2"
        output:
            pipe("{tmpdir}/pydamage/{assembler}/{sample}.sorted.raw.sam")
        message: "Align the non-UDG reads back to the corrected contigs using BowTie2's very-sensitive setting: {wildcards.sample}"
        conda: "../envs/ENVS_bowtie2.yaml"
        group: "ref_alignment_pyd"
        resources:
            mem = 16,
            cores = 16
        params:
            index = "{tmpdir}/pydamage/{assembler}/{sample}",
            n_mismatches = lambda wildcards: config['bowtie2_seed_nmismatches'],
            pe1 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "nonUDG_R1", '1'),
            pe2 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "nonUDG_R2", '2'),
            pe0 = lambda wildcards: path_to_r(wildcards.sample, wildcards.tmpdir, "nonUDG_R0", 'U')
        threads: 16
        shell:
            """
            bowtie2 -p {threads} --very-sensitive -N {params.n_mismatches} -x {params.index} \
                    {params.pe1} {params.pe2} {params.pe0} -S {output}
            """

    rule samtools_sort_pydamage:
        input:
            bam = "{tmpdir}/pydamage/{assembler}/{sample}.sorted.raw.sam",
            fasta = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
        output:
            bam = temp("{tmpdir}/pydamage/{assembler}/{sample}.sorted.pydamage.bam"),
            bai = temp("{tmpdir}/pydamage/{assembler}/{sample}.sorted.pydamage.bam.bai")
        message: "Sort the non-UDG sequencing data: {wildcards.sample}"
        conda: "../envs/ENVS_samtools.yaml"
        group: "ref_alignment_pyd"
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

    rule pyDamage:
        input:
            bam = lambda wildcards: f"{config['tmpdir']}/pydamage/{wildcards.assembler}/{wildcards.sample}.sorted.pydamage.bam",
            bai = lambda wildcards: f"{config['tmpdir']}/pydamage/{wildcards.assembler}/{wildcards.sample}.sorted.pydamage.bam.bai"
        output:
            "{resultdir}/pydamage/{sample}-{assembler}.pydamage.csv.gz"
        message: "Analyse aDNA damage using PyDamage: {wildcards.sample}"
        conda: "../envs/ENVS_pydamage.yaml"
        resources:
            mem = 12,
            cores = 8
        params:
            outputprefix = "{resultdir}/pydamage/{sample}-{assembler}.pydamage",
            pydamage_out = "{resultdir}/pydamage/{sample}-{assembler}.pydamage.csv",
        threads: 8
        shell:
            """
            pydamage -o {params.outputprefix} analyze -w 35 \
                    -p {threads} \
                    --force \
                    {input.bam}
            mv {params.outputprefix}/pydamage_results.csv {params.pydamage_out}
            gzip -f {params.pydamage_out}
            rmdir {params.outputprefix}
            """

else:

    localrules: pyDamage_dummy_output

    rule pyDamage_dummy_output:
        output:
            "{resultdir}/pydamage/{sample}-{assembler}.pydamage.csv.gz"
        message: "Create empty dummy files for pyDamage analysis: {wildcards.sample}"
        shell:
            """
            touch {output}
            """

        

