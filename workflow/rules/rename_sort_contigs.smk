rule rename_sort_contigs:
    input:
        lambda wildcards: expand("{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards))
    output:
        touch(f"{config['tmpdir']}/rename_sort_contigs.done")


if config['assembler'] == "megahit":

    rule average_depth:
        input:
            "{tmpdir}/alignment/{assembler}/{sample}.samtools_depth"
        output:
            temp("{tmpdir}/alignment/{assembler}/{sample}.avg_depth")
        message: "Infer mean coverage per contig: {wildcards.sample}"
        resources:
            mem = 4,
            cores = 1
        run:
            with open(output[0], "wt") as outfile:
                outfile.write("contig\tdepth\n")
                prev_contig = ""
                depths = [] 
                for line in open(input[0], "rt"):
                    contig, pos, depth = line.rstrip().split()
                    if prev_contig != "" and contig != prev_contig:
                        avg_depth = sum(depths) / len(depths)
                        outfile.write(f"{prev_contig}\t{avg_depth:.2f}\n")
                        depths = []
                    prev_contig = contig
                    depths.append(int(depth))
                avg_depth = sum(depths) / len(depths)
                outfile.write(f"{contig}\t{avg_depth:.2f}\n")

    rule fix_contignames:
        input:
            fasta = "{resultdir}/consensus_correction/{assembler}/{sample}_contigs.fasta.gz",
            bam = lambda wildcards: f"{config['tmpdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.sorted.noncorr.bam",
            depth = lambda wildcards: f"{config['tmpdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.avg_depth"
        output:
            fasta = "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
        message: "Fix contig names of the MEGAHIT to make them look like the metaSPAdes ones: {wildcards.sample}"
        conda: "../envs/ENVS_Python.yaml"
        resources:
            mem = 4,
            cores = 1
        params:
            header = lambda wildcards: f"{config['tmpdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.header"
        script:
            "../scripts/megahit2metaspades_contignames.py"

    rule samtools_reheader:
        input:
            bam = "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam",
            bai = "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam.bai",
            fasta = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
        output:
            bam = temp("{tmpdir}/alignment/{assembler}/{sample}.sorted.renamed.bam"),
            bai = temp("{tmpdir}/alignment/{assembler}/{sample}.sorted.renamed.bam.bai"),
        message: "Replace the header of the BAM file with the ones with the correct contig names: {wildcards.sample}"
        conda: "../envs/ENVS_samtools.yaml"
        resources:
            mem = 4,
            cores = 1
        params:
            header = "{tmpdir}/alignment/{assembler}/{sample}.header"
        shell:
            """
            samtools reheader {params.header} {input.bam} > {output.bam}
            samtools index {output.bam}
            """

elif config['assembler'] == "metaspades":

    rule copy_fasta:
        input:
            "{resultdir}/consensus_correction/{assembler}/{sample}_contigs.fasta.gz"
        output:
            "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
        message: "Copy FastA file: {wildcards.sample}"
        conda: "../envs/ENVS_samtools.yaml"
        params:
            prefix = lambda wildcards: f"{os.getcwd()}/" if wildcards.resultdir[0] != "/" else "" 
        shell:
            """
            ln -s {params.prefix}{input} {output}
            """

    rule link_bam:
        input:
            bam = "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam",
            bai = "{tmpdir}/alignment/{assembler}/{sample}.sorted.noncorr.bam.bai",
            fasta = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
        output:
            bam = temp("{tmpdir}/alignment/{assembler}/{sample}.sorted.renamed.bam"),
            bai = temp("{tmpdir}/alignment/{assembler}/{sample}.sorted.renamed.bam.bai"),
        message: "Link the sorted BAMs and copy the FastA with the metaSPAdes contigs: {wildcards.sample}" 
        resources:
            mem = 4,
            cores = 1
        params:
            prefix = lambda wildcards: f"{os.getcwd()}/" if wildcards.tmpdir[0] != "/" else "" 
        shell:
            """
            ln -s {params.prefix}/{input.bam} {output.bam}
            ln -s {params.prefix}/{input.bai} {output.bai}
            """
