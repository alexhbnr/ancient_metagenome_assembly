if config['magbinning']:

    rule create_sequence_dict:
        input:
            lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
        output:
            temp("{tmpdir}/binning/{sample}-{assembler}.dict")
        message: "Create sequence dictionary to reorder the BAM file: {wildcards.sample}"
        resources:
            mem = 8,
            mem_gb = 4
        wrapper:
            "v1.3.2/bio/picard/createsequencedictionary"

    rule reorder_sam:
        input:
            bam = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.sorted.dedup.bam",
            seq_dict = "{tmpdir}/binning/{sample}-{assembler}.dict"
        output:
            temp("{tmpdir}/binning/{sample}-{assembler}.reorder.bam")
        message: "Sort BAM file according to FastA file: {wildcards.sample}"
        resources:
            mem = 48,
            mem_gb = 32
        wrapper:
            "https://github.com/alexhbnr/snakemake-wrappers-public/raw/picard_reordersam/bio/picard/reordersam"
    
    rule index_reordered_bam:
        input:
            "{tmpdir}/binning/{sample}-{assembler}.reorder.bam"
        output:
            temp("{tmpdir}/binning/{sample}-{assembler}.reorder.bam.bai")
        message: "Index reordered BAM file: {wildcards.sample}"
        threads: 4
        wrapper:
            "v1.3.2/bio/samtools/index"

    if config['magbinning_type'] == "metawrap":

        rule prepare_metawrap:
            input:
                fasta = "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz",
                bam = lambda wildcards: f"{config['tmpdir']}/binning/{wildcards.sample}-{wildcards.assembler}.reorder.bam",
                bai = lambda wildcards: f"{config['tmpdir']}/binning/{wildcards.sample}-{wildcards.assembler}.reorder.bam.bai"
            output:
                bam = temp("{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files/{sample}-{assembler}.bam"),
                bai = temp("{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files/{sample}-{assembler}.bam.bai"),
                fasta = temp("{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files/{sample}-{assembler}.fasta"),
                fq1 = temp("{resultdir}/binning/metawrap/faux_reads/{sample}-{assembler}_1.fastq"),
                fq2 = temp("{resultdir}/binning/metawrap/faux_reads/{sample}-{assembler}_2.fastq")
            message: "Moving BAM files and creating false fq files to trick the metaWRAP binning module"
            resources:
                cores = 1
            params:
                prefix = lambda wildcards: f"{os.getcwd()}" if not config['resultdir'].startswith("/") else ""
            shell:
                """
                gunzip -c {input.fasta} > {output.fasta}
                ln -s {params.prefix}/{input.bam} {output.bam}
                ln -s {params.prefix}/{input.bai} {output.bai}
                echo "@" > {output.fq1}
                echo "@" > {output.fq2}
                """

        rule metaWRAP_binning:
            input:
                contigs = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files/{sample}-{assembler}.fasta",
                fq1 = "{resultdir}/binning/metawrap/faux_reads/{sample}-{assembler}_1.fastq",
                fq2 = "{resultdir}/binning/metawrap/faux_reads/{sample}-{assembler}_2.fastq",
                bam = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files/{sample}-{assembler}.bam",
                bai = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files/{sample}-{assembler}.bam.bai",
                original_bam = lambda wildcards: f"{config['tmpdir']}/binning/{wildcards.sample}-{wildcards.assembler}.reorder.bam",
                original_bai = lambda wildcards: f"{config['tmpdir']}/binning/{wildcards.sample}-{wildcards.assembler}.reorder.bam.bai"
            output:
                touch("{resultdir}/magbinning/{sample}-{assembler}_binning.done") 
            message: "Running the metaWRAP binning module on {wildcards.sample}"
            resources:
                mem = 60,
                cores = 16
            params:
                outdir = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}",
                workfiles = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/work_files"
            log: "{resultdir}/binning/metawrap/{sample}-{assembler}.binning.log"
            threads: 16
            wrapper:
                "file:///home/alexander_huebner/github/snakemake-wrappers/bio/metawrap/binning"

    elif config['magbinning_type'] == "individual":

        if config['metabat2']:

            rule metabat2_depth:
                input:
                    "{tmpdir}/binning/{sample}-{assembler}.reorder.bam"
                output:
                    "{tmpdir}/binning/metabat2/{sample}-{assembler}.depth"
                message: "Determine the average depth of the contigs for metaBAT2: {wildcards.sample}"
                resources:
                    mem = 8
                params:
                    minlength = config['min_binninglength']
                wrapper:
                    "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/metabat2/depth"

            rule metabat2:
                input:
                    fasta = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz",
                    depth = "{tmpdir}/binning/metabat2/{sample}-{assembler}.depth"
                output:
                    "{tmpdir}/binning/metabat2/{sample}-{assembler}.unbinned.fa"
                message: "Bin the de-novo assembled contigs using metaBAT2: {wildcards.sample}"
                resources:
                    mem = 24
                params:
                    minlength = lambda wildcards: config['min_binninglength'] if config['min_binninglength'] > 1500 else 1500,
                    outprefix = "{tmpdir}/binning/metabat2/{sample}-{assembler}"
                threads: 8
                wrapper:
                    "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/metabat2/metabat2"

        if config['maxbin2']:

            rule maxbin2_depth:
                input:
                    "{tmpdir}/binning/{sample}-{assembler}.reorder.bam"
                output:
                    temp("{tmpdir}/binning/maxbin2/{sample}-{assembler}.metabat2_depth")
                message: "Determine the average depth of the contigs for MaxBin2: {wildcards.sample}"
                resources:
                    mem = 8
                params:
                    minlength = config['min_binninglength']
                wrapper:
                    "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/metabat2/depth"

            rule format_maxbin2_depth:
                input:
                    "{tmpdir}/binning/maxbin2/{sample}-{assembler}.metabat2_depth"
                output:
                    "{tmpdir}/binning/maxbin2/{sample}-{assembler}.depth"
                message: "Format the depth output for MaxBin2: {wildcards.sample}"
                run:
                    with open(output[0], "wt") as outfile:
                        with open(input[0], "rt") as depthfile:
                            for i, line in enumerate(depthfile):
                                if i > 0:
                                    contig, _, depth = line.split("\t")[:3]
                                    outfile.write(f"{contig}\t{depth}\n")

            rule decompress_fasta_maxbin2:
                input:
                    lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
                output:
                    "{tmpdir}/binning/maxbin2/{sample}-{assembler}.fa"
                message: "Decompress the FastA file with the contigs for MaxBin2: {wildcards.sample}"
                shell:
                    "gunzip -c {input} > {output}"

            rule maxbin2:
                input:
                    fasta = "{tmpdir}/binning/maxbin2/{sample}-{assembler}.fa",
                    depth = "{tmpdir}/binning/maxbin2/{sample}-{assembler}.depth"
                output:
                    "{tmpdir}/binning/maxbin2/{sample}-{assembler}.noclass"
                message: "Bin the de-novo assembled contigs using MaxBin2: {wildcards.sample}"
                resources:
                    mem = 16
                params:
                    minlength = config['min_binninglength'],
                    markerset = lambda wildcards: config['maxbin2_markerset'],
                    outprefix = "{tmpdir}/binning/maxbin2/{sample}-{assembler}"
                threads: 8
                wrapper:
                    "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/maxbin2"
