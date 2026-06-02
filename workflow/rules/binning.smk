if config['magbinning']:

    rule binning_workflow:
        input: 
            lambda wildcards: expand("{resultdir}/binning/{sample}-{assembler}_binning.done", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards))
        output:
            touch(f"{config['tmpdir']}/binning.done")

    rule create_sequence_dict:
        input:
            lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
        output:
            temp("{tmpdir}/binning/{sample}-{assembler}.dict")
        message: "Create sequence dictionary to reorder the BAM file: {wildcards.sample}"
        container: "https://depot.galaxyproject.org/singularity/picard%3A3.4.0--hdfd78af_0"
        resources:
            mem_gb = 4,
            mem_mb = 8000,
            runtime = 120,
            slurm_partition = "short"
        params:
            java_opts = ""
        threads: 1
        shell:
            """
            picard CreateSequenceDictionary \
                {params.java_opts} \
                --REFERENCE {input} \
                --OUTPUT {output}
            """

    if config['assembler'] == "megahit":

        rule reorder_sam:
            input:
                bam = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.sorted.dedup.bam",
                seq_dict = "{tmpdir}/binning/{sample}-{assembler}.dict"
            output:
                temp("{tmpdir}/binning/{sample}-{assembler}.reorder.bam")
            message: "Sort BAM file according to FastA file: {wildcards.sample}"
            container: "https://depot.galaxyproject.org/singularity/picard%3A3.4.0--hdfd78af_0"
            resources:
                mem_gb = lambda wildcards, attempt: 32 + attempt * 40,
                mem_mb = lambda wildcards, attempt: 48000 + attempt + 48000,
                runtime = 1440,
                slurm_partition = "standard"
            params:
                java_opts = ""
            threads: 1
            shell:
                """
                picard ReorderSam \
                    {params.java_opts} \
                    --INPUT {input.bam} \
                    --SEQUENCE_DICTIONARY {input.seq_dict} \
                    --OUTPUT {output}
                """
        
        rule index_reordered_bam:
            input:
                "{tmpdir}/binning/{sample}-{assembler}.reorder.bam"
            output:
                temp("{tmpdir}/binning/{sample}-{assembler}.reorder.bam.bai")
            message: "Index reordered BAM file: {wildcards.sample}"
            container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
            resources:
                mem_mb = 8000,
                runtime = 120,
                slurm_partition = "short"
            threads: 4
            shell:
                "samtools index {input}"

    elif config['assembler'] == "metaspades":

        localrules: link_bam_binning

        rule link_bam_binning:
            input:
                bam = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.sorted.dedup.bam",
                bai = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.sorted.dedup.bam.bai"
            output:
                bam = temp("{tmpdir}/binning/{sample}-{assembler}.reorder.bam"),
                bai = temp("{tmpdir}/binning/{sample}-{assembler}.reorder.bam.bai")
            message: "Link BAM files for binning: {wildcards.sample}"
            params:
                prefix = lambda wildcards: f"{os.getcwd()}" if not config['resultdir'].startswith("/") else ""
            shell:
                """
                ln -s {params.prefix}/{input.bam} {output.bam}
                ln -s {params.prefix}/{input.bai} {output.bai}
                """

    if config['metabat2']:

        rule metabat2_depth:
            input:
                "{tmpdir}/binning/{sample}-{assembler}.reorder.bam"
            output:
                "{tmpdir}/binning/metabat2/{sample}-{assembler}.depth"
            message: "Determine the average depth of the contigs for metaBAT2: {wildcards.sample}"
            container: "docker://quay.io/biocontainers/metabat2:2.18--h38e344b_2"
            resources:
                mem_mb = 8000,
                runtime = 1440,
                slurm_partition = "standard"
            params:
                minlength = config['min_binninglength']
            shell:
                """
                jgi_summarize_bam_contig_depths \
                    --outputDepth {output} \
                    --minContigDepth 1 \
                    --minContigLength {params.minlength} \
                    {input}
                """

        rule metabat2:
            input:
                fasta = lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz",
                depth = "{tmpdir}/binning/metabat2/{sample}-{assembler}.depth"
            output:
                "{tmpdir}/binning/metabat2/{sample}-{assembler}.unbinned.fa"
            message: "Bin the de-novo assembled contigs using metaBAT2: {wildcards.sample}"
            container: "docker://quay.io/biocontainers/metabat2:2.18--h38e344b_2"
            resources:
                mem_mb = 24000,
                runtime = 2880,
                slurm_partition = "standard"
            params:
                minlength = lambda wildcards: config['min_binninglength'] if config['min_binninglength'] > 1500 else 1500,
                outprefix = "{tmpdir}/binning/metabat2/{sample}-{assembler}"
            threads: 8
            shell:
                """
                if [[ $(cat {input.depth} | wc -l) -gt 1 ]]; then
                    metabat2 \
                        -i {input.fasta} \
                        -o {params.outprefix} \
                        -a {input.depth} \
                        --minContig {params.minlength} \
                        --minClsSize 100000 \
                        -t {threads} \
                        --unbinned --verbose
                else
                    gunzip -c {input.fasta} > {output}
                fi
                """

        rule summarise_metabat2:
            input:
                lambda wildcards: f"{config['tmpdir']}/binning/metabat2/{wildcards.sample}-{wildcards.assembler}.unbinned.fa"
            output:
                "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/metabat2_bins/bin.unbinned.fa"
            message: "Copy bins into the folder structure expected by MetaWRAP: {wildcards.sample}"
            resources:
                mem_mb = 2000,
                runtime = 120,
                slurm_partition = "short"
            params:
                metawrapdir = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/metabat2_bins"
            run:
                for fn in glob(f"{input[0].replace('.unbinned.fa', '')}.*.fa"):
                    if fn.split(".")[-2] == "unbinned" or fn.split(".")[-2].isnumeric():
                        shutil.copy(fn, f"{params.metawrapdir}/bin.{fn.split('.')[-2]}.fa")

    if config['maxbin2']:

        localrules: format_maxbin2_depth, decompress_fasta_maxbin2

        rule maxbin2_depth:
            input:
                "{tmpdir}/binning/{sample}-{assembler}.reorder.bam"
            output:
                temp("{tmpdir}/binning/maxbin2/{sample}-{assembler}.metabat2_depth")
            message: "Determine the average depth of the contigs for MaxBin2: {wildcards.sample}"
            container: "docker://quay.io/biocontainers/metabat2:2.18--h38e344b_2"
            resources:
                mem_mb = 8000,
                runtime = 1440,
                slurm_partition = "standard"
            params:
                minlength = config['min_binninglength']
            shell:
                """
                jgi_summarize_bam_contig_depths \
                    --outputDepth {output} \
                    --minContigDepth 1 \
                    --minContigLength {params.minlength} \
                    {input}
                """

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
                temp("{tmpdir}/binning/maxbin2/{sample}-{assembler}.fa")
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
            container: "docker://quay.io/biocontainers/maxbin2:2.2.7--h503566f_8"
            resources:
                mem_mb = 16000,
                runtime = 2880,
                slurm_partition = "standard"
            params:
                minlength = config['min_binninglength'],
                markerset = lambda wildcards: config['maxbin2_markerset'],
                outprefix = "{tmpdir}/binning/maxbin2/{sample}-{assembler}"
            threads: 8
            shell:
                """
                run_MaxBin.pl \
                    -contig {input.fasta} \
                    -abund {input.depth} \
                    -out {params.outprefix} \
                    -markerset {params.markerset} \
                    -min_contig_length {params.minlength} \
                    -thread {threads} || \
                if [[ $? -eq 255 ]]; then
                    touch {output}
                fi
                """

        rule summarise_maxbin2:
            input:
                lambda wildcards: f"{config['tmpdir']}/binning/maxbin2/{wildcards.sample}-{wildcards.assembler}.noclass"
            output:
                "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/maxbin2_bins/{sample}-{assembler}.summary"
            message: "Copy bins into the folder structure expected by MetaWRAP: {wildcards.sample}"
            resources:
                mem_mb = 2000,
                runtime = 120,
                slurm_partition = "short"
            params:
                metawrapdir = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/maxbin2_bins"
            run:
                for fn in glob(f"{input[0].replace('.noclass', '')}.*.fasta"):
                    shutil.copy(fn, f"{params.metawrapdir}/bin.{int(fn.split('.')[-2]) - 1}.fa")
                if os.path.isfile(input[0].replace('.noclass', '.log')):
                    shutil.copy(input[0].replace('.noclass', '.log'),
                                f"{params.metawrapdir}/{wildcards.sample}-{wildcards.assembler}.log")
                else:
                    Path(f"{params.metawrapdir}/{wildcards.sample}-{wildcards.assembler}.log").touch()
                if os.path.isfile(input[0].replace('.noclass', '.summary')):
                    shutil.copy(input[0].replace('.noclass', '.summary'),
                                f"{params.metawrapdir}/{wildcards.sample}-{wildcards.assembler}.summary")
                else:
                    Path(output[0]).touch()

    if config['concoct']:

        localrules: decompress_fasta_concot

        rule decompress_fasta_concot:
            input:
                lambda wildcards: f"{config['resultdir']}/alignment/{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.fasta.gz"
            output:
                temp("{tmpdir}/binning/concoct/{sample}-{assembler}.fa")
            message: "Decompress the FastA file with the contigs for CONCOCT: {wildcards.sample}"
            shell:
                "gunzip -c {input} > {output}"

        rule concoct_split_contigs:
            input:
                "{tmpdir}/binning/concoct/{sample}-{assembler}.fa"
            output:
                fa = temp("{tmpdir}/binning/concoct/{sample}-{assembler}_10K.fa"),
                bed = temp("{tmpdir}/binning/concoct/{sample}-{assembler}_10K.bed")
            message: "Split the contigs into chunks of max. 10 kb: {wildcards.sample}"
            container: 'docker://quay.io/biocontainers/concoct:1.1.0--py312hb1d17a5_9'
            resources:
                mem_mb = 4000,
                runtime = 120,
                slurm_partition = "short"
            shell:
                """
                cut_up_fasta.py \
                    {input} -c 10000 --merge_last \
                    -b {output.bed} -o 0 > {output.fa}
                """

        rule concoct_coverage_table:
            input:
                bed = "{tmpdir}/binning/concoct/{sample}-{assembler}_10K.bed",
                bam = "{tmpdir}/binning/{sample}-{assembler}.reorder.bam",
                bai = "{tmpdir}/binning/{sample}-{assembler}.reorder.bam.bai"
            output:
                "{tmpdir}/binning/concoct/{sample}-{assembler}.depth"
            message: "Calculate the depth along the contigs: {wildcards.sample}"
            container: 'docker://quay.io/biocontainers/concoct:1.1.0--py312hb1d17a5_9'
            resources:
                mem_mb = 4000,
                runtime = 1440,
                slurm_partition = "standard"
            shell:
                """
                concoct_coverage_table.py {input.bed} {input.bam} > {output}
                """

        rule concoct:
            input:
                fa = "{tmpdir}/binning/concoct/{sample}-{assembler}_10K.fa",
                depth = "{tmpdir}/binning/concoct/{sample}-{assembler}.depth"
            output:
                "{tmpdir}/binning/concoct/{sample}-{assembler}/{sample}-{assembler}_args.txt"
            message: "Bin the de-novo assembled contigs using CONCOCT: {wildcards.sample}"
            container: 'docker://quay.io/biocontainers/concoct:1.1.0--py312hb1d17a5_9'
            resources:
                mem_mb = lambda wildcards, attempt: 8000 + attempt * 8000,
                runtime = 2880,
                slurm_partition = "standard"
            params:
                minlength = config['min_binninglength'],
                outprefix = "{tmpdir}/binning/concoct/{sample}-{assembler}/{sample}-{assembler}"
            threads: 8
            shell:
                """
                concoct \
                    -l {params.minlength} \
                    -t {threads} \
                    --coverage_file {input.depth} \
                    --composition_file {input.fa} \
                    -b {params.outprefix} || \
                if [[ $? -eq 255 ]]; then
                    touch {output}
                fi
                """

        rule concoct_merge_contigs:
            input:
                "{tmpdir}/binning/concoct/{sample}-{assembler}/{sample}-{assembler}_args.txt"
            output:
                "{tmpdir}/binning/concoct/{sample}-{assembler}/{sample}-{assembler}_clustering.csv"
            message: "Merge the split contigs: {wildcards.sample}"
            container: 'docker://quay.io/biocontainers/concoct:1.1.0--py312hb1d17a5_9'
            resources:
                mem_mb = 4000,
                runtime = 120,
                slurm_partition = "short"
            params:
                clustering = lambda wildcards: f"{wildcards.tmpdir}/binning/concoct/{wildcards.sample}-{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}_clustering_gt{config['min_binninglength']}.csv"
            shell:
                """
                merge_cutup_clustering.py {params.clustering} > {output}
                """

        rule summarise_concoct:
            # TODO: convert to script and add Python dependecies to container
            input:
                cluster = lambda wildcards: f"{config['tmpdir']}/binning/concoct/{wildcards.sample}-{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}_clustering.csv",
                fasta = "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
            output:
                "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/{sample}-{assembler}.concoct_clustering.csv"
            message: "Copy CONCOCT bins into the folder structure expected by MetaWRAP: {wildcards.sample}"
            conda: "../envs/ENVS_Python.yaml"
            resources:
                mem_mb = 2000,
                runtime = 120,
                slurm_partition = "short"
            params:
                metawrapdir = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/concoct_bins",
                min_contiglength = config['min_binninglength']
            script:
                "../scripts/summarise_concoct.py"


    rule checkpoint_binning:
        input:
            metabat2 = lambda wildcards: f"{config['resultdir']}/binning/metawrap/INITIAL_BINNING/{wildcards.sample}-{wildcards.assembler}/metabat2_bins/bin.unbinned.fa" if config['metabat2'] else [],
            maxbin2 = lambda wildcards: f"{config['resultdir']}/binning/metawrap/INITIAL_BINNING/{wildcards.sample}-{wildcards.assembler}/maxbin2_bins/{wildcards.sample}-{wildcards.assembler}.summary" if config['maxbin2'] else [],
            concoct = lambda wildcards: f"{config['resultdir']}/binning/metawrap/INITIAL_BINNING/{wildcards.sample}-{wildcards.assembler}/{wildcards.sample}-{wildcards.assembler}.concoct_clustering.csv" if config['concoct'] else []
        output:
            "{resultdir}/binning/{sample}-{assembler}_binning.done"
        shell:
            "touch {output}"

else:

    localrules: checkpoint_binning

    rule checkpoint_binning:
        input:
            "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
        output:
            "{resultdir}/binning/{sample}-{assembler}_binning.done"
        shell:
            "touch {output}"
