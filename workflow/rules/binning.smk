if config['magbinning']:

	if config['magbinning_type'] == "metawrap":

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
