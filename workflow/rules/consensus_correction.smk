rule consensus_correction:
    input:
        expand("{resultdir}/consensus_correction/{assembler}/{sample}_contigs.fasta.gz", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=SAMPLES)
    output:
        touch(f"{config['tmpdir']}/consensus_correction.done")

if config['assembler'] == "megahit":

    rule uncompress_reffasta:
        input:
            "{tmpdir}/assembly/{sample}_megahit.done"
        output:
            temp("{tmpdir}/contig_correction/{sample}.fasta")
        message: "Decompress the FastA file with the de-novo assembled contigs of MEGAHIT: {wildcards.sample}"
        conda: "../envs/ENVS_unixEssentials.yaml"
        resources:
            mem = 4,
            cores = 2
        params:
            fasta = lambda wildcards: f"{config['resultdir']}/assembly/{wildcards.sample}-megahit.fa.gz"
        shell:
            """
            pigz -d -c {params.fasta} > {output}
            """

    rule faidx_reffasta:
        input:
            "{tmpdir}/contig_correction/{sample}.fasta"
        output:
            temp("{tmpdir}/contig_correction/{sample}.fasta.fai")
        message: "Generate FastA index for de-novo assembled contigs: {wildcards.sample}"
        conda: "../envs/ENVS_samtools.yaml"
        resources:
            mem = 4,
            cores = 1
        shell:
            """
            samtools faidx {input}
            """

    rule determine_chunks_freebayes:
        input:
            fai = "{tmpdir}/contig_correction/{sample}.fasta.fai",
            depth = lambda wildcards: f"{wildcards.tmpdir}/alignment/{config['assembler']}/{wildcards.sample}.samtools_depth" 
        output:
            temp("{tmpdir}/contig_correction/{sample}.chunks")
        message: "Determine regions with approx. equal coverage to have 100 chuncks: {wildcards.sample}"
        conda: "../envs/ENVS_freebayes.yaml"
        resources:
            mem = 8,
            cores = 1
        threads: 1
        shell:
            """
            cat {input.depth} | \
            coverage_to_regions.py {input.fai} 100 > {output}
            """

    rule freebayes:
        input:
            chunks = "{tmpdir}/contig_correction/{sample}.chunks",
            fa = "{tmpdir}/contig_correction/{sample}.fasta",
            fai = "{tmpdir}/contig_correction/{sample}.fasta.fai",
            bam = lambda wildcards: f"{wildcards.tmpdir}/alignment/{config['assembler']}/{wildcards.sample}.sorted.noncorr.bam"
        output:
            pipe("{tmpdir}/contig_correction/{sample}.vcf")
        message: "Genotype the contigs using freeBayes in parallel mode: {wildcards.sample}"
        conda: "../envs/ENVS_freebayes.yaml"
        resources:
            mem = 32,
            cores = 16
        threads: 16
        shell:
            """
            freebayes-parallel {input.chunks} \
                    {threads} -f {input.fa} -F 0.33 -p 1 -q 20 {input.bam} > {output}
            """

    rule compress_vcf:
        input:
            "{tmpdir}/contig_correction/{sample}.vcf"
        output:
            "{tmpdir}/contig_correction/{sample}.vcf.gz"
        message: "Compress the VCF file produced by freebayes: {wildcards.sample}"
        conda: "../envs/ENVS_samtools.yaml"
        resources:
            mem = 4,
            cores = 1
        threads: 1
        shell:
            """
            bgzip -c {input} > {output}
            """

    rule bcftools_filter:
        input:
            lambda wildcards: f"{config['tmpdir']}/contig_correction/{wildcards.sample}.vcf.gz"
        output:
            vcf = "{resultdir}/consensus_correction/{assembler}/{sample}.filter.vcf.gz",
            tbi = "{resultdir}/consensus_correction/{assembler}/{sample}.filter.vcf.gz.tbi"
        message: "Discard low-quality differences between MEGAHIT and freebayes consensus: {wildcards.sample}"
        conda: "../envs/ENVS_bcftools.yaml"
        resources:
            mem = 4,
            cores = 1
        shell:
            """
            bcftools view \
                -v snps,mnps \
                -i 'QUAL >= 30 || (QUAL >= 20 && INFO/AO >= 3)' {input} | \
            bgzip > {output.vcf}
            bcftools index -t {output.vcf}
            """

    rule bcftools_consensus:
        input:
            fasta = lambda wildcards: f"{config['tmpdir']}/contig_correction/{wildcards.sample}.fasta",
            vcf = "{resultdir}/consensus_correction/{assembler}/{sample}.filter.vcf.gz",
            tbi = "{resultdir}/consensus_correction/{assembler}/{sample}.filter.vcf.gz.tbi"
        output:
            "{resultdir}/consensus_correction/{assembler}/{sample}_contigs.fasta.gz"
        message: "Correct the consensus sequence of the contigs: {wildcards.sample}"
        conda: "../envs/ENVS_bcftools.yaml"
        resources:
            mem = 8,
            cores = 2
        threads: 2
        shell:
            """
            cat {input.fasta} | bcftools consensus {input.vcf} | bgzip > {output}
            """

elif config['assembler'] == "metaspades":

    rule link_reffasta:
        input:
            lambda wildcards: f"{config['tmpdir']}/alignment/{wildcards.assembler}/{wildcards.sample}.raw.fasta"
        output:
            "{resultdir}/consensus_correction/{assembler}/{sample}_contigs.fasta.gz"
        message: "Link the FastA file with the de-novo assembled contigs of metaSPAdes without corrections: {wildcards.sample}"
        conda: "../envs/ENVS_samtools.yaml"
        resources:
            mem = 2,
            cores = 1
        threads: 1
        shell:
            """
            bgzip -c {input} > {output}
            """
