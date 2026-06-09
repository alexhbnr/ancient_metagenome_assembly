if config['contig_annotation']:

    rule genome_annotation:
        input:
            lambda wildcards: expand("{resultdir}/bakta/{sample}-{assembler}.gff3.gz", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards))
        output:
            touch(f"{config['tmpdir']}/genome_annotation.done")

    rule bakta:
        input:
            "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
        output:
            "{resultdir}/bakta/{sample}-{assembler}.gff3.gz"
        message: "Run Bakta on contigs: {wildcards.sample}"
        # container: "docker://quay.io/biocontainers/bakta:1.12.0--pyhdfd78af_0"
        conda: "../envs/ENVS_bakta_v1.12.0.yaml"
        resources:
            mem_mb = 128000,
            runtime = 20160,
            slurm_partition = "standard",
            cores = 8
        params:
            dbdir = f"{config['resourcedir']}/bakta",
            tmpdir = lambda wildcards: f"{config['tmpdir']}/bakta_{wildcards.sample}_{wildcards.assembler}",
            extra = "--keep-contig-headers --meta --skip-crispr",
            outdir = "{resultdir}/bakta"
        threads: 8
        shell:
            """
            mkdir -p {params.tmpdir}
            bakta -p {wildcards.sample}-{wildcards.assembler} \
                --db {params.dbdir} \
                --force \
                --output {params.tmpdir} \
                {params.extra} \
                --threads {threads} \
                {input}
            cp -r {params.tmpdir}/{wildcards.sample}-{wildcards.assembler}.{{embl,faa,ffn,fna,gff3,hypotheticals.faa,hypotheticals.tsv,inference.tsv,json,tsv,txt}} {params.outdir}/
            pigz -f -p 4 {params.outdir}/{wildcards.sample}-{wildcards.assembler}.{{embl,faa,ffn,fna,gff3,hypotheticals.faa,hypotheticals.tsv,inference.tsv,json,tsv}}
            rm -r {params.tmpdir}/
            """

else:

    localrules: bakta_dummy_output

    rule bakta_dummy_output:
        output:
            "{resultdir}/bakta/{sample}-{assembler}.gff3.gz"
        message: "Create empty dummy files for Prokka analysis: {wildcards.sample}"
        shell:
            """
            touch {output}
            """
