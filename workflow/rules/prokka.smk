rule prokka:
    input:
        "{resultdir}/alignment/{assembler}/{sample}-{assembler}.fasta.gz"
    output:
        "{resultdir}/prokka/{sample}-{assembler}.gff.gz"
    message: "Run Prokka on contigs: {wildcards.sample}"
    conda: "../envs/ENVS_prokka.yaml"
    resources:
        mem = 16,
        cores = 8
    params:
        tmpdir = lambda wildcards: f"{config['tmpdir']}/prokka_{wildcards.sample}_{wildcards.assembler}",
        outdir = "{resultdir}/prokka",
    threads: 8
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input} > {params.tmpdir}/{wildcards.sample}-{wildcards.assembler}.fasta
        prokka --outdir {params.tmpdir} \
               --prefix {wildcards.sample}-{wildcards.assembler} \
               --force \
               --compliant \
               --metagenome \
               --cpus {threads} \
               --debug \
               {params.tmpdir}/{wildcards.sample}-{wildcards.assembler}.fasta
        cp -r {params.tmpdir}/{wildcards.sample}-{wildcards.assembler}.{{faa,ffn,fna,gbk,gff,tsv,txt}} {params.outdir}/
        pigz -f -p 4 {params.outdir}/{wildcards.sample}-{wildcards.assembler}.{{faa,ffn,fna,gbk,gff,tsv,txt}}
        rm -r {params.tmpdir}/
        """
