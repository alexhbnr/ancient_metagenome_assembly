def expected_checkm_samples(wildcards):
    if "checkm" in config['quality_evaluation']:
        dn = f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins"
        if len(glob(f"{dn}/*.fa")) > 0:
            return f"{config['resultdir']}/stats/checkM/{wildcards.sample}-{wildcards.assembler}.checkM.txt"
        else:
            return []
    else:
        return []


def expected_busco_samples(wildcards):
    if "busco" in config['quality_evaluation']:
        dn = f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins"
        if len(glob(f"{dn}/*.fa")) > 0:
            return [f"{config['resultdir']}/stats/busco/{wildcards.sample}-{wildcards.assembler}/{int(os.path.basename(b).split('.')[1]):03d}"
                    for b in glob(f"{dn}/*.fa")]
        else:
            return []
    else:
        return []


def expected_gunc_samples(wildcards):
    if "gunc" in config['quality_evaluation']:
        dn = f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins"
        if len(glob(f"{dn}/*.fa")) > 0:
            return f"{config['resultdir']}/stats/gunc/{wildcards.sample}-{wildcards.assembler}/GUNC.progenomes_2.1.maxCSS_level.tsv"
        else:
            return []
    else:
        return []

wildcard_constraints:
    b = "[0-9]+"

rule quality_evaluation_workflow:
    input: 
        lambda wildcards: expand("{resultdir}/stats/bin_quality/{sample}-{assembler}_binqualities.tsv", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards))
    output:
        touch(f"{config['tmpdir']}/quality_evaluation.done")

rule quality_evaluation:
    input:
        checkm = lambda wildcards: expected_checkm_samples(wildcards),
        gunc = lambda wildcards: expected_gunc_samples(wildcards),
        busco = lambda wildcards: expected_busco_samples(wildcards)
    output:
        "{resultdir}/stats/bin_quality/{sample}-{assembler}_binqualities.tsv"
    message: "Summarise the bin qualities: {wildcards.sample}"
    script:
        "../scripts/summarise_checkM_gunc_busco.py"

if "checkm" in config['quality_evaluation'] or config['magrefinement']:

    localrules: checkM_prepareDatabase, checkM_setRoot

    rule checkM_prepareDatabase:
        output:
            f"{config['resourcedir']}/checkM/.dmanifest"
        message: "Download the checkM databases and extract the tarballs"
        container: "docker://quay.io/biocontainers/checkm-genome:1.2.5--pyhdfd78af_0"
        params:
            outdir = f"{config['resourcedir']}/checkM",
            url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
        shell:
            """
            cd {params.outdir} && \
            wget {params.url} && \
            tar xvf $(basename {params.url})
            """

    rule checkM_setRoot:
        input:
            f"{config['resourcedir']}/checkM/.dmanifest"
        output:
            touch(f"{config['resourcedir']}/checkM/setRoot.done")
        message: "Specify the database folder in checkM"
        container: "docker://quay.io/biocontainers/checkm-genome:1.2.5--pyhdfd78af_0"
        params:
            dbdir = f"{config['resourcedir']}/checkM"
        shell:
            """
            echo {params.dbdir} |
            checkm data setRoot {params.dbdir}
            """

if "checkm" in config['quality_evaluation']:

    rule checkM_lineage_wf:
        input:
            binning = lambda wildcards: f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins.stats",
            db = lambda wildcards: f"{config['resourcedir']}/checkM/setRoot.done" 
        output:
            "{tmpdir}/checkM/{sample}-{assembler}/storage/bin_stats.analyze.tsv"
        message: "Run checkM using lineage-specific workflow on sample {wildcards.sample}"
        container: "docker://quay.io/biocontainers/checkm-genome:1.2.5--pyhdfd78af_0"
        resources:
            mem_mb = 80000,
            runtime = 1440,
            slurm_partition = "standard",
            cores = 8
        params:
            fadir = lambda wildcards: f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins",
            outputfd = "{tmpdir}/checkM/{sample}-{assembler}"
        log: "{tmpdir}/checkM/{sample}-{assembler}.checkM.log"
        threads: 8
        shell:
            """
            checkm lineage_wf \
                -x fa \
                -t {threads} \
                {params.fadir} {params.outputfd} > {log}
            """

    rule checkM_qa:
        input:
            lambda wildcards: f"{config['tmpdir']}/checkM/{wildcards.sample}-{wildcards.assembler}/storage/bin_stats.analyze.tsv"
        output:
            "{resultdir}/stats/checkM/{sample}-{assembler}.checkM.txt"
        message: "Generate extended checkM report for sample {wildcards.sample}"
        container: "docker://quay.io/biocontainers/checkm-genome:1.2.5--pyhdfd78af_0"
        resources:
            mem_mb = 20000,
            runtime = 1440,
            slurm_partition = "standard",
            cores = 1
        params:
            outputfd = lambda wildcards: f"{config['tmpdir']}/checkM/{wildcards.sample}-{wildcards.assembler}"
        threads: 1
        shell:
            """
            checkm qa -o 2 --tab_table -f {output} \
                {params.outputfd}/lineage.ms {params.outputfd}
            """


if "busco" in config['quality_evaluation']:

    rule busco:
        input:
            "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins.stats"
        output:
            directory("{resultdir}/stats/busco/{sample}-{assembler}/{b}")
        message: "Run BUSCO on the contigs: {wildcards.sample}"
        container: "docker://quay.io/biocontainers/busco:6.0.0--pyhdfd78af_2"
        resources:
            mem_mb = 36000,
            runtime = 1440,
            slurm_partition = "standard",
            cores = 8
        params:
            fasta =  lambda wildcards: f"{wildcards.resultdir}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins/bin.{int(wildcards.b)}.fa",
            outprefix = lambda wildcards: wildcards.b,
            outdir = "{resultdir}/stats/busco/{sample}-{assembler}",
            mode = "genome",
            lineage = "--auto-lineage-prok",
            downloads_path = f"{config['resourcedir']}/busco_downloads"
        threads: 8
        shell:
            """
            busco --in {params.fasta} --out {params.outprefix} \
                --force --out_path {params.outdir} \
                --cpu {threads} --mode {params.mode} \
                {params.lineage} --tar --download_path {params.downloads_path}
            """

if "gunc" in config['quality_evaluation']:

    rule install_gunc_database:
        output:
            directory("{resourcedir}/GUNC/db")
        message: "Install GUNC database"
        container: "docker://quay.io/biocontainers/gunc:1.1.0--pyhdfd78af_0"
        resources:
            mem_mb = 8000,
            runtime = 1440,
            slurm_partition = "standard",
        params:
            dir = "{resourcedir}/GUNC/db"
        shell:
            """
            gunc download_db {params.dir}
            """

    rule gunc_run:
        input:
            binning = "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins.stats",
            db = lambda wildcards: f"{config['resourcedir']}/GUNC/db"
        output:
            "{resultdir}/stats/gunc/{sample}-{assembler}/GUNC.progenomes_2.1.maxCSS_level.tsv"
        message: "Run GUNC on bins: {wildcards.sample}"
        container: "docker://quay.io/biocontainers/gunc:1.1.0--pyhdfd78af_0"
        resources:
            mem_mb = 32000,
            runtime = 1440,
            slurm_partition = "standard",
            cores = 8
        params:
            db_file = lambda wildcards: f"{config['resourcedir']}/GUNC/db/gunc_db_progenomes2.1.dmnd",
            fadir = "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins",
            outdir = "{resultdir}/stats/gunc/{sample}-{assembler}",
            tmpdir = lambda wildcards: f"{config['tmpdir']}/gunc/{wildcards.sample}-{wildcards.assembler}"
        threads: 8
        shell:
            """
            mkdir -p {params.tmpdir} && \
            gunc run -r {params.db_file} \
                --input_dir {params.fadir} \
                -t {threads} \
                -o {params.outdir} \
                --temp_dir {params.tmpdir} \
                --detailed_output
            """
