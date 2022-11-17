def expected_checkm_samples(wildcards):
    if "checkm" in config['quality_evaluation']:
        dn = f"{os.path.dirname(checkpoints.metaWRAP_refinement.get(**wildcards).output[0])}/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins"
        if len(glob(f"{dn}/*.fa")) > 0:
            return f"{config['resultdir']}/stats/checkM/{wildcards.sample}-{wildcards.assembler}.checkM.txt"
        else:
            return []
    else:
        return []


def expected_busco_samples(wildcards):
    if "busco" in config['quality_evaluation']:
        dn = f"{os.path.dirname(checkpoints.metaWRAP_refinement.get(**wildcards).output[0])}/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins"
        if len(glob(f"{dn}/*.fa")) > 0:
            return [f"{config['resultdir']}/stats/busco/{wildcards.sample}-{wildcards.assembler}/{int(b.split('.')[1]):03d}"
                    for b in glob(f"{dn}/*.fa")]
        else:
            return []
    else:
        return []


def expected_gunc_samples(wildcards):
    if "gunc" in config['quality_evaluation']:
        dn = f"{os.path.dirname(checkpoints.metaWRAP_refinement.get(**wildcards).output[0])}/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins"
        if len(glob(f"{dn}/*.fa")) > 0:
            return f"{config['resultdir']}/stats/gunc/{wildcards.sample}-{wildcards.assembler}/GUNC.progenomes_2.1.maxCSS_level.tsv"
        else:
            return []
    else:
        return []

wildcard_constraints:
    b = "[0-9]+"

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

if "checkm" in config['quality_evaluation']:

    rule checkM_prepareDatabase:
        output:
            f"{config['resourcedir']}/checkM/.dmanifest"
        message: "Download the checkM databases and extract the tarballs"
        params:
            outdir = f"{config['resourcedir']}/checkM",
            url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/checkm/preparedatabase"

    rule checkM_setRoot:
        input:
            f"{config['resourcedir']}/checkM/.dmanifest"
        output:
            touch(f"{config['resourcedir']}/checkM/setRoot.done")
        message: "Specify the database folder in checkM"
        params:
            dbdir = f"{config['resourcedir']}/checkM"
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/checkm/setroot"

    rule checkM_lineage_wf:
        input:
            binning = lambda wildcards: f"{config['resultdir']}/binning/{wildcards.sample}-{wildcards.assembler}_refinement.done",
            db = lambda wildcards: f"{config['resourcedir']}/checkM/setRoot.done" 
        output:
            "{tmpdir}/checkM/{sample}-{assembler}/storage/bin_stats.analyze.tsv"
        message: "Run checkM using lineage-specific workflow on sample {wildcards.sample}"
        resources:
            mem = 80,
            cores = 8
        params:
            fadir = lambda wildcards: f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins",
            outputfd = "{tmpdir}/checkM/{sample}-{assembler}"
        log: "{tmpdir}/checkM/{sample}-{assembler}.checkM.log"
        threads: 8
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/checkm/lineage_wf"

    rule checkM_qa:
        input:
            lambda wildcards: f"{config['tmpdir']}/checkM/{wildcards.sample}-{wildcards.assembler}/storage/bin_stats.analyze.tsv"
        output:
            "{resultdir}/stats/checkM/{sample}-{assembler}.checkM.txt"
        message: "Generate extended checkM report for sample {wildcards.sample}"
        resources:
            mem = 20,
            cores = 1
        params:
            outputfd = lambda wildcards: f"{config['tmpdir']}/checkM/{wildcards.sample}-{wildcards.assembler}"
        threads: 1
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/checkm/qa"


if "busco" in config['quality_evaluation']:

    rule busco:
        input:
            lambda wildcards: f"{wildcards.resultdir}/binning/metawrap/BIN_REFINEMENT/{wildcards.sample}-{wildcards.assembler}/metawrap_50_10_bins/bin.{int(wildcards.b)}.fa"
        output:
            directory("{resultdir}/stats/busco/{sample}-{assembler}/{b}")
        message: "Run BUSCO on the contigs: {wildcards.sample}"
        resources:
            mem = 36,
            cores = 8
        params:
            mode = "genome",
            lineage = "auto-lineage-prok",
            downloads_path = f"{config['resourcedir']}/busco_downloads"
        threads: 8
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers-public/raw/master/bio/busco"

if "gunc" in config['quality_evaluation']:

    rule install_gunc_database:
        output:
            directory("{resourcedir}/GUNC/db")
        message: "Install GUNC database"
        params:
            dir = "{resourcedir}/GUNC/db"
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/gunc/download_db"

    rule gunc_run:
        input:
            binning = "{resultdir}/binning/{sample}-{assembler}_refinement.done",
            db = lambda wildcards: f"{config['resourcedir']}/GUNC/db"
        output:
            "{resultdir}/stats/gunc/{sample}-{assembler}/GUNC.progenomes_2.1.maxCSS_level.tsv"
        message: "Run GUNC on bins: {wildcards.sample}"
        resources:
            mem = 32,
            cores = 8
        params:
            db_file = lambda wildcards: f"{config['resourcedir']}/GUNC/db/gunc_db_progenomes2.1.dmnd",
            fadir = "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins",
            outdir = "{resultdir}/stats/gunc/{sample}-{assembler}",
            tmpdir = lambda wildcards: f"{config['tmpdir']}/gunc/{wildcards.sample}-{wildcards.assembler}"
        threads: 8
        wrapper:
            "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/gunc/run"
