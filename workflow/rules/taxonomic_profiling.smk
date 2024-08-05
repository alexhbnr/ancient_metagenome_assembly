def expected_bins(wildcards):
    if config['taxonomic_profiling']:
        return [f"{config['tmpdir']}/tax_profiling/bins/{wildcards.sample}-{wildcards.assembler}_{os.path.basename(fafn)}"
                for s in successful_samples(wildcards)
                for fafn in glob(f"{dn}/{s}-{config['assembler']}/metawrap_50_10_bins/*.fa")]
    else:
        return []


wildcard_constraints:
    b = "[0-9]+"

if config['taxonomic_profiling']:

    rule link_bin_fastas:
        input:
            lambda wildcards: [f"{config['resultdir']}/binning/{sample}-{config['assembler']}_refinement.done" for sample in successful_samples(wildcards)]
        output:
            touch("{tmpdir}/tax_profiling/fas_linked")
        message: "Softlink the bins into a temporary directory"
        params:
            outdir = "{tmpdir}/tax_profiling/bins",
            samples = lambda wildcards: successful_samples(wildcards)
        run:
            os.makedirs(params.outdir, exist_ok=True)
            dn = f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT"
            if not dn.startswith("/"):
                dn = f"{os.getcwd()}/{dn}"
            for s in params.samples:
                for fafn in glob(f"{dn}/{s}-{config['assembler']}/metawrap_50_10_bins/*.fa"):
                    if not os.path.islink(f"{params.outdir}/{s}-{config['assembler']}_{os.path.basename(fafn)}"):
                        os.symlink(fafn, f"{params.outdir}/{s}-{config['assembler']}_{os.path.basename(fafn)}")
            Path(output[0]).touch()

    if "gtdbtk" in config['taxprofilers']:

        rule gtdbtk_download_db:
            output:
                f"{config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}/metadata/metadata.txt"
            message: "Download and set-up the GTDBTK database"
            container: "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0"
            params:
                url = {'r207_v2': "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz",
                       'r214.1': "https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz",
                       'r220.0': "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz",
                       }[config['gtdb_version']],
                resourcesdir = config['resourcedir']
            wrapper:
                "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/gtdbtk/download_db"

        rule gtdbtk_classify:
            input:
                db = f"{config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}/metadata/metadata.txt",
                fas = f"{config['tmpdir']}/tax_profiling/fas_linked"
            output:
                "{resultdir}/stats/gtdbtk/gtdbtk.bac120.summary.tsv"
            message: "Run the GTDBTK's classify workflow"
            resources:
                mem = 80,
                cores = 32
            params:
                fadir = f"{config['tmpdir']}/tax_profiling/bins",
                outdir = "{resultdir}/stats/gtdbtk",
                dbdir = f"{config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}"
            threads: 32
            wrapper:
                "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/gtdbtk/classify_wf"

    if "phylophlan3" in config['taxprofilers']:

        rule phylophlan3:
            input:
                fas = f"{config['tmpdir']}/tax_profiling/fas_linked"
            output:
                "{resultdir}/stats/phylophlan3/phylophlan3_mash.tsv"
            message: "Assign the MAGs to a taxonomy based on MASH distances"
            resources:
                mem = lambda wildcards, attempt: 40 + attempt * 40,
                cores = 32
            params:
                dbdir = f"{config['resourcedir']}/phylophlan3",
                db = "SGB.Dec20",
                fadir = f"{config['tmpdir']}/tax_profiling/bins",
                prefix = "{resultdir}/stats/phylophlan3/phylophlan"
            threads: 32
            wrapper:
                "https://github.com/alexhbnr/snakemake-wrappers/raw/main/bio/phylophlan/metagenomic"
