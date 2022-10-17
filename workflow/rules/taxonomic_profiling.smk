def expected_bins(wildcards):
    if config['taxonomic_profiling']:
        return [f"{config['tmpdir']}/tax_profiling/bins/{wildcards.sample}-{wildcards.assembler}_{os.path.basename(fafn)}"
                for s in SAMPLES
                for fafn in glob(f"{dn}/{s}-{config['assembler']}/metawrap_50_10_bins/*.fa")]
    else:
        return []


wildcard_constraints:
    b = "[0-9]+"

if config['taxonomic_profiling']:

    rule link_bin_fastas:
        input:
            [f"{config['resultdir']}/binning/{sample}-{config['assembler']}_refinement.done" for sample in SAMPLES]
        output:
            touch("{tmpdir}/tax_profiling/fas_linked")
        message: "Softlink the bins into a temporary directory"
        params:
            outdir = "{tmpdir}/tax_profiling/bins"
        run:
            os.makedirs(params.outdir, exist_ok=True)
            dn = f"{config['resultdir']}/binning/metawrap/BIN_REFINEMENT"
            if not dn.startswith("/"):
                dn = f"{os.getcwd()}/{dn}"
            for s in SAMPLES:
                for fafn in glob(f"{dn}/{s}-{config['assembler']}/metawrap_50_10_bins/*.fa"):
                    if not os.path.islink(f"{params.outdir}/{s}-{config['assembler']}_{os.path.basename(fafn)}"):
                        os.symlink(fafn, f"{params.outdir}/{s}-{config['assembler']}_{os.path.basename(fafn)}")
            Path(output[0]).touch()

    if "gtdbtk" in config['taxprofilers']:

        rule gtdbtk_download_db:
            output:
                "{resourcedir}/gtdbtk/gtdbtk_r207_v2/metadata/metadata.txt"
            message: "Download and set-up the GTDBTK database"
            params:
                url = "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz",
                resourcedir = config['resourcedir']
            wrapper:
                "file:///home/alexander_huebner/github/snakemake-wrappers/bio/gtdbtk/download_db"

        rule gtdbtk_classify:
            input:
                db = f"{config['resourcedir']}/gtdbtk/gtdbtk_r207_v2/metadata/metadata.txt",
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
                dbdir = f"{config['resourcedir']}/gtdbtk/gtdbtk_r207_v2"
            threads: 32
            wrapper:
                "file:///home/alexander_huebner/github/snakemake-wrappers/bio/gtdbtk/classify_wf"

    if "phylophlan3" in config['taxprofilers']:

        rule phylophlan3:
            input:
                fas = f"{config['tmpdir']}/tax_profiling/fas_linked"
            output:
                "{resultdir}/stats/phylophlan3/phylophlan3_mash.tsv"
            message: "Assign the MAGs to a taxonomy based on MASH distances"
            resources:
                mem = 40,
                cores = 32
            params:
                dbdir = f"{config['resourcedir']}/phylophlan3",
                db = "SGB.Dec20",
                fadir = f"{config['tmpdir']}/tax_profiling/bins",
                prefix = "{resultdir}/stats/phylophlan3/phylophlan"
            wrapper:
                "file:///home/alexander_huebner/github/snakemake-wrappers/bio/phylophlan/metagenomic"
