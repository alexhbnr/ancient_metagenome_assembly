wildcard_constraints:
    b = "[0-9]+"

if config['taxonomic_profiling']:

    rule taxonomic_profiling_workflow:
        input: 
            lambda wildcards: expand("{resultdir}/stats/gtdbtk/gtdbtk.bac120.summary.tsv", resultdir=[config['resultdir']]),
            lambda wildcards: expand("{resultdir}/stats/phylophlan3/genome_assignment.tsv", resultdir=[config['resultdir']])
        output:
            touch(f"{config['tmpdir']}/taxonomic_profiling.done")

    localrules: link_bin_fastas

    rule link_bin_fastas:
        input:
            samples = "{tmpdir}/successful_samples.txt",
            refinement = "{tmpdir}/refinement.done"
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

        localrules: gtdbtk_download_db

        rule gtdbtk_download_db:
            output:
                f"{config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}/metadata/metadata.txt"
            message: "Download and set-up the GTDBTK database"
            container:
                lambda wildcards: {
                        'r207_v2': "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0",
                        'r214.1': "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0",
                        'r220.0': "https://depot.galaxyproject.org/singularity/gtdbtk%3A2.4.1--pyhdfd78af_1",
                        'r226.0': "https://depot.galaxyproject.org/singularity/gtdbtk%3A2.4.1--pyhdfd78af_1",
                        'r232.0': "docker://quay.io/biocontainers/gtdbtk:2.7.1--pyhdfd78af_1",
                        }[config['gtdb_version']]
            params:
                url = {'r207_v2': "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz",
                       'r214.1': "https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz",
                       'r220.0': "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz",
                       'r226.0': "https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz",
                       'r232.0': "https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz",
                }[config['gtdb_version']],
                resourcedir = config['resourcedir']
            shell:
                """
                wget -O {params.resourcedir}/gtdbtk/$(basename {params.url}) {params.url} && \
                tar -xvzf {params.resourcedir}/gtdbtk/$(basename {params.url}) \
                    -C '{params.resourcedir}/gtdbtk/$(basename {params.url} _data.tar.gz)' \
                    --strip 1 > /dev/null && \
                rm {params.resourcedir}/gtdbtk/$(basename {params.url})
                """

        rule gtdbtk_classify:
            input:
                db = f"{config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}/metadata/metadata.txt",
                fas = f"{config['tmpdir']}/tax_profiling/fas_linked"
            output:
                "{resultdir}/stats/gtdbtk/gtdbtk.bac120.summary.tsv"
            message: "Run the GTDBTK's classify workflow"
            container:
                lambda wildcards: {
                        'r207_v2': "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0",
                        'r214.1': "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0",
                        'r220.0': "https://depot.galaxyproject.org/singularity/gtdbtk%3A2.4.1--pyhdfd78af_1",
                        'r226.0': "https://depot.galaxyproject.org/singularity/gtdbtk%3A2.4.1--pyhdfd78af_1",
                        'r232.0': "docker://quay.io/biocontainers/gtdbtk:2.7.1--pyhdfd78af_1",
                        }[config['gtdb_version']]
            resources:
                mem_mb = 80000,
                runtime = 1440,
                slurm_partition = "standard",
                cores = 32
            params:
                fadir = f"{config['tmpdir']}/tax_profiling/bins",
                mash = lambda wildcards: f"--mash_db {config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}" if config['gtdb_version'] in ['r207_v2', 'r214.1', 'r220.0', 'r226.0'] else "",
                outdir = "{resultdir}/stats/gtdbtk",
                dbdir = f"{config['resourcedir']}/gtdbtk/gtdbtk_{config['gtdb_version']}"
            threads: 32
            shell:
                """
                GTDBTK_DATA_PATH={params.dbdir} \
                gtdbtk classify_wf --cpu {threads} \
                    {params.mash} \
                    --extension fa --genome_dir {params.fadir} --out_dir {params.outdir}
                """

    if "phylophlan3" in config['taxprofilers']:

        rule phylophlan3:
            input:
                fas = f"{config['tmpdir']}/tax_profiling/fas_linked"
            output:
                "{resultdir}/stats/phylophlan3/genome_assignment.tsv"
            message: "Assign the MAGs to a taxonomy based on MASH distances"
            conda: "../envs/ENVS_phylophlan3.yaml"
            resources:
                mem_mb = lambda wildcards, attempt: 40000 + attempt * 40000,
                runtime = 1440,
                slurm_partition = "standard",
                cores = 32
            params:
                dbdir = f"{config['resourcedir']}/phylophlan3",
                db = "SGB.Jan25",
                fadir = f"{config['tmpdir']}/tax_profiling/bins",
                outdir = "{resultdir}/stats/phylophlan3"
            threads: 32
            shell:
                """
                phylophlan_assign_sgbs \
                    -i {params.fadir} \
                    -e fa \
                    -o {params.outdir} \
                    --database_folder {params.dbdir} \
                    --database {params.db} \
                    --nproc_cpu {threads}
                """
