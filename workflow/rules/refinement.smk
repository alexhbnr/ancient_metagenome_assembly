if config['magrefinement']:

    if config['magrefiner'] == "metawrap":
    
        checkpoint metaWRAP_refinement:
            input:
                checkm = f"{config['resourcedir']}/checkM/setRoot.done",
                binning = "{resultdir}/binning/{sample}-{assembler}_binning.done"
            output:
                touch("{resultdir}/binning/{sample}-{assembler}_refinement.done")
            message: "Running the metaWRAP bin refinement module on {wildcards.sample}"
            resources:
                mem = 80,
                cores = 16
            params:
                outdir = "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}",
                min_completeness = config['min_completeness'],
                max_contamination = config['max_contamination'],
                maxbin2 = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/maxbin2_bins",
                metabat2 = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/metabat2_bins",
                concoct = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/concoct_bins"
            threads: 16
            wrapper:
                "https://www.github.com/alexhbnr/snakemake-wrappers/raw/main/bio/metawrap/bin_refinement"
