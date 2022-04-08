if config['magrefinement']:

    if config['magrefiner'] == "metawrap":
    
        rule metaWRAP_refinement:
            input:
                checkm = f"{config['resourcedir']}/checkM/setRoot.done",
                binning = "{resultdir}/magbinning/{sample}-{assembler}_binning.done"
            output:
                touch("{resultdir}/magbinning/{sample}-{assembler}_refinement.done")
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
                "file:///home/alexander_huebner/github/snakemake-wrappers/bio/metawrap/bin_refinement"
