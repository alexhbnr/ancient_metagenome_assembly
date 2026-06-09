if config['magrefinement']:

    rule refinement_workflow:
        input: 
            lambda wildcards: expand("{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins.stats", resultdir=[config['resultdir']], assembler=[config['assembler']], sample=successful_samples(wildcards))
        output:
            touch(f"{config['tmpdir']}/refinement.done")

    if config['magrefiner'] == "metawrap":
    
        checkpoint metaWRAP_refinement:
            input:
                checkm = f"{config['resourcedir']}/checkM/setRoot.done",
                binning = "{resultdir}/binning/{sample}-{assembler}_binning.done"
            output:
                # TODO: incorporate completion and contamination parameters into filename
                "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins.stats"
            message: "Running the metaWRAP bin refinement module on {wildcards.sample}"
            container: "docker://alexhbnr/metawrap-refinement:1.3.0--hdfd78af_3"
            resources:
                mem_mb = 80000,
                runtime = 2880,
                slurm_partition = "standard",
                cores = 16
            params:
                outdir = "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}",
                min_completeness = config['min_completeness'],
                max_contamination = config['max_contamination'],
                maxbin2 = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/maxbin2_bins",
                metabat2 = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/metabat2_bins",
                concoct = "{resultdir}/binning/metawrap/INITIAL_BINNING/{sample}-{assembler}/concoct_bins",
            threads: 16
            shell:
                """
                mkdir -p {params.outdir}
                metawrap bin_refinement -o {params.outdir} \
                    -t {threads} \
                    -c {params.min_completeness} \
                    -x {params.max_contamination} \
                    -A {params.maxbin2} \
                    -B {params.metabat2} \
                    -C {params.concoct} \
                    --skip-plotting || \
                touch {output}
                if [[ ! -d "{params.outdir}/metawrap_50_10_bins" ]]; then
                    mkdir {params.outdir}/metawrap_50_10_bins
                fi
                """

else:

    localrules: metaWRAP_refinement

    checkpoint metaWRAP_refinement:
        input:
            binning = "{resultdir}/binning/{sample}-{assembler}_binning.done"
        output:
            "{resultdir}/binning/metawrap/BIN_REFINEMENT/{sample}-{assembler}/metawrap_50_10_bins.stats"
        message: "Create dummy output of metaWRAP: {wildcards.sample}"
        shell:
            "touch {output}"
