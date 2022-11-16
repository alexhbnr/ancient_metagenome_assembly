from glob import glob
import os
from pathlib import Path
import re

import pandas as pd

if len(snakemake.input.checkm) > 0:
    # CheckM
    checkm = pd.read_csv(snakemake.input.checkm, sep="\t")

    # GUNC
    gunc = pd.read_csv(snakemake.input.gunc, sep="\t")


    # BUSCO
    def parse_busco(fn):
        with open(fn, "rt") as buscofile:
            busco_counts = []
            for i, line in enumerate(buscofile):
                if i == 1:
                    lineage = re.search(r'# The lineage dataset is: ([a-z0-9_]+) \(Creation date: .+\)',
                                        line.rstrip()).group(1)
                elif i > 9 and i < 15:
                    busco_counts.append(tuple(line.rstrip().split("\t")[1:3]))
            return pd.DataFrame(busco_counts, columns=['nBUSCOs', 'typeBUSCOs']) \
                .assign(lineage=lineage,
                        binID=f"bin.{int(os.path.basename(os.path.dirname(fn)))}")

    busco = pd.concat([parse_busco(fn)
                       for b in snakemake.input.busco
                       for fn in glob(f"{b}/short_summary.*")])
    busco['nBUSCOs'] = busco['nBUSCOs'].astype(int)
    busco = pd.pivot_table(busco, values="nBUSCOs", index=['binID', 'lineage'],
                           columns=['typeBUSCOs'])
    busco.columns = ['D', 'S', 'F', 'M', 'T']
    busco = busco.reset_index()

    # Merge
    ## checkM & GUNC following gunc checkm_merge
    gunc = gunc[['genome', 'n_contigs', 'n_genes_called', 'n_genes_mapped',
                 'taxonomic_level', 'contamination_portion',
                 'n_effective_surplus_clades', 'reference_representation_score',
                 'clade_separation_score', 'pass.GUNC']] \
        .rename({'genome': 'binID',
                 'taxonomic_level': 'GUNC.divergence_level',
                 'reference_representation_score': 'GUNC.BRS',
                 'clade_separation_score': 'GUNC.CSS'}, axis=1)
    gunc.columns = [f"GUNC.{n}" if i in [1, 2, 3, 5, 6] else n
                    for i, n in enumerate(gunc.columns)]

    checkm = checkm[['Bin Id', 'Marker lineage', 'Genome size (bp)', 'GC',
                     'Coding density', 'N50 (contigs)', 'Longest contig (bp)',
                     'Mean contig length (bp)', 'Translation table',
                     'Completeness', 'Contamination', 'Strain heterogeneity']]
    checkm.columns = ['binID', 'checkM.lineage', 'checkM.genome_size',
                      'checkM.GC', 'checkM.coding_density', 'checkM.N50',
                      'checkM.longest_contig', 'checkM.mean_contig_length',
                      'checkM.translation_table', 'checkM.completeness',
                      'checkM.contamination', 'checkM.strain_heterogeneity']
    ## BUSCO
    busco_generic = busco.loc[busco['lineage'].isin(['bacteria_odb10', 'archaea_odb10'])] \
        .drop(['lineage'], axis=1)
    busco_generic.columns = [f"BUSCO_generic.{v}" if i > 0 else v
                             for i, v in enumerate(busco_generic.columns)]
    busco_specific = busco.loc[~busco['lineage'].isin(['bacteria_odb10', 'archaea_odb10'])]
    busco_specific.columns = [f"BUSCO_specific.{v}" if i > 0 else v
                              for i, v in enumerate(busco_specific.columns)]

    bin_summary = gunc.merge(checkm, how="left", on="binID") \
        .merge(busco_generic, how="left", on="binID") \
        .merge(busco_specific, how="left", on="binID") \
        .assign(sample=snakemake.wildcards.sample,
                assembler=snakemake.wildcards.assembler)
    bin_summary['pass.MIMAG_high'] = (bin_summary['checkM.completeness'] >= 90) & (bin_summary['checkM.contamination'] < 5)
    bin_summary['pass.MIMAG_medium'] = (bin_summary['checkM.completeness'] >= 50) & (bin_summary['checkM.contamination'] < 10)
    bin_summary.iloc[:, [32, 33, 0, 34, 35, 9] + list(range(1, 9)) + list(range(10, 32))] \
        .to_csv(snakemake.output[0], sep="\t", index=False)

else:
    Path(snakemake.output[0]).touch()
