import sys

import pandas as pd

# Read sample table
sampletsv = pd.read_csv(snakemake.params.sampletsv, sep="\t", index_col=[0])
sampletsv.index = sampletsv.index.astype(str)

# Test for uniqueness of samples
if sampletsv.shape[0] > sampletsv.index.unique().shape[0]:
    print(f"Your sample TSV has {sampletsv.shape[0]} entries but only "
          f"{sampletsv.index.unique().shape[0]} unique sample identifiers. "
          "It seems that some sample identifiers are duplicate. Only unique "
          "sample identifiers are allow.")
    sys.exit(1)
# Test for paired-end data when assembler metaSPAdes was selected
elif snakemake.params.assembler == "metaspades" and sampletsv['R2'].isnull().sum() > 0:
    print("The selected assembler metaSPAdes requires paired-end sequencing "
          f"data for assembly but {sampletsv['R2'].isnull().sum() > 0} of your "
          "samples do not seem to have a second sequencing data file specified. "
          "Specify always both sequencing data files when selecting metaSPAdes "
          "for assembly.")
    sys.exit(1)
# Test for paired-end data when read correction via metaSPAdes was selected
elif snakemake.params.readcorrection == "true" and sampletsv['R2'].isnull().sum() > 0:
    print("The selected correction of reads uses SPAdes-hammer from the "
          "assembler metaSPAdes, which requires paired-end sequencing data for "
          f" assembly. However, {sampletsv['R2'].isnull().sum() > 0} of your "
          "samples do not seem to have a second sequencing data file specified. "
          "Specify always both sequencing data files when selecting read "
          "correction prior to assembly.")
    sys.exit(1)
