# Ancient metagenome assembly pipeline

## Introduction

This is basic pipeline to *de-novo* assemble metagenomic sequencing data. It is written in
[Snakemake](https://snakemake.readthedocs.io/) and it was adapted to include ancient DNA specific
steps, such as the evaluation of the presence of ancient DNA damage using
[PyDamage](https://github.com/maxibor/pydamage) that are usually absent in metagenome assembly
pipelines. This pipeline heavily borrows ideas from the [Nextflow](https://www.nextflow.io/)
pipeline [MADMAN](https://github.com/maxibor/madman), which was deprecated and its ideas further
developed in the Nextflow pipeline [nf-core/mag](https://github.com/nf-core/mag).

## Quick start

To be able to run the pipeline, [Snakemake](https://snakemake.readthedocs.io/) with a minimal
version of 6.0 is necessary. The easiest way to install the dependencies of this program and to have
reproducible results is to create a new [conda](https://docs.conda.io/en/latest/) environment using
the environment file provided with it.

```
wget https://raw.githubusercontent.com/alexhbnr/ancient_metagenome_assembly/main/environment.yml
conda env create -f environment.yaml
```

After activating the environment using `conda activate ancient_metagenome_assembly`, the necessary
Python environment has been created.
