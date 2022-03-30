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
conda env create -f environment.yml
```

After activating the environment using

```
conda activate ancient_metagenome_assembly
```

the necessary Python environment has been created.

Next, the pipeline itself needs to be downloaded. This can be easily down by cloning this repository
to your computer via git

```
git clone https://github.com/alexhbnr/ancient_metagenome_assembly.git
```

or by downloading the zip file and extracting this

```
wget -O ancient_metagenome_assembly.zip https://github.com/alexhbnr/ancient_metagenome_assembly/archive/refs/heads/main.zip
unzip ancient_metagenome_assembly.zip
```

Finally, the configuration file and sample table have to be provided. Templates for these can be
found in `config/config.yaml` for the configuration file and in `test/samples.tsv` for the sample
table.

To run a test case using the sequencing data of the Chimpanzee dental calculus sample CDC010
published by James Fellow Yates *et al.* (PNAS, 2021;
[doi:10.1073/pnas.2021655118](https://doi.org/10.1073/pnas.2021655118)), we will need to download
the FastQ files from ENA:

```
wget -O test/ERR3579712_1.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR357/002/ERR3579712/ERR3579712_1.fastq.gz
wget -O test/ERR3579712_2.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR357/002/ERR3579712/ERR3579712_2.fastq.gz
```

To start the pipeline, we run

```
snakemake --use-conda --conda-prefix conda -j 8
```

This will automatically evaluate the entries in the configuration file `config/config.yaml` and use
the sample `CDC10` as input. The temporary files are written into the folder `tmp` and the results
in the folder `results`.

By activating the options `--use-conda`, `snakemake` will download and install the programs
necessary to run the pipeline via conda and store them into the folder `conda`. This step will only
happen once, at the first time or when you change the folder for storing these conda environments
using `--conda-prefix`.
