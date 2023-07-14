import gzip
import os
import re

import pandas as pd

# CalN50
total_length = []
number_contigs = []
nx_lx = []
## Process individual files
for s in snakemake.params.samples:
    fn = f"{snakemake.params.caln50_dir}/{s}-{snakemake.params.assembler}.caln"
    nl = []
    with open(fn, "rt") as caln50file:
        for line in caln50file:
            infotype = line.split("\t")[0]
            if infotype == "SZ":
                total_length.append((s, line.rstrip().split("\t")[1]))
            elif infotype == "NN":
                number_contigs.append((s, line.rstrip().split("\t")[1]))
            elif infotype == "NL":
                nl.append(line.rstrip().split("\t")[1:])
            nx_lx.append(pd.DataFrame(nl, columns=['x', 'Nx', 'Lx']) \
                .assign(sample=s))
## Summarise total length
total_length_df = pd.DataFrame(total_length, columns=["sample", "length"])
total_length_df['length'] = total_length_df['length'].astype(int)
## Summarise number of contigs
n_contigs_df = pd.DataFrame(number_contigs, columns=["sample", "nContigs"])
n_contigs_df['nContigs'] = n_contigs_df['nContigs'].astype(int)
## Summarise N50
nx_lx_df = pd.concat(nx_lx)
nx_lx_df['Nx'] = nx_lx_df['Nx'].astype(int)
nx_lx_df['Lx'] = nx_lx_df['Lx'].astype(int)
nx_lx_df['x'] = nx_lx_df['x'].astype(int)
nx_lx_df = nx_lx_df.drop_duplicates() \
    .sort_values(['sample']) \
    [['sample', 'x', 'Nx', 'Lx']]

nx_df = pd.pivot_table(nx_lx_df[['sample', 'x', 'Nx']],
                       index="sample", columns=["x"])
nx_df.columns = [f"N{i}" for i in range(0, 101, 10)]
nx_df = nx_df.reset_index()

# Flagstat
reads = []
for s in snakemake.params.samples:
    for i, line in enumerate(open(f"{snakemake.params.flagstat_dir}/{s}-{snakemake.params.assembler}.flagstat", "rt")):
        if i == 1:
            total_reads = int(line.split(" ")[0])
        elif i == 7:
            mapped_reads = int(line.split(" ")[0])
    reads.append((s, total_reads, mapped_reads)) 
reads_df = pd.DataFrame(reads, columns=['sample', 'totalReads', 'mappedReads'])
reads_df['fracMapped'] = reads_df['mappedReads'] / reads_df['totalReads']

# MetaQUAST
metaquast = []
for s in snakemake.params.samples:
    if os.path.isfile(f"{snakemake.params.metaquast_dir}/{s}-{snakemake.params.assembler}/combined_reference/transposed_report.tsv"):
        fn = f"{snakemake.params.metaquast_dir}/{s}-{snakemake.params.assembler}/combined_reference/transposed_report.tsv"
    elif os.path.isfile(f"{snakemake.params.metaquast_dir}/{s}-{snakemake.params.assembler}/transposed_report.tsv"):
        fn = f"{snakemake.params.metaquast_dir}/{s}-{snakemake.params.assembler}/transposed_report.tsv"
    else:
        fn = ""
    if fn != "":
        report = pd.read_csv(fn,
                            sep="\t", usecols=['# contigs (>= 1000 bp)', '# contigs (>= 5000 bp)', '# contigs (>= 10000 bp)',
                                                '# contigs (>= 25000 bp)', '# contigs (>= 50000 bp)']) \
            .assign(sample=s)
        metaquast.append(report)
metaquast_df = pd.concat(metaquast)
metaquast_df.columns = [f"nContigs_{size}bp"
                        for c in metaquast_df.columns[:-1]
                        for size in re.search(r"# contigs \(>= ([0-9]+) bp\)", c).groups(1)] + ['sample']

# Prokka
prokka = []
for s in snakemake.params.samples:
    if os.stat(f"{snakemake.params.prokka_dir}/{s}-{snakemake.params.assembler}.gff.gz").st_size > 50:
        annotations = dict(line.rstrip().split(": ")
                        for i, line in enumerate(gzip.open(f"{snakemake.params.prokka_dir}/{s}-{snakemake.params.assembler}.txt.gz", "rt"))
                        if i > 0)
        for k in ['CDS', 'gene', 'rRNA', 'tRNA', 'tmRNA']:
            annotations.setdefault(k, 0)
        prokka.append((s, int(annotations['CDS']), int(annotations['gene']),
                    int(annotations['rRNA']), int(annotations['tRNA']),
                    int(annotations['tmRNA'])))
    else:
        prokka.append((s, "NA", "NA", "NA", "NA", "NA"))
prokka_df = pd.DataFrame(prokka, columns=['sample', 'CDS', 'genes', 'rRNA', 'tRNA', 'tmRNA'])

# PyDamage
if snakemake.params.pydamage:
    pydamage = []
    for s in snakemake.params.samples:
        pyd_res = pd.read_csv(f"{snakemake.params.pydamage_dir}/{s}-{snakemake.params.assembler}.pydamage.csv.gz",
                            sep=",")
        pydamage.append(pyd_res[['predicted_accuracy', 'qvalue', 'nb_reads_aligned', 'coverage',
                                'reflen', 'CtoT-0', 'CtoT-1', 'CtoT-2']]
                        .median()
                        .to_frame()
                        .transpose()
                        .assign(sample=s))
    pydamage_df = pd.concat(pydamage)
    pydamage_df.columns = ['PyDamage_predAccuracy', 'PyDamage_Qvalue', 'PyDamage_nReads_contig',
                        'PyDamage_coverage', 'PyDamage_contiglength', 'PyDamage_CtoT-1',
                        'PyDamage_CtoT-2', 'PyDamage_CtoT-3', 'sample']

    # Combine
    df = total_length_df \
        .merge(n_contigs_df, how="left", on="sample") \
        .merge(metaquast_df, how="left", on="sample") \
        .merge(nx_df, how="left", on="sample") \
        .merge(prokka_df, how="left", on="sample") \
        .merge(pydamage_df, how="left", on="sample") \
        .fillna(0)
    df['length'] = df['length'].astype(int)
    df.iloc[:,3:8] = df.iloc[:,3:8].astype(int)
    df['PyDamage_nReads_contig'] = df['PyDamage_nReads_contig'].astype(int)
    df['PyDamage_contiglength'] = df['PyDamage_contiglength'].astype(int)
else:
    df = total_length_df \
        .merge(n_contigs_df, how="left", on="sample") \
        .merge(metaquast_df, how="left", on="sample") \
        .merge(nx_df, how="left", on="sample") \
        .merge(prokka_df, how="left", on="sample") \
        .fillna(0)
    df['length'] = df['length'].astype(int)
    df.iloc[:,3:8] = df.iloc[:,3:8].astype(int)

df.to_csv(snakemake.output[0], sep="\t", index=False)
