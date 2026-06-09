import os
import shutil

import pandas as pd
import pyfastx
from tqdm import tqdm

metawrapdir = snakemake.params.get("metawrapdir", "")

if os.stat(snakemake.input[0]).st_size > 0:
    clusters = pd.read_csv(snakemake.input[0], sep=",")
    contig_lengths = {name: len(seq)
                      for name, seq in pyfastx.Fasta(snakemake.input[1], build_index=False)}

    # Filter clusters that surpass the minimal length of 50 kb
    clusters['length'] = [contig_lengths[c] for c in clusters['contig_id']]
    total_length = clusters.groupby(['cluster_id'])['length'].agg(sum)
    final_clusters = total_length.loc[total_length > 50000].index.tolist()
    binned_contigs = clusters.loc[clusters['cluster_id'].isin(final_clusters)] \
        .set_index(['contig_id'])

    # Open files for writing
    os.makedirs(metawrapdir, exist_ok=True)
    outfiles = [open(f"{metawrapdir}/bin.{i}.fa", "wt") for i in final_clusters]
    unbinned_outfile = open(f"{metawrapdir}/unbinned.fa", "wt")
    outfile_map = {c: i for i, c in enumerate(final_clusters)}
    for name, seq in tqdm(pyfastx.Fasta(snakemake.input[1], build_index=False)):
        if name in binned_contigs.index:
            outfiles[outfile_map[binned_contigs.at[name, 'cluster_id']]].write(f">{name}\n{seq}\n")
        else:
            unbinned_outfile.write(f">{name}\n{seq}\n")
        if len(seq) < snakemake.params.min_contiglength:
            break

    # Close files for writing
    for f in outfiles:
        f.close()
    unbinned_outfile.close()
else:
    clustered_contigs = []

shutil.copy(snakemake.input[0],
            f"{metawrapdir}/../{snakemake.wildcards.sample}-{snakemake.wildcards.assembler}.concoct_clustering.csv")
