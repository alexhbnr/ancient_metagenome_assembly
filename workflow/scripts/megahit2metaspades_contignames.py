import bgzip
import pandas as pd
import pysam
import pyfastx

# Read the average depth per contig
depth = pd.read_csv(snakemake.input.depth, sep="\t", index_col=[0])
# Read FastA file and store in memory
seqs = {name: seq
        for name, seq in pyfastx.Fasta(snakemake.input.fasta, build_index=False)}
# Determine lengths
lengths = {name: len(seq)
           for name, seq in seqs.items()}
# Sort contig names by length descending
contigs = [name for name, length in sorted(lengths.items(), key=lambda kv: kv[1])[::-1]]
# Generate alternative names:
contigs_metaspades_names = {name: f"NODE_{i+1}_length_{lengths[name]}_cov_{depth.at[name, 'depth']:.1f}"
                            for i, name in enumerate(contigs)}

# Write FastA output
with open(snakemake.output.fasta, "wb") as raw:
    with bgzip.BGZipWriter(raw) as outfile:
        for name in contigs:
            outfile.write(str(f">{contigs_metaspades_names[name]}\n{seqs[name]}\n").encode("utf-8"))
# Fix BAM header
bamfile = pysam.AlignmentFile(snakemake.input.bam)
header = bamfile.header.to_dict()
header['SQ'] = [{'SN': contigs_metaspades_names[contig['SN']],
                 'LN': contig['LN']}
                for contig in header['SQ']]
with pysam.AlignmentFile(snakemake.params.header, "wh", header=header):
    pass
