import h5py
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import functools
from tqdm import tqdm

binary_to_bp = {
    '00': 'A',
    '01': 'C',
    '10': 'G',
    '11': 'T'
}

bp_to_binary = {v: k for k, v in binary_to_bp.items()}


def int_to_sequence(val, length):
    bstring = format(val, '0' + str(length) + 'b')
    bstring = bstring[-1 * length:]
    twobits = [bstring[i:i + 2] for i in range(0, length, 2)]
    seq = [binary_to_bp[x] for x in twobits]
    return "".join(seq)


def sequence_to_int(seq):
    units = [bp_to_binary[x] for x in seq]
    bseq = ''.join(units)
    return int(bseq, 2)


@functools.lru_cache(maxsize=None)
def int_to_barcode(val):
    return int_to_sequence(val, 32)


def int_to_umi(val):
    return int_to_sequence(val, 20)


raw_molecule_file = snakemake.input["raw_molecule_file"]
filtered_molecule_file = snakemake.input["filtered_molecule_file"]
# e.g. molecule_file = "cellranger_out/outs/filtered_molecules.h5"

out_file = snakemake.output["out"]
out_pdf_file = snakemake.output["out_pdf"]

f = h5py.File(raw_molecule_file)

# Groups in the file (via f.visit(print))
#   barcode
#   barcode_corrected_reads
#   conf_mapped_uniq_read_pos
#   gem_group
#   gene
#   gene_ids
#   gene_names
#   genome
#   genome_ids
#   metrics
#   nonconf_mapped_reads
#   reads
#   umi
#   umi_corrected_reads
#   unmapped_reads

# barcode = f['barcode']  # now barcode is a dataset
# all_barcode = barcode[()] # Selects all data from barcode

# Throw it all into a giant dataframe??  Heck yeah!

barcode = f['barcode'][()]
gem = f['gem_group'][()]

# convert barcode+gem into barcode-gem
print("Converting barcodes")
barcode = [int_to_barcode(x)+'-'+str(y) for x, y in tqdm(zip(barcode, gem))]
barcode = pd.Series(barcode, name='barcode')

print("Loading data from raw file")
umi = pd.Series(f['umi'][()], name='umi')
gene = pd.Series(f['gene'][()], name='gene')
bcr = pd.Series(f['barcode_corrected_reads'][()], name='barcode_corrected')
ucr = pd.Series(f['umi_corrected_reads'][()], name='umi_corrected')
unmapped = pd.Series(f['unmapped_reads'][()], name='unmapped_reads')
nonconf = pd.Series(f['nonconf_mapped_reads'][()], name='genome_not_gene')
reads = pd.Series(f['reads'][()], name='reads')

f.close()

data = pd.concat((barcode, umi, gene, bcr, ucr,
                  unmapped, nonconf, reads), axis=1)

ff = h5py.File(filtered_molecule_file)
fbar = ff['barcode'][()]
fgem = ff['gem_group'][()]
ff.close()

# convert barcode+gem into barcode-gem
print("Converting filtered barcodes")
ffbarcode = {int_to_barcode(x)+'-'+str(y) for x, y in zip(fbar, fgem)}
valid_barcodes = [x in ffbarcode for x in data['barcode']]

data = data.loc[valid_barcodes]

print("Computing metrics")
# reads mapped to gene (proportion proportion per barcode)
mapped_per_cell = data.groupby("barcode")["reads"].sum()
mapped_per_cell.name = "mapped_reads"

# reads mapped to genome (proportion per barcode)
genome_not_gene = data.groupby("barcode")["genome_not_gene"].sum()

# Unmapped reads per cell (proportion per barcode)
unmapped = data.groupby("barcode")["unmapped_reads"].sum()

# Total reads per barcode
reads_per_barcode = mapped_per_cell + genome_not_gene + unmapped
reads_per_barcode.name = "num_reads"

# Corrected UMIs per cell (proportion)
umi_corrected = data.groupby("barcode")["umi_corrected"].sum()

# Corrected Cell barcode per cell (proportion)
barcode_corrected = data.groupby("barcode")["barcode_corrected"].sum()

# Turn raw numbers into proportion:
sub_barcode_data = pd.concat(
    (mapped_per_cell, genome_not_gene, unmapped,
     umi_corrected, barcode_corrected),
    axis=1)

sub_barcode_data = sub_barcode_data.divide(reads_per_barcode, axis=0)

# UMIs per cell
umi_per_cell = data.groupby("barcode")["umi"].nunique()
umi_per_cell.name = "num_umi"

# Mean reads per UMI and Std reads per UMI
reads_per_umi = data.groupby(["barcode", "umi"])[
    "reads", "genome_not_gene", "unmapped_reads"].sum() \
    .sum(axis=1).reset_index()

mean_reads_per_umi = reads_per_umi.groupby("barcode")[0].mean()
mean_reads_per_umi.name = "mean_reads_per_umi"

std_reads_per_umi = reads_per_umi.groupby("barcode")[0].std()
std_reads_per_umi.name = "std_reads_per_umi"


x = (data.groupby(["barcode", "gene"])["reads"]).sum()
x = x[x > 0]
num_genes_detected = x.groupby("barcode").size()
num_genes_detected.name = "num_genes_detected"

barcode_data = pd.concat((umi_per_cell, reads_per_barcode,
                          mean_reads_per_umi, std_reads_per_umi,
                          sub_barcode_data, num_genes_detected), axis=1)


print("Generating plot")
# Plot histograms of all the things
plt.figure(figsize=(10, 10))
for i, col in enumerate(barcode_data.columns):
    ax = plt.subplot(4, 3, i + 1)
    rmin = barcode_data[col].quantile(.05)
    rmax = barcode_data[col].quantile(.95)
    rrange = rmax - rmin
    rmin -= rrange * .1
    rmax += rrange * .1

    plt.hist(barcode_data[col], bins=20, range=(rmin, rmax))

    plt.title(col)

    if i % 3 == 0:
        plt.ylabel("# Barcodes")

plt.subplots_adjust(wspace=.5, hspace=.4)
plt.suptitle("QC Metrics")
plt.savefig(out_pdf_file)


# Write initial set of QC to file
print("Saving metrics")
barcode_data.to_csv(out_file, sep="\t", compression="gzip")


# Add some more quality metrics
# -- Mean % GC
# -- Mean Gene Length
# -- Skewness

# Data is in sparse format - convert to regular pandas dataframe

#  m_file = os.path.join(data_dir, "filtered_gene_bc_matrices/mm10/matrix.mtx")
#  b_file = os.path.join(data_dir, "filtered_gene_bc_matrices/mm10/barcodes.tsv")
#  g_file = os.path.join(data_dir, "filtered_gene_bc_matrices/mm10/genes.tsv")
#
#
#  gx = scipy.io.mmread(m_file)
#  gx = gx.toarray()
#
#  barcodes = pd.read_table(b_file, header=None)
#  barcodes = [x for x in barcodes[0]]
#
#  genes = pd.read_table(g_file, header=None)
#  gene_ids = genes[0].tolist()
#
#  gx = pd.DataFrame(gx, index=gene_ids, columns=barcodes)
#  barcode_data = barcode_data.loc[gx.columns]
#
#
#  # Load the reference gene_meta table for this purpose
#  gene_meta = pd.read_table("gene_meta.txt", index_col=0)
#
#  # collapse expression to gene symbols and save
#  gx_sym = gx.join(gene_meta.Name.str.upper()).groupby('Name').sum()
#  gx_sym.to_csv("gene_expression.txt", sep="\t")
#
#  # Compute weighted average (by count) of gene length per sample
#  num = gx.multiply(gene_meta['Length'], axis=0).sum(axis=0)
#  denom = gx.sum(axis=0)
#  gene_length = num / denom
#
#  # Compute weighted average (by count) of gene GC content per sample
#  num = gx.multiply(gene_meta['GC'], axis=0).sum(axis=0)
#  denom = gx.sum(axis=0)
#  gene_gc = num / denom
#
#  # Compute the Skewness of each sample
#
#
#  def skew_fun(sample, N):
#      total = sample.sum()
#      topN = sample.sort_values()[-1 * N:].sum()
#      return topN / total
#
#
#  sample_skew = gx.apply(skew_fun, axis=0, reduce=True, args=(50,))
#
#  # add these metrics to the barcode data we have already
#  gene_length.name = 'gene_length'
#  gene_gc.name = 'gc'
#  sample_skew.name = 'skew'
#
#  barcode_data2 = pd.concat(
#      (barcode_data, gene_length, gene_gc, sample_skew), axis=1)
#
#
#  # Plot the new metrics
#  plt.figure(figsize=(4, 10))
#  for i, col in enumerate(['gene_length', 'gc', 'skew']):
#      ax = plt.subplot(3, 1, i + 1)
#      rmin = barcode_data2[col].quantile(.05)
#      rmax = barcode_data2[col].quantile(.95)
#      rrange = rmax - rmin
#      rmin -= rrange * .1
#      rmax += rrange * .1
#
#      plt.hist(barcode_data2[col], bins=20, range=(rmin, rmax))
#
#      plt.title(col)
#
#      if i % 3 == 0:
#          plt.ylabel("# Barcodes")
#
#  plt.subplots_adjust(wspace=.5, hspace=.4)
#  plt.suptitle("Extra QC Metrics")
#  plt.savefig("molecule_qc_extra.pdf")
#
#
#  barcode_data2.to_csv("molecule_qc2.txt", sep="\t")
