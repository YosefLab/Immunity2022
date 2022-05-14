import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

tpm_in = snakemake.input['tpm']
counts_in = snakemake.input['counts']
qc_in = snakemake.input['qc']
meta_in = snakemake.input['meta']

tpm_out = snakemake.output['tpm']
try:
    tpm_out_all = snakemake.output['tpm_all']
except AttributeError:
    tpm_out_all = None

counts_out = snakemake.output['counts']
qc_out = snakemake.output['qc']
meta_out = snakemake.output['meta']

# Load data
tpm = pd.read_table(tpm_in, index_col=0)
counts = pd.read_table(counts_in, index_col=0)
qc = pd.read_table(qc_in, index_col=0)
meta = pd.read_table(meta_in, index_col=0)


# Filter data
def filter_data(tpm, counts, qc, meta, valid_samples=None, valid_genes=None):

    if valid_samples is not None:
        tpm = tpm[valid_samples]
        counts = counts[valid_samples]
        qc = qc.loc[valid_samples]
        meta = meta.loc[valid_samples]

    if valid_genes is not None:
        tpm = tpm.loc[valid_genes]
        counts = counts.loc[valid_genes]

    return tpm, counts, qc, meta


try:
    stimulus = snakemake.params['stimulus']
    valid_samples = meta.index[meta.Stimulus == stimulus]
    tpm, counts, qc, meta = filter_data(
        tpm, counts, qc, meta, valid_samples=valid_samples)
except AttributeError:
    pass


# First filter out wells G1 and H1
# G1 tends to be no samples, and H1 is, maybe, population data
valid_samples = meta.index[
    (meta.Well != "H1") &
    (meta.Well != "G1")
]
tpm, counts, qc, meta = filter_data(
    tpm, counts, qc, meta, valid_samples=valid_samples)

# Threshold NREADS at 800,000 (min) and 3M (max) reads

valid_samples = qc.index[
    (qc['NREADS'] > 800e3) & (qc['NREADS'] < 3e6)
]
tpm, counts, qc, meta = filter_data(
    tpm, counts, qc, meta, valid_samples=valid_samples)

# RALIGN at 40%
valid_samples = qc.index[(qc['RALIGN'] > 40)]
tpm, counts, qc, meta = filter_data(
    tpm, counts, qc, meta, valid_samples=valid_samples)

# tpm out all
if tpm_out_all is not None:
    tpm.to_csv(tpm_out_all, sep="\t", compression="gzip")

# filter genes, need detection (count > 0) in at least 5 samples
valid_genes = counts.index[(counts > 0).sum(axis=1) >= 5]
tpm, counts, qc, meta = filter_data(
    tpm, counts, qc, meta, valid_genes=valid_genes)

# Save results
meta.to_csv(meta_out, sep="\t", compression="gzip")
counts.to_csv(counts_out, sep="\t", compression="gzip")
tpm.to_csv(tpm_out, sep="\t", compression="gzip")
qc.to_csv(qc_out, sep="\t", compression="gzip")
