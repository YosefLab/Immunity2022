"""
Assumes that de_results input file has columns
    - 'logFC', 'PValue', and 'name'

If split_by parameter is given, then input 'de_results' is run on this
and every split is processed separately.

Outputs are:
    out_dir/prerank_res_<split_name>.pkl   -- this is needed for plots later
    out_dir/res_table_<split_name>.txt
"""
import os
import gseapy
import numpy as np
import pandas as pd
import pickle

from pathlib import Path

de_results_file = snakemake.input["de_results"]
out_file_txt = snakemake.output['txt']


def read_gmt(gmt_file):
    gs = {}
    for line in open(gmt_file, "r"):
        sline = line.strip().split("\t")
        gs_name = sline[0]
        gs_genes = sline[2:]
        gs_genes = [x.upper() for x in gs_genes]
        gs[gs_name] = gs_genes
    return gs


de_results = pd.read_table(de_results_file, index_col=0)
de_results.index.name = 'name'
de_results = de_results.reset_index()

sig_files = [
    "../Signatures/MSIGDB/C7_IMMSIG_TH.gmt",
    "../Signatures/MSIGDB/H_Hallmark.gmt",
    "../Signatures/WikiPathways/WikiPathways_Mus_Musculus.gmt",
    "../Signatures/GO/GO_biological_process.gmt",
    "../Signatures/IBD_Colitis/Venema_Voskuil_2018.gmt",
]

gs_all = {}
for file in sig_files:
    gs = read_gmt(file)
    gs_all.update(gs.items())


rank_data = de_results.copy()
rank_data["rank"] = (
    np.sign(rank_data["logFC"]) * np.log10(rank_data["PValue"]) * -1
)
rank_data = rank_data[["name", "rank"]].sort_values("rank")

prerank_res = gseapy.prerank(
    rank_data,
    gene_sets=gs_all,
    processes=8,
    permutation_num=1000,
    no_plot=True,
    seed=5000,
)


# with open(out_file_pkl, 'wb') as fout:
#     pickle.dump(prerank_res, fout, pickle.HIGHEST_PROTOCOL)

res_table = prerank_res.res2d.sort_values("nes")

res_table.to_csv(out_file_txt, sep="\t")
