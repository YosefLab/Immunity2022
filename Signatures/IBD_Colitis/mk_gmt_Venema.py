"""
Parses signatures from the Venema_2019 supplement table
"""
import pandas as pd
from gene_enrich import GeneSet, write_gmt_file


sig_desc = "Source: PMID 30391472, Supp Table 1"
data = pd.read_table("Venema_2019_SuppTable1.txt")

sigs_dict = {}
groups = data.groupby(["Compartment", "Cell Type"])
for (comp, ct), dsub in groups:
    genes = set(dsub['Gene'])
    name = comp + "_" + ct
    name = name.replace("/", ";")
    sigs_dict[name] = genes


sigs_gs = [GeneSet(name, genes, description=sig_desc)
           for name, genes in sigs_dict.items()]

write_gmt_file(sigs_gs, "Venema_Voskuil_2018.gmt")
