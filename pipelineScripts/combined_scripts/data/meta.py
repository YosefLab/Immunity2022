import pandas as pd

agg_file = snakemake.input["agg_file"]
meta_file = snakemake.input["meta_file"]
bc_file = snakemake.input["expression"]
out_file = snakemake.output["out"]

agg = pd.read_csv(agg_file, index_col=0)

sample_map = {str(i+1): e for i, e in enumerate(agg.index)}

meta = pd.read_table(meta_file)


def bc_to_gem(bc):
    return bc[bc.index('-')+1:]


bc = pd.read_table(bc_file, nrows=1).columns[1:]

meta_out = pd.DataFrame({"Barcode": bc})
meta_out["gem"] = [bc_to_gem(x) for x in meta_out["Barcode"]]
meta_out["Sample"] = [sample_map[x] for x in meta_out["gem"]]
meta_out = meta_out.merge(meta, on="Sample")

meta_out = meta_out.drop("gem", axis=1).set_index("Barcode")

meta_out.to_csv(out_file, sep="\t", compression="gzip")
