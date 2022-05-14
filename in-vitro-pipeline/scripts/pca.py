import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import os
import subprocess

repo_dir = os.path.dirname(subprocess.check_output(['git', 'rev-parse', '--git-dir'])).decode()


exp_file = snakemake.input['exp']
out_file = snakemake.output['out']

os.makedirs(os.path.dirname(out_file), exist_ok=True)

exp = pd.read_table(exp_file, index_col=0)

# Remove genes that are not protein coding

gb = pd.read_table(
    os.path.join(repo_dir, "Signatures/mm10_il23/transcript_dict.txt"),
    comment='#'
)

gb_map = {x: y for x, y in zip(gb['Gene Symbol'], gb['Biotype'])}

keep_gene = [x in gb_map and gb_map[x] == 'protein_coding' for x in exp.index]

print(exp.shape)
exp = exp.loc[keep_gene, :]
print(exp.shape)


# Drop variables with all zeros

exp = exp.loc[
    exp.std(axis=1) > 0, :
]

exp = np.log2(exp + 1)

exp_std = exp.subtract(
    exp.mean(axis=1), axis=0
).divide(
    exp.std(axis=1), axis=0
)

model = PCA(n_components=30)
model.fit(exp_std.values.T)

out = model.transform(exp_std.values.T)

out = pd.DataFrame(
    out, index=exp.columns,
    columns=['PC{}'.format(i+1) for i in range(out.shape[1])]
)

# Only take the top 15 components

out = out.iloc[:, 0:15]

out.to_csv(out_file, sep="\t", compression="gzip")
