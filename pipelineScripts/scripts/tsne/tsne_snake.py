import pandas as pd
from MulticoreTSNE import MulticoreTSNE as TSNE
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pca_file = snakemake.input['pca']
output_file = snakemake.output['tsne']
output_image = snakemake.output['image']

out_dir = os.path.dirname(output_file)
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

N_JOBS = 4

pc_data = pd.read_table(pca_file, index_col=0)


model = TSNE(n_components=2, verbose=1, n_jobs=N_JOBS)

result = model.fit_transform(pc_data)

result = pd.DataFrame(result, index=pc_data.index, columns=['tsne1', 'tsne2'])

result.to_csv(output_file, sep='\t')

plt.figure()
plt.plot(result.tsne1, result.tsne2, 'o', ms=3)
plt.xlabel('TSNE-1')
plt.ylabel('TSNE-2')
plt.savefig(output_image)
