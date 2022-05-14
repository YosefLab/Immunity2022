import pandas as pd
from sklearn.manifold import TSNE


pc_data = pd.read_table('../pca/PCS.txt', index_col=0)


model = TSNE(n_components=2, verbose=1, random_state=0)

result = model.fit_transform(pc_data)

result = pd.DataFrame(result, index=pc_data.index, columns=['tsne1', 'tsne2'])

result.to_csv('tsne.txt', sep='\t')
