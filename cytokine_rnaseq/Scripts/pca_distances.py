import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations, product


meta = pd.read_table("rnaseq/meta.txt", index_col=0)
pca = pd.read_csv("pca.csv", index_col=0)

N_COMPONENTS = 2

pca = pca.iloc[:, :N_COMPONENTS]
meta = meta.loc[pca.index]


def dist_for_samples(sample1, sample2):
    s1_vec = pca.loc[sample1]
    s2_vec = pca.loc[sample2]

    return ((s1_vec - s2_vec)**2).sum()**.5


results = []
for condition in meta.Condition.unique():
    samples = meta.loc[lambda x: x['Condition'] == condition].index
    combos = combinations(samples, 2)
    for sample1, sample2 in combos:
        results.append(
            [
                condition + ' Replicates',
                dist_for_samples(sample1, sample2)
            ]
        )

for conditionA, conditionB in combinations(meta.Condition.unique(), 2):
    samplesA = meta.loc[lambda x: x['Condition'] == conditionA].index
    samplesB = meta.loc[lambda x: x['Condition'] == conditionB].index
    combos = product(samplesA, samplesB)
    for sample1, sample2 in combos:
        results.append(
            [
                conditionA + ' vs. ' + conditionB,
                dist_for_samples(sample1, sample2)
            ]
        )

results_df = pd.DataFrame(results, columns=['Group', 'Dist'])

# %%

plt.figure()
sns.boxplot(
    data=results_df,
    y='Group', x='Dist'
)
plt.ylabel('')
plt.subplots_adjust(left=0.4)
plt.savefig('out.png')
