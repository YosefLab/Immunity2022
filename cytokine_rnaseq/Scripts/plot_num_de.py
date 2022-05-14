import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations, product

plt.rcParams['svg.fonttype'] = 'none'

data = pd.read_table("pairwise_de.txt")

de_counts = (data
        .loc[lambda x: (x['FDR'] < .1) & (x['logFC'].abs() > .5)]
        .groupby(['NumComparison', 'DenomComparison'])
        .size()
        .reset_index()
)

groups = [
    ("IL-1b+IL-6+IL-23", "TGFb+IL-6"),

    ("TGFb+IL-6", "IL-12"),
    ("TGFb+IL-6", "IL-12+IL-21"),

    ("IL-1b+IL-6+IL-23", "IL-12"),
    ("IL-1b+IL-6+IL-23", "IL-12+IL-21"),

    ("IL-12+IL-21", "IL-12"),
]


plot_data = []
for group_num, group_denom in groups:
    label = group_num + " vs. " + group_denom
    val = de_counts.loc[
        lambda x: (x['NumComparison'] == group_num) & (x['DenomComparison'] == group_denom),
        0
    ].iloc[0]

    plot_data.append(
        [group_num, group_denom, label, val]
    )

plot_data = pd.DataFrame(plot_data, columns=['Num', 'Denom', 'Label', 'Value'])

# %%

plt.figure(figsize=(5, 4))
plt.barh(
    y=plot_data['Label'],
    width=plot_data['Value'],
    height=0.5,
    color='#333333',
)
plt.xlabel('# Differential Genes\n(FDR < .1 and $|log_2FC|$ > .5)')
plt.yticks(size=9)
plt.grid(True, axis='x', color='#CCCCCC')
plt.gca().set_axisbelow(True)
plt.subplots_adjust(left=0.5, bottom=0.2)

plt.savefig('Figures/cytokine_de_counts.svg')
