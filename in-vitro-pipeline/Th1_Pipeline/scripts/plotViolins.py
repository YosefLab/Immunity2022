import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np

plt.rcParams['svg.fonttype'] = 'none'

tpm = pd.read_table("data_filtered/tpm.txt.gz", index_col=0)
meta = pd.read_table("data_filtered/meta.txt.gz", index_col=0)

# Re-code some things to make them look nicer
meta.GFP = meta.GFP.map({'eGFP_pos': 'GFP+', 'eGFP_neg': 'GFP-'})

logtpm = np.log2(tpm + 1)

all_dat = logtpm.T.join(meta)

all_dat_ctrl = all_dat.loc[all_dat.Genotype == "ctrl"]



# %%
gene = 'TIMP1'

plt.figure()

palette = sns.color_palette()
del palette[1]

sns.violinplot(data=all_dat_ctrl, x='GFP', y=gene, cut=0, inner=None,
               palette=palette, linewidth=.5)
sns.stripplot(data=all_dat_ctrl, x='GFP', y=gene,
              color='#111111', dodge=True, size=2)

# Lines for mean
mw = .1
pos_mean = all_dat_ctrl[gene][all_dat_ctrl.GFP == "GFP+"].mean()
neg_mean = all_dat_ctrl[gene][all_dat_ctrl.GFP == "GFP-"].mean()
plt.hlines(y=pos_mean, xmin=1-mw, xmax=1+mw)
plt.hlines(y=neg_mean, xmin=0-mw, xmax=0+mw)

plt.ylabel('$log_2(TPM+1)$')
plt.xlabel('')
plt.title(gene)
plt.show()

# %% Ok, now do this in a grid
genes = [
    'IL23R',
    'IFNG',
    'CSF1',
    'IL22',
    'IL18R1',
    'IL2RA',
    'IL1R2',
    'SOCS3',
    'CD3E',
    'LGALS1',
    'CD81',
    'S1PR1',
    'STAT5A',
    'BATF',
    'RORC',
    'RORA',
]


with mpl.rc_context(rc={
        'axes.titlesize': 12,
        'axes.titlepad': 4,
        'axes.labelsize': 9,
        'xtick.labelsize': 7,
        'xtick.direction': 'inout',
        'ytick.direction': 'in',
        'ytick.major.size': 4,
        'ytick.major.pad': 1,
        'ytick.labelsize': 7,
}):

    fig, axs = plt.subplots(4, 4, figsize=(9, 9))

    palette = sns.color_palette()
    del palette[1]

    i = 0
    for ax, gene in zip(axs.ravel(), genes):

        plt.sca(ax)

        sns.violinplot(data=all_dat_ctrl, x='GFP', y=gene, cut=0, inner=None,
                       palette=palette, linewidth=.5)
        sns.stripplot(data=all_dat_ctrl, x='GFP', y=gene,
                      color='#222222', dodge=True, size=1, rasterized=True)

        data_max = all_dat_ctrl[gene].max()
        # Lines for mean
        mw = .1
        pos_mean = all_dat_ctrl[gene][all_dat_ctrl.GFP == "GFP+"].mean()
        neg_mean = all_dat_ctrl[gene][all_dat_ctrl.GFP == "GFP-"].mean()
        plt.hlines(y=pos_mean, xmin=1-mw, xmax=1+mw)
        plt.hlines(y=neg_mean, xmin=0-mw, xmax=0+mw)

        plt.xlabel('')
        plt.ylabel('')
        plt.title(gene)

        if data_max < 10:
            plt.yticks(list(range(0, int(data_max+1), 2)))
        else:
            plt.yticks(list(range(0, int(data_max+1), 3)))

        if i < 12:
            ax.set_xticklabels(["", ""])

        if i % 4 == 0:
            plt.ylabel('$log_2(TPM+1)$')

        i = i + 1

    plt.subplots_adjust(wspace=.3, hspace=.4)
    #plt.show()
    plt.savefig("Figures/GFP_Violins.svg", dpi=300)

# %% Or do it using the two colors as a hue variable


def plot_genes(genes, ax):
    dsub = all_dat_ctrl.loc[:, genes + ['GFP']]
    dsub = dsub.melt(id_vars='GFP', var_name='Gene', value_name='Expression')

    sns.violinplot(
        data=dsub, x='Gene', y='Expression', hue='GFP',
        palette=palette, cut=0, inner=None
    )

    sns.stripplot(data=dsub, x='Gene', y='Expression', hue='GFP',
                  color='#222222', dodge=True, size=1)

    plt.ylabel('$log_2(TPM+1)$')
    plt.xlabel('Gene')
    plt.xticks(rotation=75)

    # Get the handles and labels. For this example it'll be 2 tuples
    # of length 4 each.
    handles, labels = ax.get_legend_handles_labels()

    # When creating the legend, only use the first two elements
    # to effectively remove the last two.
    l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1),
                   loc='upper left', borderaxespad=0.)

    ax.set_axisbelow(True)
    plt.grid(color='#888888', linestyle='--', linewidth=.5, axis='y')

fig, axs = plt.subplots(2, 1, figsize=(6, 6))
plt.sca(axs[0])
plot_genes(genes[:8], axs[0])
plt.sca(axs[1])
plot_genes(genes[8:], axs[1])

plt.tight_layout()
plt.show()

# %% Same thing but horizontal

def plot_genes(genes, ax):
    dsub = all_dat_ctrl.loc[:, genes + ['GFP']]
    dsub = dsub.melt(id_vars='GFP', var_name='Gene', value_name='Expression')

    sns.violinplot(
        data=dsub, y='Gene', x='Expression', hue='GFP',
        palette=palette, cut=0, inner=None, orient='h'
    )

    sns.stripplot(data=dsub, y='Gene', x='Expression', hue='GFP',
                  color='#222222', dodge=True, size=1, orient='h')

    plt.xlabel('$log_2(TPM+1)$')
    plt.ylabel('Gene')

    # Get the handles and labels. For this example it'll be 2 tuples
    # of length 4 each.
    handles, labels = ax.get_legend_handles_labels()

    # When creating the legend, only use the first two elements
    # to effectively remove the last two.
    l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(0, 1.05),
                   loc='lower left', borderaxespad=0.)

    ax.set_axisbelow(True)
    plt.grid(color='#888888', linestyle='--', linewidth=.5, axis='x')

fig, axs = plt.subplots(1, 2, figsize=(6, 6))
plt.sca(axs[0])
plot_genes(genes[:8], axs[0])
plt.sca(axs[1])
plot_genes(genes[8:], axs[1])

plt.tight_layout()
plt.show()
