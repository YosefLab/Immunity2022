import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if not os.path.exists('gene_plots'):
    os.mkdir('gene_plots')


# Load expression

expression = pd.read_table(
    '../expression/expression_sym.txt.gz',
    index_col=0
)


# Load tsne coordinates
tsne = pd.read_table(
    'tsne.txt',
    index_col=0
)


key_genes = [
    'TBX21',
    'IL12RB2',
    'CCL5',
    'IFNG',
    'IL23R',
    'CXCR3',
    'STAT4',
    'RORC',
    'IL17A',
    'IL17F',
    'IL22',
    'CCR6',
    'STAT3',
    'CSF2',
    'FOXP3',
    'IL10',
    'ENTPD1'
]


def plot_gene(gene):

    x = tsne['tsne1']
    y = tsne['tsne2']
    color = np.log2(expression.loc[gene]+1)

    plt.figure()
    plt.scatter(x, y, c=color, s=4)
    cb = plt.colorbar()
    cb.set_label('$Log_2(CPM+1)$')
    plt.xlabel('TSNE-1')
    plt.ylabel('TSNE-2')
    plt.title(gene)

    plt.savefig('gene_plots/{}.png'.format(gene))
    plt.close()


for g in key_genes:
    plot_gene(g)
