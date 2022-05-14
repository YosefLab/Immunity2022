import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib import gridspec, cm
import numpy as np
import pandas as pd
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

# %% Load Data
de_res_Xf = pd.read_table("DE_Th1_full2_noqpc_relevel/out.txt.gz")
de_res_X = (
    de_res_Xf.query('contrast == "Genotypectrl:GFPeGFP_pos_sub"')
             .set_index('primerid')
)

# %% Make Plot

genes = [
    'IL22', 'GZMB', 'CCR3', 'LAMP1',
    'ITGB1', 'GPR18', 'RAC1', 'PELI1',
    'MAPK1', 'MAP2K2', 'DUSP2', 'TRIM27',
    'JUND', 'LTA', 'BHLHE40',
]

sig_genes = de_res_X.index[de_res_X.FDR < 0.1]
nonsig_genes = de_res_X.index.difference(sig_genes)

fig = plt.figure()
plt.plot(
    de_res_X.logFC[sig_genes],
    np.log10(de_res_X.FDR[sig_genes])*-1,
    'o', ms=2
)
plt.plot(
    de_res_X.logFC[nonsig_genes],
    np.log10(de_res_X.FDR[nonsig_genes])*-1,
    'o', ms=3, color='#AAAAAA', alpha=0.5, mew=0
)
plt.hlines(1, -1000, 1000, colors='#000000', alpha=0.5, linestyles='dashed')
plt.xlim(-3, 3.5)
plt.ylim(-.5, 8)
plt.xlabel('Log Fold-Change')
plt.ylabel('FDR ($-log_{10}$)')

# Add annotations
ax = plt.gca()
for g in genes:
    x = de_res_X.logFC[g]
    y = np.log10(de_res_X.FDR[g])*-1
    ax.text(x, y, g.capitalize())

# plt.show()
fig.patch.set_visible(False)
plt.savefig('Figures/Interaction_Volcano.svg')

# %% Make an interactive version
from plotly.offline import plot, iplot
import plotly.graph_objs as go

trace = go.Scatter(
    x=de_res_X.logFC,
    y=np.log10(de_res_X.FDR)*-1,
    text=de_res_X.index,
    mode='markers'
)

layout = dict(
    title=go.layout.Title(
            text="GFP+ vs. GFP- in KO and Ctrl",
            xref="paper",
            x=0
        ),
    xaxis=go.layout.XAxis(
        title=go.layout.xaxis.Title(
            text="IL23RxGFP+ Effect (log2FC)",
        )
    ),
    yaxis=go.layout.YAxis(
        title=go.layout.yaxis.Title(
            text="IL23RxGFP+ FDR (-log10)"
        )
    ),
    shapes=[
        go.layout.Shape(
            type="line",
            x0=-2,
            y0=1,
            x1=3,
            y1=1,
            line=dict(
                color="#cccccc",
                width=1
            )
        ),
    ]
)

fig = go.Figure([trace], layout)

plot(fig, filename='Figures/Interaction_Volcano.html', auto_open=False)
