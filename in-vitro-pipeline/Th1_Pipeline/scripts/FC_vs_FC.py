import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ctrl = pd.read_table("DE_Th1_ctrl_noqc/results.txt.gz") \
    .groupby("contrast").get_group("GFPeGFP_pos") \
    .set_index("primerid")

ko = pd.read_table("DE_Th1_KO_noqc/results.txt.gz") \
    .groupby("contrast").get_group("GFPeGFP_pos") \
    .set_index("primerid")


to_keep = (
    ctrl.index[ctrl.FDR < .05] |
    ko.index[ko.FDR < .05]
)

plt.figure()

plt.plot(
    ctrl.loc[to_keep].logFC,
    ko.loc[to_keep].logFC,
    'o', ms=2
)
plt.show()


from bio_utils.plots import hover_plot
hover_plot(
    ctrl.loc[to_keep].logFC,
    ko.loc[to_keep].logFC,
    to_keep,
    'o', ms=2
)
plt.plot(
    [-3, 3], [-3, 3], '--'
)
plt.show()


# What if we color by interaction effect???

ix = pd.read_table("DE_Th1_full2_noqpc_relevel/out.txt.gz") \
    .groupby("contrast").get_group("Genotypectrl:GFPeGFP_pos") \
    .set_index("primerid")


plt.figure()
plt.scatter(
    ctrl.loc[to_keep].logFC,
    ko.loc[to_keep].logFC,
    s=4,
    c=np.log10(ix.loc[to_keep].FDR),
    vmin=-4, vmax=0, cmap='RdBu',
)
plt.colorbar()
plt.plot(
    [-3, 5], [-3, 5], '--', color='#bbbbbb'
)
plt.plot(
    [0, 0], [-3, 5], '-', color='black'
)
plt.plot(
    [-3, 5], [0, 0], '-', color='black'
)
plt.show()


# It doesn't look as good with the significance - maybe the color is better
# When using the interaction logFC


# %% Ok, make a plotly version?

from plotly.offline import plot, iplot
import plotly.graph_objs as go


trace = go.Scatter(
    x=ctrl.loc[to_keep].logFC,
    y=ko.loc[to_keep].logFC,
    text=to_keep,
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
            text="Ctrl: GFP+ vs GFP- (log2FC)",
        )
    ),
    yaxis=go.layout.YAxis(
        title=go.layout.yaxis.Title(
            text="KO: GFP+ vs GFP- (log2FC)",
        )
    ),
    shapes=[
        go.layout.Shape(
            type="line",
            x0=-3,
            y0=-3,
            x1=5,
            y1=5,
            line=dict(
                color="#cccccc",
                width=1
            )
        ),
    ]
)

fig = go.Figure([trace], layout)

plot(fig, filename='GFP+_vs_GFP-.html', auto_open=False)
