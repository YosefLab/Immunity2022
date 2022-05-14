import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib import gridspec, cm
import numpy as np
import pandas as pd
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

# Make the effect-size plot

de_ctrl_file = "DE_Th1_ctrl_noqc/results.txt.gz"
de_ko_file = "DE_Th1_KO_noqc/results.txt.gz"

de_ctrl = pd.read_table(de_ctrl_file) \
    .query("contrast == 'GFPeGFP_pos'") \
    .set_index('primerid')

de_ko = pd.read_table(de_ko_file) \
    .query("contrast == 'GFPeGFP_pos'") \
    .set_index('primerid')

print((de_ctrl.FDR < 0.1).sum())
print((de_ko.FDR < 0.1).sum())

x1 = de_ctrl.FDR.sort_values().values
x2 = de_ko.FDR.sort_values().values

y = np.arange(x1.size)+1

# %%

fig = plt.figure(figsize=(5, 4))

plt.plot(
    np.log10(x1), y,
    label='$IL23R^{WT/eGFP}$'
)

plt.plot(
    np.log10(x2), y,
    label='$IL23R^{eGFP/eGFP}$'
)

plt.xlim(-3, -.5)
plt.ylim(-10, 1500)

plt.vlines(-1, ymin=-1e6, ymax=1e6, colors='#444444')
plt.legend(loc='upper left')

plt.xlabel('FDR ($log_{10}$)')
plt.ylabel('# of Differential Genes')

plt.tick_params(labelsize=8)

plt.subplots_adjust(bottom=0.15, left=0.2)

fig.patch.set_visible(False)

ax = plt.gca()
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.grid(linestyle='--', color='#888888', linewidth=.5)

# plt.show()

plt.savefig('Figures/EffectSize.svg')
