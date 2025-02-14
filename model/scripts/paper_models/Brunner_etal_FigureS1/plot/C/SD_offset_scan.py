import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks


hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
path = "/extra2/brunner/paper_models/statics/offset_scan/test/"

cell_df = pd.read_hdf(path + "cell_df" + "" + ".h5", mode="r")
global_df = pd.read_hdf(path + "global_df" + "" + ".h5", mode="r")
cell_df["IL-2_surf_c"] *= 1e3
fig_path = "/home/brunner/Documents/Current work/2023_11_03/"  + "FigS1D_offset_scan.pdf"


# global_df = activation_to_global_df(cell_df, global_df)
#%%
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

sns.lineplot(data=cell_df.loc[cell_df.type_name == "Th"], x = "scan_value", y = "IL-2_surf_c", color="black")
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("offset layers")
plt.ylabel("surface conc. (pM)")
plt.yscale("linear")
# plt.xscale("log")
# plt.gca().invert_xaxis()
# plt.xticks([10, 6, 2], ["10", "6", "2"])
plt.xlim(0,3)
plt.ylim(3,9)
plt.yticks([3,6,9])

fig.savefig(fig_path, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

# sns.lineplot(data=cell_df.loc[cell_df.type_name == "Th"], x = "scan_value", y = "IL-2_surf_c", legend=False, estimator="mean")
# plt.axhline(1.6917890859614957, color="red")
# plt.tight_layout()
# plt.show()

