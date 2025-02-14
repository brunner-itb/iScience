import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks


hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
path = "/extra2/brunner/paper_models/boxed_static/mesh_resolution/"
timing_df = pd.read_hdf(path + "timing_df" + "" + ".h5", mode="r")
cell_df = pd.read_hdf(path + "cell_df" + "" + ".h5", mode="r")
global_df = pd.read_hdf(path + "global_df" + "" + ".h5", mode="r")
cell_df["IL-2_surf_c"] *= 1e3
fig_path = "/home/brunner/Documents/Current work/2023_11_03/"  + "FigS1D_mesh_scaling.pdf"


# global_df = activation_to_global_df(cell_df, global_df)
#%%
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
run_df = timing_df.loc[timing_df["name"] == "step"]
# run_df["duration"] /= 60
resolution = [1630021, 620812, 324041, 168439, 85864, 73970, 63917, 34758, 33504, 29973]
run_df["resolution"] = np.array([[a] * 10 for a in resolution]).flatten()
sns.lineplot(data=run_df, x = "resolution", y = "duration", color="black")
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("mesh vertices")
plt.ylabel("computation time (s)")
plt.yscale("linear")
plt.xscale("log")
# plt.gca().invert_xaxis()
# plt.xticks([10, 6, 2], ["10", "6", "2"])
plt.yticks([0, 75, 150])


fig.savefig(fig_path, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

# sns.lineplot(data=cell_df.loc[cell_df.type_name == "Th"], x = "scan_value", y = "IL-2_surf_c", legend=False, estimator="mean")
# plt.axhline(1.6917890859614957, color="red")
# plt.tight_layout()
# plt.show()

#%%
path_2 = "/extra2/brunner/paper_models/statics/cells_scan/0.1_Tsecs/"
timing_df_2 = pd.read_hdf(path_2 + "timing_df" + "" + ".h5", mode="r")

fig_path = "/home/brunner/Documents/Current work/2023_11_03/"  + "FigS1D_mesh_size_scaling.pdf"

# rc_ticks["figure.figsize"] = [1.5, 1.11]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
run_df = timing_df_2.loc[timing_df_2["name"] == "step"]
run_df["cells"] = np.array([[a] * 10 for a in np.arange(2,13)**3]).flatten()
# run_df["duration"] /= 60
sns.lineplot(data=run_df, x = "cells", y = "duration", estimator="median", color="black")
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel(r"system size (# cells)")
plt.ylabel("computation time (s)")
plt.yscale("linear")
# plt.ylim((0, 2.5))
# plt.gca().invert_xaxis()
plt.ylim([0,30])
plt.yticks([0, 15, 30])
# plt.xticks([0,5, 10], ["small", "standard", "large"])


fig.savefig(fig_path, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
