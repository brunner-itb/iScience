import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

from thesis.scripts.paper_models.utilities.plot_helper import EC50_calculation

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/" + "surf_c_over_q_lin_satu" + ".pdf"


path = "/extra2/brunner/paper_models/statics/saturated/Figure_2C/dataframes_q_scan_linear_steady_state/"
lin_cell_df = pd.read_hdf(path + "cell_df" + "" + ".h5", mode="r")
lin_cell_df["IL-2_surf_c"] *= 1e3

path = "/extra2/brunner/paper_models/statics/saturated/Figure_2C/dataframes_q_scan_satu_steady_state/q_scan_satu_sec_q_scanfth_0.9_fsec_0.1/"
satu_cell_df = pd.read_hdf(path + "cell_df" + "" + ".h5", mode="r")
satu_cell_df["IL-2_surf_c"] *= 1e3


#%%
rc_ticks["figure.figsize"] = (1.7, 1.2)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

color1 = "blue"
color2 = "black"
sns.lineplot(data=satu_cell_df, x = "scan_value", y = "IL-2_surf_c", label="saturated", color=color1)
sns.lineplot(data=lin_cell_df, x = "scan_value", y = "IL-2_surf_c", label="linear", color=color2)
plt.yscale("log")
plt.xscale("log")
ax.set(ylabel=r"surface conc. avg. (pM)", xlabel="cytokine secretion rate (1/s)")

import matplotlib
locmaj = matplotlib.ticker.LogLocator(base=10,numticks=1000)
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10,subs=np.arange(1,10),numticks=1000)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.yticks([1e0, 1e1, 1e2, 1e3, 1e4], [r"10$^0$", "", r"10$^2$", "", r"10$^4$"])
plt.legend().remove()

ax.tick_params(axis='y', colors=color1, which="both")
ax.yaxis.label.set_color(color1)

fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
