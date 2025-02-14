import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
from copy import deepcopy
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.statics.figure_2.plotting.funcs_for_plotting import get_path, scale_dataframes
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation

'''
define which static sim you want to plot. 
bc defines the boundary condition, scan variable the scan.
"IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan and IL-2_KD the KD scan.
'''

# bc = "standard"
bc = "saturated"

scan_variables = ["IL-2_Tsec_fraction", "IL-2_sigma", "IL-2_KD"]

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/" + "full_SD_{m}_static_{sv}_{bc}".format(m = "box", sv=scan_variables[0], bc = bc) + ".pdf"

sv_styles = ['-','-.',':']
plot_names = ["sec. scan", "R scan", "KD scan"]
# define the standard values (fold-change = 0) for each of the scan_variables
standards = [0.05, 1, 7.4]

# which models to plot

# models = ["loc. q and R", "loc. q", "well_mixed"]
models = ["loc. q and R"]
model_styles = ["-", "--"]
# Which colors to use for each scan_variable.
model_colors = [["red", "green", "blue"]]
model_alphas = [1,0.5]

plot_legend = True


# plotting parameters

yscale = "linear"
yscale_base = None
xscale = "log"
xscale_base = 2
# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum value
max_value = 1e5
min_value = -1

xlim = (2 ** -2.1, 2 ** 2.1)
ylim = (None, None)

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0
# which cell type to plot in which colour. The second colour is the SD colour
cell_types = ["Th"]

scan_measure = "IL-2_surf_c"

replicat = "replicat_index"

# if line smoothing is desired, sfs are the smoothing factors
interpolate = False
sfs = [0.002, 0.01, 0.01]

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

# load runs, apply offset, merge dataframes
dataframes = [[] for x in models]
print("loading data")
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
for m,model in enumerate(models):
    for sv, scan_variable in enumerate(scan_variables):
        m_sv_path = get_path(bc, hdd, model, scan_variable)
        c_df, g_df = my_load_df(m_sv_path, offset = 0, run_range = [0], custom_ending = "_combined")
        try:
            c_df["IL-2_Tsec_fraction"]
        except KeyError:
            c_df["IL-2_Tsec_fraction"] = c_df["fractions_Tsec"]
            try:
                g_df["IL-2_Tsec_fraction"] = g_df["fractions_Tsec"]
            except:
                pass
        try:
            c_df["IL-2_sigma"]
        except KeyError:
            c_df["IL-2_sigma"] = c_df["misc_sigma"]
            try:
                g_df["IL-2_sigma"] = g_df["scan_value"]
            except:
                pass
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, min_value, max_value, standards)
        # if scan_variable == 'IL-2_KD':
        #     c_df = c_df.loc[c_df["IL-2_KD"].isin(c_df["IL-2_KD"].unique()[5:])]
        #     g_df = g_df.loc[g_df["IL-2_KD"].isin(c_df["IL-2_KD"].unique()[5:])]
        # print(c_df["IL-2_KD"].unique())
        dataframes[m].append([c_df, g_df])

#%%
for c_df, g_df in dataframes[0]:
    print(len(c_df["replicat_index"].unique()))
#%%
########################################################################################################################
print("plotting")
########################################################################################################################
# fc_to_plot = [2**-1, 2**0, 2**1]
fc_to_plot = [1/10, 1/5, 1/2, 1, 2, 5, 10]
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()
df_list = []
for m,model in enumerate(models):
    x = np.arange(0, len(scan_variables) * 2)
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]

        ste = []
        cv = []
        for f, fc in enumerate(fc_to_plot):
            scan = c_df[scan_variable].unique()[(np.abs(c_df[scan_variable].unique() - (standards[sv] * fc))).argmin()]
            sliced_df = c_df.loc[(c_df["type_name"] == "Th") & (c_df[scan_variable] == scan)]
            run_cv = []
            mean = []
            std = []
            for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                std.append(sliced_df.loc[(sliced_df[replicat] == run), scan_measure].std())
                mean.append(sliced_df.loc[(sliced_df[replicat] == run), scan_measure].mean())
                run_cv.append(std[-1] / mean[-1])
                df_list.append({"model": model, "scan_variable": scan_variable, "fold_change": fc, "replicat": run,
                                "cv": run_cv[-1], "std":std[-1]})
            ste.append(np.std(run_cv))# / np.sqrt(len(g_df[replicat].unique())))
            cv.append(np.mean(run_cv))
plot_df = pd.DataFrame(df_list)
#%%
rc_ticks["figure.figsize"] = (1.7, 1.2)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()
from matplotlib.colors import to_rgba
base_colors = ["black", "black"]
base_alphas = np.logspace(np.log10(0.15),0,len(fc_to_plot))
model_colors = []
for bc in base_colors:
    model_colors.append([to_rgba(bc, alpha) for alpha in base_alphas])

xs = np.array(np.arange(len(plot_df.fold_change.unique())))
for sv, scan in enumerate(plot_df.scan_variable.unique()):
    scan_df = plot_df.loc[(plot_df.scan_variable == scan)]
    moved_xs = xs + (np.max(xs) + 2) * sv
    ys = scan_df.groupby("fold_change").mean()["std"].values
    yerr = scan_df.groupby("fold_change").std()["std"].values
    plt.bar(moved_xs, ys, yerr=yerr, linewidth=0.2, width=0.8, color=model_colors[m], edgecolor="black", error_kw=dict(lw=0.3))


# ax.set(ylabel=r"log$_2$ FC (surf. conc. CV)", xlabel="", ylim=ylim)
ax.set(ylabel=r"surface conc. s.d. (pM)", xlabel="", yticks=[0, 10, 20], yticklabels=[0, 10, 20], ylim=(0, 20))
plt.xticks([np.mean(xs),np.mean([np.mean(moved_xs), np.mean(xs)]),np.mean(moved_xs)],
           [r"f$_{\rm sec}$", r"$\sigma$", r"K$_{\rm D}$"])
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()