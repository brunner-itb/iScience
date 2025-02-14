import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
from brokenaxes import brokenaxes
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
rc_ticks["xtick.labelsize"] = 7
rc_ticks["ytick.labelsize"] = 6
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation
from thesis.scripts.paper_models.statics.figure_2.plotting.funcs_for_plotting import get_path, scale_dataframes

'''
define which static sim you want to plot. 
bc defines the boundary condition, scan variable the scan.
"IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan and IL-2_KD the KD scan.
'''

# bc = "standard"
bc = "saturated"

scan_variables = ["IL-2_Tsec_fraction"]
# scan_variables = ["IL-2_sigma"]


sv_styles = ['-','-.',':']
# plot_names = ["sec. scan", "R scan", "KD scan"]
plot_names = ["sec. scan"]
# plot_names = ["R scan"]
# define the standard values (fold-change = 0) for each of the scan_variables
standards = [0.1]
# standards = [1]

# which models to plot

# models = ["loc. q and R", "loc. q", "well_mixed"]
# models = ["well_mixed"]
models = ["well_mixed", "loc. q and R"]
model_styles = ["-", "--"]
# Which colors to use for each scan_variable.
model_colors = [["red", "green", "blue"]]
model_alphas = [1,0.5]

plot_legend = True

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/" + "pSTAT5_{model}_{m}_static_{sv}_{bc}".format(model = "ODE" if models[0] == "well_mixed" else "spatial", m = "box", sv=scan_variables[0], bc = bc) + ".pdf"

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
ylim = (-5, 5)

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
        # scales the dataframes. Depends on your desired units. Currently everything is in pM
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, min_value, max_value, standards)
        if scan_variable == 'IL-2_KD' and model != "well_mixed":
            c_df = c_df.loc[c_df["IL-2_KD"].isin(c_df["IL-2_KD"].unique()[5:])]
            g_df = g_df.loc[g_df["IL-2_KD"].isin(c_df["IL-2_KD"].unique()[5:])]
        dataframes[m].append([c_df, g_df])

#%%
# test_df = dataframes[1][0][0]
# # print(test_df["misc_KD"].unique())
# print(test_df["IL-2_KD"].unique())
# print(test_df.loc[test_df["scan_value"] == 1, "IL-2_R"].mean())
# print(test_df["IL-2_q"].unique())
# print(test_df["scan_value"].unique())
# print(test_df["fractions_Tsec"].unique())
#%%
########################################################################################################################
print("plotting")
########################################################################################################################
fc_to_plot = [1/10, 1/5, 1/2, 1, 2, 5, 10]
# fc_to_plot = [1/5, 1, 5]
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

df_list = []
for m,model in enumerate(models):
    x = np.arange(0, len(scan_variables) * 2)
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]

        ste = []
        act = []
        surf_c = []
        cv = []
        for f, fc in enumerate(fc_to_plot):
            scan = c_df[scan_variable].unique()[(np.abs(c_df[scan_variable].unique() - (standards[sv] * fc))).argmin()]
            sliced_df = c_df.loc[(c_df["type_name"] == "Th") & (c_df[scan_variable] == scan)]
            sliced_df["pSTAT5"] = sliced_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=0.8, R=sliced_df["IL-2_R"]) * 1e12) ** 3 +
                    sliced_df["IL-2_surf_c"] ** 3).values
            run_act = []
            try:
                for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                    run_df = sliced_df.loc[(sliced_df[replicat] == run)]
                    run_act.append(len(run_df.loc[(run_df["pSTAT5"] >= 0.5)].values)/len(run_df.values))
                    df_list.append({"model": model, "scan_variable": scan_variable, "fold_change": fc, "replicat": run, "activation": run_act[-1]})
            except KeyError: # ODE system does not have a gradient
                run_act = [0]

            ste.append(np.std(run_act) / np.sqrt(len(sliced_df[replicat].unique())))
            act.append(np.mean(run_act))
            surf_c.append(sliced_df["IL-2_surf_c"].mean())
            # print(scan_variable, surf_c)
            cv.append(sliced_df["IL-2_R"].std()/sliced_df["IL-2_R"].mean())
plot_df = pd.DataFrame(df_list)
# PROPS = {
#     'boxprops': {'facecolor': 'red', 'edgecolor': '0.1', 'alpha': 0.5},
#     'medianprops': {'color': '0.1'},
#     'whiskerprops': {'color': '0.1'},
#     'capprops': {'color': '0.1'}
# }
from copy import deepcopy
scaled_df = deepcopy(plot_df)
scaled_df["activation"] *= 100
scaled_df["activation_standard"] = scaled_df.loc[scaled_df["fold_change"] == 1, "activation"].mean()

scaled_df["act_change"] = (scaled_df["activation"] - scaled_df["activation_standard"])
scaled_df.loc[scaled_df["act_change"] == -np.inf, "act_change"] = 0



#%%
# rc_ticks["figure.figsize"] = (0.55, 1.2)
rc_ticks['xtick.major.size'] = 2.5
rc_ticks['xtick.major.pad'] = 2.5
rc_ticks['axes.titlesize'] = 6
rc_ticks['xtick.labelsize'] = 5.5
# rc_ticks["figure.figsize"] = (5, 5)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()
from matplotlib.colors import to_rgba
base_colors = ["black", "black"]
base_alphas = np.logspace(np.log10(0.15),0,len(fc_to_plot))
model_colors = []
for bc in base_colors:
    model_colors.append([to_rgba(bc, alpha) for alpha in base_alphas])

for sv, scan in enumerate(scaled_df.scan_variable.unique()):
    scan_df = scaled_df.loc[(scaled_df.scan_variable == scan)]
    xs = np.array(np.arange(len(scan_df.fold_change.unique())))
    for m, model in enumerate(scan_df.model.unique()):
        model_df = scan_df.loc[(scan_df.model == model)]
        xs += (np.max(xs) + 2)*m
        ys = model_df.groupby("fold_change").mean()["activation"].values
        yerr = model_df.groupby("fold_change").std()["activation"].values
        plt.bar(xs[:len(ys)], ys, yerr=yerr, linewidth=0.2, width=0.8, color=model_colors[m], edgecolor="black", error_kw=dict(lw=0.3))

# ax.set(ylim=(0, 3), ylabel=r"fold change %pSTAT5$^+$", xlabel="")
# ax.set(ylim=(-2.5, 4.5), ylabel=r"log$_2$ FC (%pSTAT5$^+$)", xlabel="")
# ax.set(ylabel=r"log$_2$ FC (%pSTAT$^+$)", xlabel="", ylim=ylim)
ax.set(ylabel=r"pSTAT$^+$ (%)", xlabel="", title=r"f$_{\rm sec}$: secreting cells", ylim=(0,None))
# plt.xticks([0.5,3.25], ["sec. cells", "R hetero."], rotation=-45, ha="left")
first = np.arange(len(scan_df.fold_change.unique()))[int(np.floor(len(scan_df.fold_change.unique())/2))]
second = xs[int(np.floor(len(xs)/2))]
plt.xticks([first, second], [r"well" "\n" r"mixed", "RD"])
# ax.get_legend().remove()
plt.yscale("linear")
ax.set_yticks([0, 6, 12])
ax.set_ylim(0, 12)
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
