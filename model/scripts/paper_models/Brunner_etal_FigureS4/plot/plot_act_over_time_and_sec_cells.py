import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df
from scipy.interpolate import UnivariateSpline
from thesis.main.my_debug import message

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
prefix = "neg"
# prefix = "pos"
# prefix = "ODE_s pos"
if prefix == "pos":
    model_name = "Tsec_scan_7"
    name = "dataframes_positive_fig3B_timeseries" #michaelis
    gamma = 12
    xlim = (None, 0.5)
elif prefix == "neg":
    model_name = "Tsec_scan_6"
    # name = "gamma_0.5"  # compute4
    # model_name = "feedback_scan_4"
    name = "dataframes_negative_timeseries" # menten
    gamma = 0.01


saving_string =r"/home/brunner/Documents/Current work/2023_06_16/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

filename = prefix + "_act_over_time_and_sec_cells"


plot_ODE = True
if prefix == "neg":
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/Tsec_scan_5/gamma_0.01/"
else:
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/Tsec_scan_5/gamma_100/"
ODE_cell_df =  pd.read_hdf(ODE_path + "cell_df" + ".h5", mode="r")
ODE_cell_df["IL-2_surf_c"] *= 1e3

yscale = "linear"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (None,None)
ylim = (-1, 101)

# colours, first = positive, second = negative
# colours = ["black", "blue"]
# opacity for the ODE
ODE_alpha = 0.25

c_interpolate = False
c_sf = 0
c_ste_sf = 0
sd_interpolate = False
sd_sf = 0
sd_ste_sf = 0


plot_every_cell = False
plot_std_area = True
# define which scan to plot
x_axis = "time"

# the first time index does not have pSTAT calculated yet.
skip_first_time_index = True
show_Ths = True
show_Tregs = False

plot_legend = True

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

base_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, extra=hdd, n=name)
# base_path = "/{extra}/{u}/paper_models/kinetics/{mn}/dataframes_positive_fig3B_timeseries/positive_fig3B_fsec_scangamma_12_ftreg_0_fth_0.999/".format(u=user, mn=model_name, extra=hdd, n=name)

# cell_df, global_df = my_load_df(base_path, offset=0, custom_ending="")
# base_path = "/extra2/brunner/paper_models/ODE/saturated/kinetics/Tsec_scan_5/gamma_100/"
spatial_cell_df, global_df =  my_load_df(base_path, offset=0, custom_ending = "_combined")
try:
    spatial_cell_df["IL-2_gamma"]
except KeyError:
    spatial_cell_df["IL-2_gamma"] = spatial_cell_df["misc_gamma"]
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]

spatial_cell_df["IL-2_surf_c"] *= 1e3
if skip_first_time_index == True:
    spatial_cell_df = spatial_cell_df.loc[spatial_cell_df["time_index"] != 0]

#%%
print(spatial_cell_df["time"].unique()/(3600 * 24))
print("plotting")
rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
g_cell_df = spatial_cell_df.loc[(spatial_cell_df["IL-2_gamma"] == gamma)]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()

for c,cell_df in enumerate([g_cell_df, ODE_cell_df]):
# for c,cell_df in enumerate([g_cell_df]):
    fractions = cell_df["IL-2_Tsec_fraction"].unique()[::2] if c == 0 else [0.102]
    # print(fractions)
    if prefix == "neg":
        cmap_name = "Blues"
        color = "blue"
        reverse = False
    else:
        cmap_name = "Reds"
        color = "red"
        reverse = False
    if c == 1:
        cmap_name = "Greys"
        color = "black"
        reverse = False

    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 +
                cell_df[
                    "IL-2_surf_c"] ** 3).values

    alphas = np.linspace(0.1, 1, len(fractions) + 1)
    for f, f_Tsec in enumerate(fractions):
        frac_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == f_Tsec) & (cell_df["type_name"] == "Th")]

        frac_df["time"] /= (3600)

        sns.set_style("ticks")

        rep_act = []
        for rep in frac_df["replicat_index"].unique():
            act = []
            for t, time in enumerate(frac_df["time"].unique()):
                timed_df = frac_df.loc[(frac_df["replicat_index"] == rep) & (frac_df["time"] == time)]
                act.append(len(timed_df.loc[(timed_df["pSTAT5"] > 0.5)]) / len(timed_df) * 100)
            rep_act.append(act)
            # act = np.array(act) * 100
        sns.lineplot(x = frac_df["time"].unique(), y = np.mean(rep_act, axis = 0), label=f_Tsec, color=color, alpha=alphas[f], legend=False)
        # my_interpolation(x = cell_df["IL-2_gamma"].unique(), y = np.mean(rep_act, axis = 0), smoothing_factor = 0.01,
        #                  plot=True, label=f_Tsec, color=palette[f], fill = True, std = np.std(rep_act, axis=0), fill_smoothing_factor = 0.05)
        plt.fill_between(frac_df["time"].unique(), np.mean(rep_act, axis = 0) - np.std(rep_act, axis = 0),
                         np.mean(rep_act, axis = 0) + np.std(rep_act, axis = 0), alpha=0.3 * alphas[f], color=color, linewidth=0)

        if plot_every_cell == True:
            plt.ylabel("pSTAT5")
        else:
            plt.ylabel("%pSTAT5$^+$")
        plt.xlabel("time (d)")
        plt.yscale(yscale)
        plt.xscale(xscale)
        plt.ylim(ylim)
        plt.yticks([0, 50, 100])
        if prefix == "pos":
            plt.xticks([0, 48, 96], [0,2,4])
            plt.xlim(xlim)
        else:
            plt.xticks([0, 24 * 4, 24 * 8], [0,4,8])
            plt.xlim((0, 24 * 8))

if save_plot == True:
    fig.savefig(saving_string + filename + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
