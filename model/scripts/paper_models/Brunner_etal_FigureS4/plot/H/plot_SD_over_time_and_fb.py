import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

# prefix = "test"
# prefix = "pos"
prefix = "neg"
# prefix = "ODE_pos"
if prefix == "pos":
    # model_name = "Tsec_scan_5"
    # name = "gamma_12" #compute2 & compute3
    # fraction = 7
    # model_name = "feedback_scan_for_fig3_D"
    model_name = "feedback_scan_5"
    name = "dataframes_positive_0.176_time_till_4_days_timeseries" #michaelis
    fraction = 0
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/gamma_12/"
elif prefix == "neg":
    model_name = "feedback_scan_4"
    name = "dataframes_negative_timeseries"  # compute4
    fraction = 4
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/Tsec_scan_5/gamma_0.01/"


saving_string =r"/home/brunner/Documents/Current work/2024_09_06/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

filename = prefix + "_SD_over_time_and_fb"

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


# the first time index does not have pSTAT calculated yet.
skip_first_time_index = True
########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

base_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, extra=hdd, n=name)

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
# print(spatial_cell_df["time"]/(3600 * 24))
# print(spatial_cell_df["IL-2_Tsec_fraction"].unique())
#%%
print("plotting")
f_cell_df = spatial_cell_df.loc[(spatial_cell_df["IL-2_Tsec_fraction"] == np.sort(spatial_cell_df["IL-2_Tsec_fraction"].unique())[fraction])]
f_ODE_cell_df = ODE_cell_df.loc[(ODE_cell_df["IL-2_Tsec_fraction"] == (0.1 if prefix == "pos" else 0.102))]

rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
# rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 0.5)
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
# amount_of_reps = 10
for c,cell_df in enumerate([f_cell_df]):
    gammas = cell_df["IL-2_gamma"].unique()
    gammas = np.sort(gammas[np.where(gammas != 1.)])

    if prefix == "neg":
        cmap_name = "Blues"
        color = "blue"
        gammas = gammas[::-1] # reversed
    else:
        cmap_name = "Reds"
        color = "red"
        reverse = False
        gamma_range = np.around(np.logspace(-2, np.log10(12), 30), 3)
        gamma_range = gamma_range[gamma_range > 1]
        gammas = gamma_range
    if c == 1:
        cmap_name = "Greys"
        color = "black"
        reverse = False
        gammas = cell_df["IL-2_gamma"].unique()

    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 +
                cell_df[
                    "IL-2_surf_c"] ** 3).values

    alphas = np.logspace(-1, 0, len(gammas) + 1) if c != 1 else [1] #1 for ODE
    for g, gamma in enumerate(gammas):
        frac_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
        frac_df["time"] /= 3600

        sns.set_style("ticks")
        rep_SD = []
        for rep in frac_df["replicat_index"].unique():
            SD = []
            rep_df = frac_df.loc[(frac_df["replicat_index"] == rep) & (frac_df["type_name"] == "Th")]
            if len(rep_df["time"].unique()) == len(frac_df["time"].unique()):
                for t, time in enumerate(rep_df["time"].unique()):
                    timed_df = rep_df.loc[(rep_df["time"] == time)]
                    SD.append(timed_df["IL-2_surf_c"].std())
                rep_SD.append(SD)
        sns.lineplot(x = frac_df["time"].unique(), y = np.mean(rep_SD, axis = 0), label=gamma, color=color, alpha=alphas[g], legend=False)
        # my_interpolation(x = cell_df["IL-2_gamma"].unique(), y = np.mean(rep_SD, axis = 0), smoothing_factor = 0.01,
        #                  plot=True, label=f_Tsec, color=palette[f], fill = True, std = np.std(rep_SD, axis=0), fill_smoothing_factor = 0.05)
        plt.fill_between(frac_df["time"].unique(), np.mean(rep_SD, axis = 0) - np.std(rep_SD, axis = 0)/np.sqrt(len(frac_df["replicat_index"].unique())),
                         np.mean(rep_SD, axis = 0) + np.std(rep_SD, axis = 0)/np.sqrt(len(frac_df["replicat_index"].unique())), alpha=0.3 * alphas[g], color=color, linewidth=0)


        # plt.ylabel(r"pSTAT$^+$ T$_{\rm resp}$ cells (%)")
        plt.ylabel(r"surface conc. s.d. (pM)")
        plt.xlabel("time (d)")
        plt.yscale("log")
        plt.xscale(xscale)
        plt.ylim((1e0, 3e2))
        # plt.yticks([0, 50, 100])
        if prefix == "pos":
            plt.xticks([0, 48, 96], [0,2,4])
            plt.xlim((0, 4 * 24))
        else:
            plt.xticks([0, 4 * 24, 8 * 24], [0,4,8])
            plt.xlim((0, 8 * 24))

if save_plot == True:
    fig.savefig(saving_string + filename + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()