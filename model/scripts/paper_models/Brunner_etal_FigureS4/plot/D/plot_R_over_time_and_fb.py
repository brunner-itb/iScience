import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
rc_ticks["figure.figsize"] = (1.67475, 1.386)
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df
from thesis.main.my_debug import message

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

# prefix = "test"
# prefix = "neg"
prefix = "pos"
# prefix = "ODE_pos"
if prefix == "pos":
    # model_name = "Tsec_scan_5"
    # name = "gamma_12" #compute2 & compute3
    # fraction = 7
    # model_name = "feedback_scan_for_fig3_D"
    # model_name = "feedback_scan_5"
    # name = "dataframes_positive_0.176_time_till_4_days_timeseries" #michaelis
    model_name = "feedback_scan_4"
    name = "dataframes_positive_for_Fig3C_act_plot_timeseries/" #menten
    fraction = 0
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/gamma_10_v2/"
elif prefix == "neg":
    model_name = "feedback_scan_4"
    name = "dataframes_negative_timeseries"  # compute4
    fraction = 4
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/gamma_0.1_v2/"

saving_string =r"/home/brunner/Documents/Current work/2025/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

filename = prefix + "_R_over_time_and_fb"

yscale = "linear"
xscale = "linear"

# plotting limits
xlim = (None,None)
ylim = (-1, 101)

c_interpolate = False
c_sf = 0


plot_every_cell = False
plot_std_area = True
# define which scan to plot
x_axis = "time"

# the first time index does not have pSTAT calculated yet.
skip_first_time_index = True

plot_legend = True

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

ODE_cell_df =  pd.read_hdf(ODE_path + "cell_df" + ".h5", mode="r")
ODE_cell_df.time /= 3600

# print(spatial_cell_df["time"]/(3600 * 24))
# print(spatial_cell_df["IL-2_Tsec_fraction"].unique())
f_cell_df = spatial_cell_df.loc[(spatial_cell_df["IL-2_Tsec_fraction"] == np.sort(spatial_cell_df["IL-2_Tsec_fraction"].unique())[fraction])]
#%%
print("plotting")
# rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 0.65)
rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
rc_ticks["xtick.labelsize"] = 6
rc_ticks["ytick.labelsize"] = 6
rc_ticks["axes.labelsize"] = 7
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
# for c,cell_df in enumerate([f_cell_df, f_ODE_cell_df]):
for c,cell_df in enumerate([f_cell_df]):
    gammas = cell_df["IL-2_gamma"].unique()
    gammas = np.sort(gammas[np.where(gammas != 1.)])

    if prefix == "neg":
        cmap_name = "Blues"
        color = "blue"
        gammas = gammas[::-1] # reversed
        gammas = np.delete(gammas, np.argwhere(gammas == 0.5)) #bugged
    else:
        cmap_name = "Reds"
        color = "red"
        reverse = False
        # gamma_range = np.around(np.logspace(-2, np.log10(12), 30), 3)
        # gamma_range = gamma_range[gamma_range > 1]
        # gammas = gamma_range
        gammas = cell_df["IL-2_gamma"].unique()
        gammas = np.sort(gammas[np.where(gammas != 1.)])
        gamma_range = gammas[gammas < 15]
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
        frac_df = frac_df.loc[(frac_df["type_name"] == "Th")]
        frac_df["time"] /= 3600
        if frac_df.loc[frac_df["time_index"] == 1, "IL-2_R"].mean() > 1700:
            print(gamma)
        Rs = []
        for time in frac_df.time.unique():
            Rs.append(frac_df.loc[frac_df.time == time, "IL-2_R"].mean())
        plt.plot(frac_df.time.unique(), Rs, alpha=alphas[g], color=color, linewidth=0.9) # right plot now. SD is so small anyways, didnt bother to implement it here
        # sns.lineplot(data=frac_df,x="time",y="IL-2_R", estimator='mean', sort=False)
        sns.set_style("ticks")
        rep_act = []
        # sns.lineplot(data = frac_df, x="time", y="EC50_EC50", label=gamma, color=color, alpha=alphas[g],
        #              legend=False, errorbar="se", err_kws={"linewidth":0, "alpha": alphas[g]*0.5}) # Produced the wrong plot somehow!
        if plot_every_cell == True:
            plt.ylabel("pSTAT5")
        else:
            # plt.ylabel(r"pSTAT$^+$ T$_{\rm resp}$ cells (%)")
            plt.ylabel(r"EC50_EC50")
        plt.xlabel("time (d)")
        plt.yscale("log")
        plt.xscale(xscale)
        # plt.yticks([0, 50, 100])
        if prefix == "pos":
            plt.xticks([0, 48, 96], [0,2,4])
            plt.xlim((0, 4 * 24))
        else:
            #
            plt.xticks([0, 4 * 24, 8 * 24], [0,4,8])
            plt.xlim((0, 8 * 24))
        # plt.ylim(6e2, 1e5)
        # plt.ylim(1e3, 2e3)
# sns.lineplot(data=ODE_cell_df, x = "time", y = "IL-2_R", color="black", linewidth=1, errorbar="se", err_kws={"linewidth":0})
# ODE_Rs = []
# for time in ODE_cell_df.time.unique():
#     ODE_Rs.append(ODE_cell_df.loc[ODE_cell_df.time == time, "IL-2_R"].mean())
# plt.plot(ODE_cell_df.time.unique(), ODE_Rs, color="black", linewidth=0.9)
if save_plot == True:
    fig.savefig(saving_string + filename + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
