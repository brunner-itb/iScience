import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
rc_ticks["figure.figsize"] = (1.67475, 1.386)
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df
from thesis.cellBehaviourUtilities.halftime_estimation import halftime_estimation
import pandas as pd


def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

saving_string = r"/home/brunner/Documents/Current work/2024_09_06/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)

# feedback_type = "pos"
feedback_type = "neg"


gamma = None
scan_index = None
if feedback_type == "pos":
    # model_name = "Tsec_scan_6"
    # name = "dataframes_positive_timeseries" #michaelis
    # gamma = 12.
    # scan_index = 6
    model_name = "feedback_scan_5"
    name = "dataframes_positive_0.176_time_till_4_days_timeseries" #michaelis
    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
    gamma = 12.
    fraction = 0
    filename = "pos_surf_c_over_time_every_cell_and_means.pdf"
    x_ticks = [0, 2, 4]
    xlim = (None, None)
elif feedback_type == "neg":
    # model_name = "Tsec_scan_6"
    # name = "dataframes_negative_timeseries" # menten
    # gamma = 0.01
    # scan_index = 5
    model_name = "feedback_scan_4"
    name = "dataframes_negative_timeseries"  # compute4
    gamma = 0.01
    fraction = 3
    path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
    filename = "neg_surf_c_over_time_every_cell_and_means.pdf"
    x_ticks = [0, 4, 8]
    xlim = (None, 8)


replicat_index = 0

yscale = "log"
xscale = "linear"

# plotting limits
ylim = (0, 1)

# if line smoothing is desired, sfs are the smoothing factors
c_interpolate = False
c_sf = 0

try:
    cell_df, global_df =  my_load_df(path, offset=0, custom_ending = "_combined")
except FileNotFoundError:
    cell_df = pd.read_hdf(path + "cell_df" + ".h5", mode="r")
try:
    cell_df["id_id"]
except KeyError:
    cell_df["id_id"] = cell_df["id"]

cell_df["IL-2_surf_c"] *= 1e3
# assert 1 == 0
# plot only a random fraction of the cells due to visibility
picked_fraction = 1/10

# define which x_axis to plot
x_axis = "time"

# the first time index does not have pSTAT5 calculated yet.
skip_first_time_index = True

plot_legend = False

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################

#%%
print("plotting")
for plot_df in [cell_df]:
    rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()


    if gamma in plot_df["misc_gamma"].unique():
        plot_df = plot_df.loc[(plot_df["misc_gamma"] == gamma)]
    elif scan_index != None:
        plot_df = plot_df.loc[(plot_df["misc_gamma"] == np.sort(plot_df["misc_gamma"].unique())[scan_index])]
        # plot_df = plot_df.loc[(plot_df["scan_name_scan_name"] == 'fsec_scangamma_100_ftreg_0_fth_0.99')]
    plot_df = plot_df.loc[
        (plot_df["fractions_Tsec"] == np.sort(plot_df["fractions_Tsec"].unique())[fraction])]

    if gamma in plot_df["misc_gamma"].unique():
        plot_df = plot_df.loc[(plot_df["misc_gamma"] == gamma)]
    elif gamma != None:
        pass

    plot_df = plot_df.loc[plot_df["replicat_index"] == replicat_index]
    plot_df = plot_df.loc[plot_df["type_name"] == "Th"]

    if skip_first_time_index == True:
        plot_df = plot_df.loc[plot_df["time_index"] > 0]

    # if plot_df["IL-2_gamma"].unique()[0] >= 1:
    #     ht_tilde = halftime_estimation(plot_df["misc_eta"].unique()[0], plot_df["misc_pos_half_time"].unique()[0], plot_df["misc_R_start_pos"].unique()[0])
    plot_df["time"] /= (3600 * 24)
    plot_df["time"] -= plot_df["time"].min()

    try:
        plot_df["pSTAT5"] = plot_df["IL-2_pSTAT5"]
    except:
        pass
        # plot_df["pSTAT5"] = plot_df["IL-2_surf_c"] ** 3 / (
        #         (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=plot_df["IL-2_R"]) * 1e12) ** 3 + plot_df[
        #     "IL-2_surf_c"] ** 3).values

    plot_df["activated"] = "pSTAT5$^-$"
    act_ids = plot_df.loc[(plot_df["time"] == plot_df["time"].max()) & (plot_df["pSTAT5"] > 0.5) & (plot_df["type_name"] == "Th"), "id"]
    plot_df.loc[plot_df["id"].isin(act_ids), "activated"] = "pSTAT5$^+$"

    cell_ids = plot_df["id_id"].unique()
    np.random.seed(1)
    picked_cell_ids = np.random.choice(cell_ids, int(len(cell_ids) * picked_fraction), replace=False)
    # if feedback_type == "pos":
    #     picked_cell_ids = np.array([229, 267, 326, 369, 291, 994, 135, 995, 594, 990, 575, 429, 379,
    #         104, 584,  44, 152, 796, 540, 117])
    # elif feedback_type == "neg":
    #     picked_cell_ids = np.array([601, 345, 448, 613, 774, 380, 734, 661,  95, 637, 527, 303, 348,
    #        434, 226, 364, 350, 702, 103, 172, 573, 942, 793])
    reduced_df = plot_df.loc[plot_df["id_id"].isin(picked_cell_ids)]

    sns.set(rc={"lines.linewidth": 0.2})
    sns.lineplot(x="time", y="IL-2_surf_c", data=reduced_df.loc[(reduced_df["type_name"] == "Th") & (~reduced_df["id"].isin(act_ids))].sort_values(by="time"), estimator=None,
                        units="id_id", color="grey", alpha=1, zorder=-10)
    sns.lineplot(x="time", y="IL-2_surf_c", data=reduced_df.loc[(reduced_df["type_name"] == "Th") & (reduced_df["id"].isin(act_ids))].sort_values(by="time"), estimator=None,
                        units="id_id", color="red" if feedback_type == "pos" else "blue", alpha=0.4, zorder=-10)

    sns.set(rc={"lines.linewidth": 0.8})
    palette = ["red" if feedback_type == "pos" else "blue", "black"]
    sns.lineplot(x="time", y="IL-2_surf_c",
                        data=plot_df.loc[(plot_df["type_name"] == "Th")].sort_values(by="time"),
                      palette=palette, ci=0, hue="activated", legend=False, hue_order=["pSTAT5$^+$", "pSTAT5$^-$"], alpha=0.5)
    plt.xlabel("time (d)")
    plt.yscale("log")
    plt.xscale("linear")
    plt.ylim((1e-1, 1e3))
    plt.xlim(xlim)
    plt.xticks(x_ticks)
    # plt.yticks([0, 20, 40])
    if feedback_type == "pos":
        plt.ylabel(r"surface conc. avg. (pM)")
        plt.yticks([1e-1, 1e0, 1e1, 1e2, 1e3])
        plt.xlim(0,4)
    else:
        plt.yticks([1e-1, 1e0, 1e1, 1e2, 1e3], ["", "", "", "", ""])
        plt.ylabel("")
        plt.xlim(0,8)
    ax.set_rasterization_zorder(0)
    if save_plot == True:
        fig.savefig(saving_string + "means_" + filename, bbox_inches='tight', transparent=True)
    plt.tight_layout()
    plt.show()