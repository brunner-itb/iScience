import getpass
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, EC50_calculation
from scipy.interpolate import UnivariateSpline

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

prefix = "pos"
# prefix = "neg"

if prefix == "pos":
    model_name = "Tsec_scan_6"
    name = "dataframes_positive_2_steady_state" #compute2 & compute3
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/Tsec_scan_5/gamma_100/"
    xlim = (1, 50)
elif prefix == "neg":
    model_name = "Tsec_scan_6"
    # name = "gamma_0.5"  # compute4
    # model_name = "feedback_scan_4"
    name = "dataframes_negative_steady_state/" # michaelis
    ODE_path = f"/{hdd}/brunner/paper_models/ODE/saturated/kinetics/Tsec_scan_5/gamma_0.01/"
    xlim = (1, 20)

# model_name = "feedback_scan_4"
# name = "dataframes_negative_timeseries"  # c4

saving_string = r"/home/brunner/Documents/Current work/2023_10_27/"
if not os.path.exists(saving_string):
    os.mkdir(saving_string)
filename = prefix + "_SD_over_tsecs.pdf"

# choose a scan and replicat index
path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, n=name, mn=model_name, extra=hdd)
spatial_cell_df, global_df =  my_load_df(path, offset=0, custom_ending = "_combined")
spatial_cell_df["IL-2_surf_c"] *= 1e3
try:
    spatial_cell_df["IL-2_gamma"]
except KeyError:
    spatial_cell_df["IL-2_gamma"] = spatial_cell_df["misc_gamma"]
    spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]
# spatial_cell_df = spatial_cell_df[spatial_cell_df["IL-2_gamma"] == 100]
print(spatial_cell_df["IL-2_gamma"].unique())

t0_df = spatial_cell_df.loc[(spatial_cell_df["time_index"] == spatial_cell_df["time_index"].max()) & (spatial_cell_df["IL-2_Tsec_fraction"] == 0.101) &
                                    (spatial_cell_df["IL-2_gamma"] == spatial_cell_df["IL-2_gamma"].unique()[-1])]


#%%
ODE_cell_df =  pd.read_hdf(ODE_path + "cell_df" + ".h5", mode="r")
ODE_cell_df["IL-2_surf_c"] *= 1e3

plot_every_cell = True
plot_std_area = True
# define which scan to plot
# the first time index does not have pSTAT calculated yet.
skip_first_time_index = True
########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################
# no_fb_cell_df = spatial_cell_df.loc[(spatial_cell_df["time_index"] == spatial_cell_df["time_index"].min() + 2) &
#                                     (spatial_cell_df["IL-2_gamma"] == spatial_cell_df["IL-2_gamma"].unique()[0])]

if skip_first_time_index == True:
    spatial_cell_df = spatial_cell_df.loc[spatial_cell_df["time_index"] != 0]

#%%
labels = ["feedback", "ODE"]
plot_x = [[], [], []]
plot_y = [[], [], []]
plot_std = [[], [], []]
loop_var = "scan_name_scan_name"
# for c,cell_df in enumerate([ODE_cell_df]):
for c,cell_df in enumerate([spatial_cell_df]):
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print(f"calculating own activation for c = {c}")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
                        "IL-2_surf_c"] ** 3).values
    cell_df = cell_df.loc[cell_df["type_name"] == "Th"]
    cell_df = cell_df.loc[cell_df["time"] == cell_df["time"].max()]
    for g, gamma in enumerate(np.sort(cell_df["IL-2_gamma"].unique())):
        gamma_df = cell_df.loc[(cell_df["IL-2_gamma"] == gamma)]
        CV = np.zeros((len(gamma_df["IL-2_Tsec_fraction"].unique()), len(gamma_df["replicat_index"].unique()))) # amount_of_components / amount_of_Tsecs
        SD = np.zeros((len(gamma_df["IL-2_Tsec_fraction"].unique()), len(gamma_df["replicat_index"].unique()))) # amount_of_components / amount_of_Tsecs
        CV[:] = np.nan
        SD[:] = np.nan
        print("gamma:", gamma)
        for f, frac in enumerate(np.sort(gamma_df["IL-2_Tsec_fraction"].unique())):
            frac_df = gamma_df.loc[(gamma_df["IL-2_Tsec_fraction"] == frac)]
            for r, rep in enumerate(np.sort(frac_df["replicat_index"].unique())):
                rep_df = frac_df.loc[(frac_df["replicat_index"] == rep)]
                CV[f][r] = np.nanmean(rep_df["IL-2_surf_c"].std()/rep_df["IL-2_surf_c"].mean()) if rep_df["IL-2_surf_c"].mean() != 0 else 0
                SD[f][r] = np.nanmean(rep_df["IL-2_surf_c"].std())
        plot_x[c].append(np.sort(gamma_df["IL-2_Tsec_fraction"].unique()) * 100)
        # plot_y[c].append(np.nanmean(CV, axis=1))
        # plot_std[c].append(np.nanstd(CV, axis=1))
        plot_y[c].append(np.nanmean(SD, axis=1))
        plot_std[c].append(np.nanstd(SD, axis=1))
# print(acts)
#%%
print("plotting")
rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
# rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 0.5)
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
to_be_plotted_y = plot_y
to_be_plotted_std = plot_std
for i in range(len(to_be_plotted_y[:1])): #cell_dfs
    alphas = np.logspace(-1, 0, len(to_be_plotted_y[i]))
    if prefix == "neg":
        alphas = [x for x in reversed(alphas)]
        print("reversed")
        colour = "Blue"
    elif prefix == "pos":
        colour = "Red"
    for e, entry in enumerate(to_be_plotted_y[i]): #gammas
        smoothing_factor = 0.5e-2 if prefix == "pos" else 5e-2
        x = plot_x[i][e]
        std = to_be_plotted_std[i][e]

        # delete weird outliers
        if np.max(entry) > 60:
            index = np.where(entry > 60)[0][0]
            x = np.delete(x, index)
            entry = np.delete(entry, index)
            std = np.delete(std, index)

        # extend the data a bit so boundary effects vanish, only relevant if we smooth via spline
        longer_y = np.concatenate([entry, [entry[-1] for x in range(10)]])
        longer_std = np.concatenate([std, [std[-1] for x in range(10)]])
        longer_x = np.concatenate([x, [x[-1]+ 0.01*bla for bla in range(10)]])

        xnew = np.linspace(longer_x.min(), np.max([sublist[-1] for sublist in plot_x[i]]), 200)
        spl = UnivariateSpline(longer_x, np.log10(longer_y))
        spl.set_smoothing_factor(smoothing_factor)

        std_spl = UnivariateSpline(longer_x, np.log10(longer_std))
        std_spl.set_smoothing_factor(smoothing_factor)

        # plt.plot(xnew, 10 ** (spl(xnew)), label=labels[i], color=colour if i != 1 else "black",
        #          alpha=alphas[e] if i != 1 else 1)
        # plt.fill_between(xnew,
        #              np.clip(10 ** (spl(xnew)) - 10 ** (std_spl(xnew))/np.sqrt(10), 0, None),
        #              np.clip(10 ** (spl(xnew)) + 10 ** (std_spl(xnew))/np.sqrt(10), 0, None), color=colour if i != 1 else "black",
        #                  alpha=alphas[e] * 0.5 if i != 1 else 0.15, linewidth=0.0)
        plt.plot(x, entry, label=labels[i], color=colour if i != 1 else "black", alpha=alphas[e] if i != 1 else 1)
        plt.fill_between(x, entry - std/np.sqrt(10), entry + std/np.sqrt(10), color=colour if i != 1 else "black",
                         alpha=alphas[e] * 0.5 if i != 1 else 0.15, linewidth=0.0)


plt.ylabel("SD (pM)")
plt.xlabel(r"fraction T$_{\rm sec}$ (%)")
plt.yscale("log")
plt.xscale("linear")
plt.ylim((1, 3e2))
plt.xlim(xlim)
# plt.yticks([0, 2, 4])
max_x_ticks = np.max(plt.xticks()[0])
plt.xticks([1, max_x_ticks/2, max_x_ticks])
# plt.title("pSTAT5", pad = title_pad)
# plt.axhline(0.5, 0, 200, color="black")
# plt.legend()


if save_plot == True:
    fig.savefig(saving_string + filename, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
