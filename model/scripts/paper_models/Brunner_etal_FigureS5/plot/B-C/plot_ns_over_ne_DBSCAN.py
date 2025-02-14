import getpass
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
import os
from scipy.optimize import curve_fit
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, EC50_calculation, my_interpolation, myDBscan


save_plot = True
hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()

# prefix = "pos"
prefix = "neg"

save_plot = True
niche_effect_variable = "pSTAT5"

yscale = "log"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# time and fraction max limit
max_time = 100
max_frac = -1
# colours, first = positive, second = negative
colours = ["black", "blue"]
# opacity for the ODE
ODE_alpha = 0.25
# if interpolation is desired
interpolate = False
sf = 0.003
if prefix == "neg":
    model_name = "feedback_scan" #c4
    # name = "test_q_100_2"
    # name = "dataframes_negative_4_steady_state"
    name = "dataframes_negative_for_FigS5C_2_steady_state"  # michaelis
    # my_hue = "IL-2_gamma"
    my_hue = "scan_name_scan_name"
    Tsec_fraction = 0
elif prefix == "pos":
    model_name = "feedback_scan_4" #c4
    name = "dataframes_positive_10_steady_state"
    Tsec_fraction = 2
    my_hue = "scan_name_scan_name"
    # model_name = "Tsec_scan_6"
    # # name = "dataframes_positive_timeseries" #compute2 & compute3
    # name = "dataframes_positive_steady_state" #compute2 & compute3
    # Tsec_fraction = 3
    # my_hue = "misc_gamma"



saving_string =f"/home/brunner/Documents/Current work/2023_11_03/{prefix}_ns_over_ne_{model_name}_{name}_with_Tsecs_pSTAT5"

pos_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{n}/".format(u=user, mn=model_name, n=name, extra=hdd)
custom_ending = "_combined_2"
try:
    spatial_cell_df = pd.read_hdf(pos_path + "cell_df" + custom_ending + ".h5", mode="r")
    print("read h5")
except:
    print("reading .h5 failed, attempting .pkl")
    spatial_cell_df = pd.read_pickle(pos_path + 'cell_df' + custom_ending + '.pkl')
    print("read pkl")
spatial_cell_df["IL-2_Tsec_fraction"] = spatial_cell_df["fractions_Tsec"]

#%%
print("calculating")
epsilon = 20
# rc_ticks["figure.figsize"] = (8/6 * 1.5, 1.5)
results = []
for c, cell_df in enumerate([spatial_cell_df]):

    cell_df["IL-2_surf_c"] *= 1e3
    try:
        cell_df["pSTAT5"] = cell_df["IL-2_pSTAT5"]
    except:
        print("pSTAT5 loading failed, recalculating")
        cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
                (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
            "IL-2_surf_c"] ** 3).values

    # cell_df["R_fc"] = None
    # cell_df.loc[cell_df["time_index"] == cell_df["time_index"].max(), "R_fc"] = (cell_df.loc[(cell_df["time_index"] == 0), "IL-2_R"].values - cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max()), "IL-2_R"].values)/cell_df.loc[(cell_df["time_index"] == 0), "IL-2_R"].values
    # t0_c_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].min() + 2)]
    if my_hue != "time":
        cell_df = cell_df.loc[(cell_df["time_index"] == cell_df["time_index"].max())]
    if Tsec_fraction != None:
        # cell_df = cell_df.loc[(cell_df["fractions_Tsec"] == np.sort(cell_df["fractions_Tsec"].unique())[Tsec_fraction])]
        cell_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == np.sort(cell_df["IL-2_Tsec_fraction"].unique())[Tsec_fraction])]

    # cell_df.loc[cell_df["R_fc"] < 0.01, "pSTAT5"] = 1
    silhouette_scores = np.zeros((len(cell_df["replicat_index"].unique()), len(cell_df[my_hue].unique())))
    calinski_harabasz_scores = np.zeros((len(cell_df["replicat_index"].unique()), len(cell_df[my_hue].unique())))
    davies_bouldin_scores = np.zeros((len(cell_df["replicat_index"].unique()), len(cell_df[my_hue].unique())))

    # frac_of_sec_cells = np.zeros((len(cell_df["replicat_index"].unique()), len(cell_df[my_hue].unique())))
    niche_concentrations = [[[] for x in cell_df[my_hue].unique()] for y in cell_df["replicat_index"].unique()]
    niche_pSTAT5 = [[[] for x in cell_df[my_hue].unique()] for y in cell_df["replicat_index"].unique()]

    niche_score = np.zeros((len(cell_df["replicat_index"].unique()),
                            len(cell_df[my_hue].unique())))  # amount_of_components / amount_of_Tsecs
    niche_effect = np.zeros((len(cell_df["replicat_index"].unique()),
                            len(cell_df[my_hue].unique())))  # amount_of_components / amount_of_Tsecs
    no_niche_activation = np.zeros((len(cell_df["replicat_index"].unique()),
                            len(cell_df[my_hue].unique())))  # amount_of_components / amount_of_Tsecs

    gammas = np.zeros((len(cell_df["replicat_index"].unique()),
                            len(cell_df[my_hue].unique())))

    for r, rep in enumerate(cell_df["replicat_index"].unique()):
        print(rep)
        rep_df = cell_df.loc[cell_df["replicat_index"] == rep]
        for idx, value in enumerate(np.sort(rep_df[my_hue].unique())[1:]):
            # print(value)
            frac_df = rep_df.loc[(rep_df[my_hue] == value)]
            # coordinates = frac_df[["x", "y", "z"]].values
            # from scipy.spatial import distance_matrix

            # dist_m = distance_matrix(coordinates, coordinates)
            Tsec_ids = frac_df.loc[frac_df["type_name"] == "Tsec", "id"].values

            DBres, sliced_df, DBcoords = myDBscan(frac_df, epsilon, with_Tsecs=True)
            sliced_df["cluster"] = DBres

            if len(np.unique(DBres)) > 1 and len(np.unique(DBres)) < len(DBcoords):
                silhouette_scores[r, idx] = silhouette_score(DBcoords, DBres)
                calinski_harabasz_scores[r, idx] = calinski_harabasz_score(DBcoords, DBres)
                davies_bouldin_scores[r, idx] = davies_bouldin_score(DBcoords, DBres)
            else:
                silhouette_scores[r, idx] = None
                calinski_harabasz_scores[r, idx] = None
                davies_bouldin_scores[r, idx] = None

            niche_score[r][idx] = np.sum(~np.isnan(np.where(np.unique(DBres, return_counts=True)[1] > 1))) / len(
                                    frac_df.loc[frac_df["type_name"] == "Tsec"]) if len(
                                    frac_df.loc[frac_df["type_name"] == "Tsec"]) != 0 else 0
            try:
                gammas[r][idx] = frac_df["misc_gamma"].unique()[frac_df["misc_gamma"].unique() != 1][0]
            except KeyError:
                gammas[r][idx] = frac_df["IL-2_gamma"].unique()[frac_df["IL-2_gamma"].unique() != 1][0]


            effects = []
            for u in np.unique(DBres):
                if len(np.where(DBres == u)[0]) >= 1:
                    ids = sliced_df.loc[sliced_df["cluster"] == u, "id"].values

                    in_cluster_values = sliced_df.loc[sliced_df["cluster"] == u, niche_effect_variable].values
                    outside_cluster_values = frac_df.loc[
                        (frac_df["pSTAT5"] < 0.5) & (frac_df["type_name"] != "Tsec"), niche_effect_variable].values

                    if len(in_cluster_values[~np.isnan(in_cluster_values)]) > 0 and len(outside_cluster_values[~np.isnan(outside_cluster_values)]) > 0:
                        effects.append(np.nanmean(in_cluster_values) / np.nanmean(outside_cluster_values))
                    # effects.append(np.nanmax(in_cluster_concs) / np.nanmax(outside_cluster_concs))
            niche_effect[r][idx] = np.nanmean(effects) if len(effects) > 0 else None

            if niche_score[r][idx] > 0.3 and niche_effect[r][idx] < 300:
                print("value = ", value, "rep = ", rep)

    results.append([niche_score, niche_effect, no_niche_activation, silhouette_scores, calinski_harabasz_scores, davies_bouldin_scores])

# rc_ticks["figure.figsize"] = (8/6 * 1.5, 1.5)
# sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
# fig, ax = plt.subplots()
#%%
print("plotting")
# factor = 1.225
# rc_ticks['figure.figsize'] = [2.8 * factor + 0.15, 2.1 * factor + 0.24]
factor = 0.89
rc_ticks['figure.figsize'] = [1.67475 * factor, 1.46 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
# fig, ax = plt.subplots(figsize = (8/6 * 1.5, 1.5))
for e, entry in enumerate([spatial_cell_df]):
    niche_score = results[e][0]
    niche_effect = results[e][1]

    if prefix == "pos":
        cmap_name = "Red"
        reverse = False
    elif prefix == "neg":
        cmap_name = "Blue"
        reverse = True
    # colors = sns.color_palette(cmap_name + "_r" if reverse == True else cmap_name, len(np.nanmean(niche_score, axis=0)))
    # cmap = plt.cm.get_cmap(cmap_name + "_r" if reverse == True else cmap_name)
    # norm = mpl.colors.SymLogNorm(2, vmin=np.min(dataframes[0][my_hue].unique()), vmax=np.max(dataframes[0][my_hue].unique()))

    no_fb_colors = sns.color_palette("Greys", len(np.nanmean(niche_score, axis=0)))


    import matplotlib as mpl
    for i in range(len(niche_score)):
        # scatterplot = plt.scatter(niche_score[i][np.argsort(gammas[0])] + np.random.normal(0, 0.003, len(niche_score[i])),
        #             niche_effect[i][np.argsort(gammas[0])], c=np.sort(gammas[0]), s = 7, linewidth=0, norm=mpl.colors.LogNorm(),
        #             cmap=cmap_name)
        ax.scatter(niche_score[i][np.argsort(gammas[0])] + np.random.normal(0, 0.002, len(niche_score[i])),
                    niche_effect[i][np.argsort(gammas[0])], c=cmap_name, s = 7, linewidth=0, norm=mpl.colors.LogNorm(),
                    alpha = np.logspace(0,-0.9, len(gammas[0])) if reverse == True else np.logspace(-0.9,0, len(gammas[0])))

    # mean_ns = np.nanmean(niche_score, axis=0)[[np.argsort(gammas[0])]]
    # mean_ne = np.nanmean(niche_effect, axis=0)[[np.argsort(gammas[0])]]

    def func(x,a,b):
        return a + b*x
    no_nans = ~np.isnan(niche_effect.flatten())
    no_zeros = np.where(niche_effect.flatten()[no_nans] != 0)[0]
    popt, pcov = curve_fit(func, niche_score.flatten()[no_nans][no_zeros], niche_effect.flatten()[no_nans][no_zeros])
    if prefix == "neg":
        if Tsec_fraction == 0:
            base = [0.3, func(0.3, *popt)]
            dx, dy = [0.6, func(0.6, *popt)]
        elif Tsec_fraction == 1:
            base = [0.03, func(0.03, *popt)]
            dx, dy = [0.25, func(0.25, *popt)]
        elif Tsec_fraction == 2:
            popt, pcov = curve_fit(func, niche_score.flatten()[no_nans],
                                   niche_effect.flatten()[no_nans])
            base = [0.001, func(0.001, *popt)]
            dx, dy = [0.024, func(0.024, *popt)]

        xoffset = 0
        yoffset = 0
        plt.annotate("", xytext=(dx + xoffset, dy + yoffset), xy=[base[0] + xoffset, base[1] + yoffset],
                     arrowprops={"linewidth": 1.5, "width": 1, "headwidth": 5, "headlength": 5}) #white line around arrow, controlled via linewidth
        plt.annotate("", xytext=(dx + xoffset, dy + yoffset), xy=[base[0] + xoffset, base[1] + yoffset],
                     arrowprops={"linewidth": 0, "width": 1, "headwidth": 5, "headlength": 5, "color": "black"})
    else:
        if Tsec_fraction != 2:
            if Tsec_fraction == 0:
                base = [0.035, func(0.035, *popt)]
                dx, dy = [0.235, func(0.235, *popt)]
            elif Tsec_fraction == 1:
                base = [0.1, func(0.1, *popt)]
                dx, dy = [0.35, func(0.35, *popt)]
            elif Tsec_fraction == 2:
                base = [0.25, func(0.25, *popt)]
                dx, dy = [0.3, func(0.3, *popt)]
            xoffset = 0
            yoffset = 0
            plt.annotate("", xy=(dx + xoffset, dy + yoffset), xytext=[base[0] + xoffset, base[1] + yoffset],
                         arrowprops={"linewidth": 1.5, "width": 1, "headwidth": 5, "headlength": 5}) #white line around arrow, controlled via linewidth
            plt.annotate("", xy=(dx + xoffset, dy + yoffset), xytext=[base[0] + xoffset, base[1] + yoffset],
                         arrowprops={"linewidth": 0, "width": 1, "headwidth": 5, "headlength": 5, "color": "black"})

    # plt.colorbar(scatterplot)
        # plt.scatter(niche_score[i] + np.random.normal(0, 0.002, len(niche_score[i])), niche_effect[i], color="red", alpha=alphas2, s = 0.8, linewidth=0)
        # plt.scatter(niche_score[i], niche_effect[i], color=cmap(norm(dataframes[0][my_hue].unique()[i])), s = 0.8)

# norm = plt.Normalize(np.min(cell_df["IL-2_Tsec_fraction"].unique()), np.max(cell_df["IL-2_Tsec_fraction"].unique()))
# sm = plt.cm.ScalarMappable(cmap=cmap_name + "_r" if reverse == True else cmap_name, norm=norm)
# sm.set_array([])
#
# cbar = fig.colorbar(sm, label="frac. of sec. cells", ticks=[np.min(cell_df["IL-2_Tsec_fraction"].unique()), np.max(cell_df["IL-2_Tsec_fraction"].unique())])
# # cbar.ax.set_yticklabels(["neg.", "pos."])
# cbar.ax.set_yticklabels([round(np.min(cell_df["IL-2_Tsec_fraction"].unique()),5), round(np.max(cell_df["IL-2_Tsec_fraction"].unique()),5)])

plt.xlabel("niche score")
plt.ylabel("niche effect")
plt.yscale("linear")
plt.xscale("linear")
if prefix == "neg":
    if Tsec_fraction == 0:
        plt.xlim((0.15, 0.75))
        plt.ylim((10, 40))
        plt.yticks([10, 25, 40])
    elif Tsec_fraction == 1:
        plt.xlim((0., 0.3))
        plt.ylim((6, 16))
        plt.yticks([6,11, 16])
    elif Tsec_fraction == 2:
        plt.xlim((0, 0.03))
        plt.ylim((0, 8))
        plt.yticks([0, 4, 8])
else:
    if Tsec_fraction == 0:
        plt.xlim((0.02, 0.25))
        plt.ylim((250, 1750))
        plt.yticks([250, 1000, 1750])
        plt.yticks()
    elif Tsec_fraction == 1:
        plt.xlim((0.05, 0.4))
        plt.ylim((0, 500))
        plt.yticks([0,250, 500])
        # plt.xticks([0.05, 0.25, 0.4])
    elif Tsec_fraction == 2:
        plt.xlim((0.14, 0.4))
        plt.ylim((0, 340))
        plt.yticks([0, 170, 340])
        # plt.xticks([0.14, 0.27, 0.4])


# plt.title(np.around(frac_df["IL-2_Tsec_fraction"].unique()[0],2))
# plt.xticks([0, 0.2, 0.4])

if save_plot == True:
    fig.savefig(saving_string + "_" + str(frac_df["IL-2_Tsec_fraction"].unique()[0]) + ".pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
