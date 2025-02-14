import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from scipy.interpolate import UnivariateSpline
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation
from thesis.scripts.paper_models.statics.figure_2.plotting.funcs_for_plotting import get_path, scale_dataframes
'''
define which static sim you want to plot. 
bc defines the boundary condition, scan variable the scan.
"IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan and IL-2_KD the KD scan.
'''
bc = "linear"
# bc = "saturated"
# scan_variable = "IL-2_Tsec_fraction" #michaelis
# scan_variable = "IL-2_sigma" #michaelis
scan_variable = "IL-2_KD" #michaelis

scan_variables = [scan_variable]

models = ["loc. q and R", "well_mixed"] if scan_variable != "IL-2_KD" else ["loc. q and R", "well_mixed", "k_on_linear"]
# models = ["loc. q and R"]

saving_string =r"/home/brunner/Documents/Current work/2023_11_03/" + "{m}_static_{sv}_{bc}".format(m = "spatial", sv=scan_variable, bc = bc) + ".pdf"

# plotting parameters
plot_legend = False

yscale = "linear"
xscale = "linear"
# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum value
max_value = 1e5
min_value = -1
# x and ylims have to be adjusted for each plot
if scan_variable == "IL-2_Tsec_fraction":
    xlim = (0.0, 0.4)
elif scan_variable == "IL-2_sigma":
    xlim = (0.1, 10)
else:
    xlim = (3, 20)

ylim = (0, 16)
yticks = [0, 8, 16]

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0
# which cell type to plot
cell_types = ["Th"]
# Which colors to use, dim-2 first entry is spatial, second is ODE. dim-1 first entry is concentration, second SD
# y_colours = [["#aa0000ff", "darkorange"] , ["black"], ["blue"]]
y_colours = [["black", "blue"] , ["black"], ["green"]]
legend_labels = ["RD-system", "well-mixed", "s.d."]
linestyles = [["-", "-"], ["--"], ["-"]]

y_variable = "surf_c"
title = "Surface concentration"

replicat = "replicat_index"

# if line smoothing is desired
c_interpolate = 0
c_sf = 1e-2
c_ste_sf = 1e-2 #if scan_variable != "IL-2_sigma" else 1e-3
std_interpolate = 0
std_sf = 5e-4 if scan_variable != "IL-2_Tsec_fraction" else 3e-3
std_ste_sf = 5e-3

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################
# enter the paths for all the scans here
dataframes = [[] for x in models]
print("loading data")
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
for m,model in enumerate(models):
    for sv, scan_variable in enumerate(scan_variables):
        m_sv_path = get_path(bc, hdd, model, scan_variable)
        print(m_sv_path)
        c_df, g_df = my_load_df(m_sv_path, offset = 0, run_range = [0], custom_ending = "_combined")
        if scan_variable == "IL-2_Tsec_fraction":
            try:
                c_df["IL-2_Tsec_fraction"]
            except KeyError:
                c_df["IL-2_Tsec_fraction"] = c_df["fractions_Tsec"]
                g_df["IL-2_Tsec_fraction"] = g_df["fractions_Tsec"]
        if "ODE" not in m_sv_path:
            try:
                c_df["IL-2_sigma"]
            except KeyError:
                c_df["IL-2_sigma"] = c_df["misc_sigma"]
                try:
                    g_df["IL-2_sigma"] = g_df["scan_value"]
                except:
                    pass
            if model == "k_on_linear":
                c_df["IL-2_KD"] = c_df["IL-2_k_off"].unique()[0]/c_df["IL-2_k_on"] * 1e3
                g_df["IL-2_KD"] = c_df["IL-2_k_off"].unique()[0]/g_df["IL-2_k_on"] * 1e3

        # scales the dataframes. Depends on your desired units. Currently everything is in pM
        c_df, g_df = scale_dataframes(c_df, g_df, model, scan_variable, sv, min_value, max_value, [1])
        dataframes[m].append([c_df, g_df])

test_df = dataframes[0][0][0]
# print(test_df["misc_KD"].unique())
print("KD:", test_df["IL-2_KD"].unique())
print("R:", test_df.loc[test_df.type_name == "Th", "IL-2_R"].mean())
print("q:", test_df["IL-2_q"].unique())
# print("sigma:", test_df["misc_sigma"].unique())
print(test_df["IL-2_sigma"].unique())
print("Tsec:", test_df["fractions_Tsec"].unique())
# print(test_df["IL-2_Tsec_fraction"].unique())
print("k_endo:", test_df["IL-2_k_endo"].unique())
# print(test_df["IL-2_Tsec_fraction"].unique())
########################################################################################################################
print("plotting")
########################################################################################################################
#%%
rc_ticks["figure.figsize"] = (1.7, 1.2)
rc_ticks["axes.labelpad"] = 0.
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()
for m,model in enumerate(models[:1]):
    x = np.arange(0, len(scan_variables) * 2)
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]
        if scan_variable == "IL-2_Tsec_fraction":
            factor = 100
        else:
            factor = 1
        for c, cell_type in enumerate(cell_types):
            tmp_df = c_df.loc[c_df["type_name"] == cell_type]
            ste = []
            std = []
            x_axis = []
            for s, scan in enumerate(np.sort(tmp_df[scan_variable].unique())):
                sliced_df = tmp_df.loc[(tmp_df[scan_variable] == scan)]
                if len(sliced_df) != 0:
                    run_std = []
                    run_mean = []
                    run_R = []
                    if len(np.sort(sliced_df[replicat].unique())) > 10:
                        x_axis.append(scan * factor)
                        for r, run in enumerate(np.sort(sliced_df[replicat].unique())):
                            run_std.append(sliced_df.loc[(sliced_df[replicat] == run), "IL-2_surf_c"].std())
                            run_mean.append(sliced_df.loc[(sliced_df[replicat] == run), "IL-2_surf_c"].mean())
                            run_R.append(sliced_df.loc[(sliced_df[replicat] == run), "IL-2_R"].mean())
                        # print(scan, np.mean(run_R))
                        run_std = np.array(run_std)
                        run_mean = np.array(run_mean)
                        run_std = run_std[~np.isnan(run_std)]
                        # run_std = run_std/run_mean
                        ste.append(np.std(run_std) / np.sqrt(len(sliced_df[replicat].unique())))
                        std.append(np.mean(run_std))

            if std_interpolate == True:
                my_interpolation(x_axis, std, std_sf, plot=True, label="surface conc. s.d. (pM)", color=y_colours[m][1], fill=True,
                                 std=np.array(ste), ci=0, linestyle = linestyles[m][1], ax=ax, legend=plot_legend,
                                 fill_sf=None, extend_r_bdry=True)
            else:
                sns.lineplot(x = x_axis, y = std,
                         color=y_colours[m][1], linestyle=linestyles[m][1], ci=0, label="surface conc. s.d. (pM)", legend=False, ax=ax)
                plt.fill_between(x_axis,
                             np.clip(np.array(std) - np.array(ste), 0, None),
                             np.clip(np.array(std) + np.array(ste), 0, None), alpha=0.3, color=y_colours[m][1], linewidth=0)

ax2 = ax.twinx()
for m,model in enumerate(models):
    if model == "well_mixed":
        c_interpolate = False
    else:
        c_interpolate = True
    for sv,scan_variable in enumerate(scan_variables):
        c_df = dataframes[m][sv][0]
        g_df = dataframes[m][sv][1]
        linew = None
        offset = 0
        if model == "well_mixed": #ODE
            linew = 0.5
            offset = -0.3 #to make the well-mixed line visible
        if scan_variable == "IL-2_Tsec_fraction":
            factor = 100
        else:
            factor = 1
        for c, cell_type in enumerate(cell_types):
            tmp_df = c_df.loc[c_df["type_name"] == cell_type]

            def EC50_calculation(E_max, E_min, k, N, R):
                return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)
            tmp_df["pSTAT5"] = tmp_df["IL-2_surf_c"] ** 3 / (
                    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=0.8, R=tmp_df["IL-2_R"]) * 1e12) ** 3 +
                    tmp_df["IL-2_surf_c"] ** 3).values

            runs_ste = []
            runs_mean = []
            runs_act = []
            for s, scan in enumerate(np.sort(tmp_df[scan_variable].unique())):
                scan_df = tmp_df.loc[tmp_df[scan_variable] == scan]
                # std = []
                mean = []
                act = []
                for r, rep in enumerate(scan_df.replicat_index.unique()):
                    rep_df = scan_df.loc[scan_df.replicat_index == rep]
                    mean.append(rep_df["IL-2_" + y_variable].mean())
                    act.append(len(rep_df.loc[(rep_df["pSTAT5"] > 0.5), "IL-2_surf_c"]) / len(rep_df["IL-2_surf_c"]))
                runs_ste.append(np.std(mean))# / np.sqrt(len(scan_df[replicat].unique())))
                runs_mean.append(np.mean(mean) + offset)
                runs_act.append(np.mean(act))

            if c_interpolate == True:
                my_interpolation(np.sort(tmp_df[scan_variable].unique()) * factor, runs_mean, c_sf, plot=True, label="mean",
                                 color=y_colours[m][0], fill=True, std=np.array(runs_ste), ci=0, linestyle = linestyles[m][0],
                                 legend=plot_legend, ax=ax2, linewidth=linew, fill_sf = c_ste_sf)
            else:
                sns.lineplot(x = np.sort(tmp_df[scan_variable].unique()) * factor, y = runs_mean,
                         color=y_colours[m][0], linestyle = linestyles[m][0], ci=0, label="mean", legend=plot_legend, linewidth=linew, ax=ax2)

                plt.fill_between(np.sort(tmp_df[scan_variable].unique()) * factor,
                             np.clip(np.array(runs_mean) - np.array(runs_ste), 0, None),
                             np.clip(np.array(runs_mean) + np.array(runs_ste), 0, None), alpha=0.3, color=y_colours[m][0], linewidth=0)

if scan_variable == "IL-2_Tsec_fraction":
    xlabel = r"secreting cells (%)"
    standard = 0.05 * 100
elif scan_variable == "IL-2_sigma":
    xlabel = "receptor heterogeneity"
    standard = 1
elif scan_variable == "IL-2_KD":
    xlabel = "saturation constant (pM)"
    standard = 7.4

plt.axvline(standard, color="red")

ax.set(ylabel="surface conc. s.d. (pM)", ylim=(0,16), yticks=yticks)
ax.tick_params(axis='y', colors=y_colours[0][1])
ax.yaxis.label.set_color(y_colours[0][1])
ax.set_xlabel(xlabel)

ax2.set(yscale=yscale, xscale=xscale, ylim=ylim, xlim=xlim)
ax2.set_ylabel(r"surface conc. avg. (pM)", labelpad=0.5)
ax2.set_yticks(yticks)
ax2.tick_params(axis='y', colors=y_colours[0][0])
ax2.yaxis.label.set_color(y_colours[0][0])

if scan_variable == "IL-2_sigma":
    ax.set_xticks([0.1, 2, 4])
    ax.set_xticklabels([0.1, 2, 4])
    ax.set_xscale("log")
elif scan_variable == "IL-2_KD":
    # ax.set_xticks([0.001, 0.2, 0.4], ["", "0.2", "0.4"])
    ax.set_xticks([3, 10, 20], ["3", "10", "20"])
    # ax.set_xticklabels([0, 5, 10])
    # ax.axvline(14.8, color="yellow")
elif scan_variable == "IL-2_Tsec_fraction":
    # ax.set_xticks([0.001, 0.2, 0.4], [r"10$^{-2}$", "0.2", "0.4"])
    # ax.set_xticks([0.01, 0.1, 0.4], ["1", "10", "20"])
    ax.set_xlim([0.5, 40])
    ax2.set_xlim([0.5, 40])
    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 5, 10, 40])
    ax.set_xticklabels([0.5, "",  5, "",  40])
# if bc == "linear":
#     ax.set_title("linear uptake")


fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()