import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

from thesis.scripts.paper_models.utilities.plot_helper import my_load_df, my_interpolation

'''
define which static sim you want to plot. 
bc defines the boundary condition, scan variable the scan.
"IL-2_Tsec_fraction" is the secreting cells scan, "IL-2_sigma" the R_lognorm scan and IL-2_KD the KD scan.
'''

# bc = "standard"
bc = "saturated"

saving_string = r"/home/brunner/Documents/Current work/2023_11_03/boxplot_static_surf_c_homogen.pdf"

cell_types = ["Th"]

scan_measure = "IL-2_surf_c"

replicat = "replicat_index"

# if line smoothing is desired, sfs are the smoothing factors
interpolate = True
sfs = [0.01, 0.005, 0.005]

########################################################################################################################
########################################################################################################################
##                                              Plotting                                                              ##
########################################################################################################################
########################################################################################################################


print("loading data")


homo_path = "/extra2/brunner/paper_models/statics/homogeneous_secretion/comparison_to_tsec_scan_halfed_2/"
homo_c_df, homo_g_df = my_load_df(homo_path, offset=0, custom_ending="_combined")

ODE_path = '/extra2/brunner/paper_models/ODE/saturated/activation_q_10_R_1e4_q_scan_1/'
# ODE_path = '/extra2/brunner/paper_models/ODE/saturated/Tsec_scan_fixed_2/'
ODE_c_df, ODE_g_df = my_load_df(ODE_path, offset = 0, run_range = [0], custom_ending = "_combined")
ODE_c_df["IL-2_surf_c"] *= 1e3
########################################################################################################################
print("plotting")
########################################################################################################################

# fracs = [0.051, 0.076, 0.126, 0.2, 0.3, 0.4]
fracs = ODE_c_df["IL-2_q"].unique()
ODE_boxplot_data = []
for frac in fracs:
    ODE_frac_df = ODE_c_df.loc[(ODE_c_df["IL-2_q"] == frac)]
    ODE_frac_df = ODE_frac_df.loc[ODE_frac_df["type_name"] == "Th"]

    act_list = []
    R_list = []
    for idx, rep in enumerate(np.sort(ODE_frac_df[replicat].unique())):
        tmp_df = ODE_frac_df.loc[(ODE_frac_df[replicat] == rep)]

        frac_act = tmp_df["IL-2_surf_c"].mean()

        act_list.append(frac_act)
    if len(act_list) == 0:
        print(frac)
    ODE_boxplot_data.append(np.array(act_list))
#%%
svs = homo_c_df["scan_value"].unique()
homo_boxplot_data = []
for q in np.sort(svs):
    # frac_df = big_c_df.loc[(np.abs(big_c_df["IL-2_Tsec_fraction"] - frac) < 1e-3) & (big_c_df["time_index"] == 12.)]
    frac_df = homo_c_df.loc[(np.abs(homo_c_df["scan_value"] - q) < 1e-5) & (homo_c_df["time_index"] == 0.)]
    print(frac_df.loc[frac_df.type_name == "Th", "IL-2_R"].mean())
    act_list_2 = []
    for rep in frac_df["replicat_index"].unique():
        rep_df = frac_df.loc[(frac_df["replicat_index"] == rep) & (frac_df.type_name == "Th")]
        rep_df["IL-2_surf_c"] *= 1e3
        act_list_2.append(rep_df["IL-2_surf_c"].mean())
        # act_list_2.append(rep_df["IL-2_surf_c"].std()/rep_df["IL-2_surf_c"].mean())
    homo_boxplot_data.append(np.array(act_list_2))

np.random.seed(1)
# rc_ticks["figure.figsize"] = (0.95 * 1.25, 1.19)
# rc_ticks['figure.figsize'] = (1.67475 * 1.32, 1.386 * 1.27)
# rc_ticks["xtick.labelsize"] = 7
# rc_ticks["ytick.labelsize"] = 6
# rc_ticks['axes.labelsize'] = 7
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

linewidth = 0.5
alpha = 0.9

color1 = "#aa0000ff"
color2 = "grey"


flierprops = dict(marker='o', markerfacecolor=color1, markersize=0,
                  linestyle='none', markeredgecolor=color1, linewidth=linewidth)
boxprops=dict(linewidth=linewidth)
whiskerprops=dict(linewidth=linewidth)
capsprops=dict(linewidth=linewidth)


# ODE_boxp = plt.boxplot(ODE_boxplot_data, flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops, widths=0.5, patch_artist=True)

# two_colours = [[color1] * len(svs)]
# from matplotlib.colors import to_rgba
# rgba_1 = list(to_rgba(color1))
# rgba_1[-1] = 0.6
# rgb_colors = [[rgba_1] * len(svs)]
# four_colours = [[color1] * len(svs) * 2]

two_colours = [[color1] * len(fracs), [color2] * len(fracs)]
from matplotlib.colors import to_rgba
rgba_1 = list(to_rgba(color1))
rgba_1[-1] = 0.3
rgba_2 = list(to_rgba(color2))
rgba_2[-1] = 0.3
rgb_colors = [[rgba_1] * len(fracs),
              [rgba_2] * len(fracs)]
four_colours = [[color1] * len(fracs) * 2, [color2] * len(fracs) * 2]

for b, bdata in enumerate([homo_boxplot_data]):
    x = homo_c_df["scan_value"].unique() * (260**3)/100
    # alphas = np.logspace(0, 0, len(bdata)) if b == 1 else np.linspace(1,1, len(bdata))
    alphas = np.linspace(1,1, len(bdata))
    whisker_alphas = np.array([[a] * 2 for a in alphas]).flatten()
    boxp = plt.boxplot(bdata, positions=x, flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops, widths=4, patch_artist=True)

    for element in ["fliers", "means", "medians"]:
        for e,entry in enumerate(boxp[element]):
            plt.setp(entry, color=two_colours[b][e], alpha=alphas[e]/2)

    for element in ['boxes']:
        for e,entry in enumerate(boxp[element]):
            plt.setp(entry, facecolor=rgb_colors[b][e], edgecolor=two_colours[b][e], alpha=alphas[e]/2)

    for element in ['whiskers', "caps"]:
        for e,entry in enumerate(boxp[element]):
            plt.setp(entry, color=four_colours[b][e], linewidth=linewidth, alpha=np.clip(whisker_alphas[e] * 2,0,1))
    for e,entry in enumerate(bdata):
        factor = 1/1
        subsample = np.random.choice(entry, int(len(entry) * factor), replace=False)
        plt.plot(x[e] + np.random.normal(0, 1, size=len(subsample)), subsample, "o", color = two_colours[b][e], alpha=alphas[e], markersize=0.5)

plt.plot(ODE_c_df["IL-2_q"].unique(), np.mean(ODE_boxplot_data, axis=1), color=color2, markersize=0.5, linewidth=0.6)
# import matplotlib.lines as mlines
# handles = [mlines.Line2D([], [], color=color1, linewidth=2), mlines.Line2D([], [], color=color2, linewidth=2)]
# labels = ["well mixed", "spatial resolution"]
# plt.legend(handles=handles, labels=labels)


# plt.ylabel("avg. surface conc. (pM)")
plt.ylabel(r"surface conc. avg. (pM)")
# plt.xlabel("%pSTAT$^+$")
plt.xscale("linear")
plt.yscale("log")
# plt.xlim((0.65, 2.35))
# plt.ylim(0, 26)
# plt.ylim(0, 315)
# plt.ylim(-4, 104)
# plt.ylim(-0.5, 15.5)

# ax.set_xticklabels(["well\nmixed", "spatial\nresolution"])
ax.set_xticks([int(x) for x in homo_c_df["scan_value"].unique() * (260**3)/100][::2])
ax.set_xticklabels([int(x) for x in homo_c_df["scan_value"].unique() * (260**3)/100][::2])
# ax.set_xticks([1, 2, 4, 6, 8, 10])
# all_ticks = [int(x) for x in homo_c_df["scan_value"].unique() * (260**3)/100]
# ax.set_xticklabels([all_ticks[1], all_ticks[2], all_ticks[4], all_ticks[6], all_ticks[8], all_ticks[10]])
# ax.set_xticks([1, 2, 3, 4])
# ax.set_xticklabels(["1", "5", "10", "15"])
plt.xlabel(r"secretion rate (100/s)")
# ax.set_xticks([])
# plt.yticks([0,12.5,25])
# plt.xlim((None, 10.5))
plt.title("homogeneous secretion")
plt.ylim(1e-1, None)

fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()