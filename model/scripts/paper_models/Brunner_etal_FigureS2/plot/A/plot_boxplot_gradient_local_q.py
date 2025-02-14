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

saving_string = r"/home/brunner/Documents/Current work/2024_03_15/boxplot_static_gradient.pdf"

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
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"


# if len(big_c_df) == 0:
import getpass
hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
# base_path = "/extra2/brunner/paper_models/statics/saturated/Figure_1C/act_over_Tsec/"
base_path = "/extra2/brunner/paper_models/kinetics/Figure_1C/dataframes_act_over_Tsec_large_2_steady_state/"
big_c_df, big_g_df = my_load_df(base_path, offset=0, custom_ending="_combined")


try:
    print(np.sort(big_c_df["IL-2_Tsec_fraction"].unique()))
except KeyError:
    big_c_df["IL-2_Tsec_fraction"] = big_c_df["fractions_Tsec"]
    big_g_df["IL-2_Tsec_fraction"] = big_g_df["fractions_Tsec"]
    print(np.sort(big_g_df["IL-2_Tsec_fraction"].unique()))

########################################################################################################################
print("plotting")
########################################################################################################################

#%%
boxplot_data = []
# fracs = [0.04, 0.076, 0.121, 0.202, 0.301, 0.4]
# fracs = [0.039, 0.077, 0.1245, 0.2005, 0.305, 0.4]
fracs = [0.02, 0.0295, 0.039, 0.058, 0.077, 0.096, 0.115, 0.1435, 0.2005, 0.248, 0.305, 0.4]
grad_boxplot_data = []
for frac in fracs:
    frac_df = big_g_df.loc[(np.abs(big_g_df["IL-2_Tsec_fraction"] - frac) < 1e-3) & (big_g_df["time_index"] == 0)]
    grad_list = []
    for rep in frac_df["replicat_index"].unique():
        rep_df = frac_df.loc[frac_df["replicat_index"] == rep]
        grad_list.append(rep_df["Gradient"].mean() * 1e3) #/ rep_df["surf_c"].mean())
    grad_boxplot_data.append(np.array(grad_list))


np.random.seed(1)
# rc_ticks["figure.figsize"] = (0.95 * 1.25, 1.19)
# rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
# rc_ticks["xtick.labelsize"] = 7
# rc_ticks["ytick.labelsize"] = 6
# rc_ticks['axes.labelsize'] = 7
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

linewidth = 0.5
alpha = 0.9

color1 = "#aa0000ff"


flierprops = dict(marker='o', markerfacecolor=color1, markersize=0,
                  linestyle='none', markeredgecolor=color1, linewidth=linewidth)
boxprops=dict(linewidth=linewidth)
whiskerprops=dict(linewidth=linewidth)
capsprops=dict(linewidth=linewidth)


# ODE_boxp = plt.boxplot(ODE_cv_cv_boxplot_data, flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops, widths=0.5, patch_artist=True)

two_colours = [[color1] * len(fracs)]
from matplotlib.colors import to_rgba
rgba_1 = list(to_rgba(color1))
rgba_1[-1] = 0.6
rgb_colors = [[rgba_1] * len(fracs)]
four_colours = [[color1] * len(fracs) * 2]

for b, bdata in enumerate([grad_boxplot_data]):
    x = np.array(fracs) * 100
    # alphas = np.logspace(0, 0, len(bdata)) if b == 1 else np.linspace(1,1, len(bdata))
    alphas = np.linspace(1,1, len(bdata))
    whisker_alphas = np.array([[a] * 2 for a in alphas]).flatten()
    w = 0.04 * 4
    width = lambda p, w: 10 ** (np.log10(p) + w / 2.) - 10 ** (np.log10(p) - w / 2.)
    boxp = plt.boxplot(bdata, positions=x, flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops, widths=width(x, w), patch_artist=True)

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

# import matplotlib.lines as mlines
# handles = [mlines.Line2D([], [], color=color1, linewidth=2), mlines.Line2D([], [], color=color2, linewidth=2)]
# labels = ["well mixed", "spatial resolution"]
# plt.legend(handles=handles, labels=labels)


plt.ylabel(r"gradient (pM/µm)")
# plt.xlabel("%pSTAT$^+$")
plt.xscale("log")
plt.yscale("linear")
# plt.xlim((0.65, 2.35))

# ax.set_xticklabels(["well\nmixed", "spatial\nresolution"])
# ax.set_xticks([int(x) for x in [5, 10, 15, 20, 30, 40]])
# ax.set_xticklabels([int(x) for x in [5, 10, 15, 20, 30, 40]])
# plt.xlabel(r"fraction T$_{\rm sec}$ (%)")
ax.set_yticks([0, 0.5, 1])

zero_offset_percent = 0
ymax = 1
ymin = -ymax/100 * zero_offset_percent
plt.ylim(ymin, ymax)

plt.xlabel(r"secreting cells (%)")
plt.xlim((3, 51))
plt.xticks([3, 10, 40], [3, 10, 40])

fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()