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

scan_variables = ["act"]

saving_string = r"/home/brunner/Documents/Current work/2024_10_25/boxplot_static_SD_local_q.pdf"
sv_styles = ['-','-.',':']
plot_names = ["sec. scan", "R scan", "KD scan"]
# define the standard values (fold-change = 0) for each of the scan_variables
standards = [0.05]


plot_legend = True


# plotting parameters

yscale = "linear"
yscale_base = 2
xscale = "log"
xscale_base = 2
# define data start and end point
startingPoint = None
stoppingPoint = None
# maximum value
max_value = 8
min_value = 1/4

xlim = (2 ** -2.1, 2 ** 2.1)
ylim = (-1, 2)

# defines a outer layer of N cells to be ignored in plotting. Used to further limit unwanted boundary effects.
offset = 0
# which cell type to plot in which colour. The second colour is the SD colour
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

# base_path = "/{extra}/{u}/paper_models/statics/saturated/{mn}/{n}/".format(u=user, mn=model_name, extra=hdd, n=name)
base_path = "/extra2/brunner/paper_models/kinetics/Figure_1C/dataframes_act_over_Tsec_large_2_steady_state/"
big_c_df, big_g_df = my_load_df(base_path, offset=0, custom_ending="_combined")
big_c_df["IL-2_surf_c"] *= 1e3

try:
    print(np.sort(big_c_df["IL-2_Tsec_fraction"].unique()))
except KeyError:
    big_c_df["IL-2_Tsec_fraction"] = big_c_df["fractions_Tsec"]
    print(np.sort(big_c_df["IL-2_Tsec_fraction"].unique()))
try:
    print(np.sort(big_c_df["IL-2_gamma"].unique()))
except KeyError:
    big_c_df["IL-2_gamma"] = big_c_df["misc_gamma"]
    print(np.sort(big_c_df["IL-2_gamma"].unique()))
if len(big_c_df["IL-2_gamma"].unique()) > 1:
    big_c_df = big_c_df.loc[big_c_df["IL-2_gamma"] == 12.]
########################################################################################################################
print("plotting")
########################################################################################################################
#%%
def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)
# x_axis_label = "log fold change"
# ODE_cv_boxplot_data = [entry for entry in reversed(ODE_cv_boxplot_data)]
#%%
cv_boxplot_data = []
# fracs = [0.05011, 0.1, 0.1584, 0.1995, 0.31622, 0.3981]
# fracs = [0.049, 0.094, 0.148, 0.202, 0.301, 0.4]
# fracs = [0.039, 0.077, 0.1245, 0.2005, 0.305, 0.4]
fracs = [0.02, 0.0295, 0.039, 0.058, 0.077, 0.096, 0.115, 0.1435, 0.2005, 0.248, 0.305, 0.4]
grad_boxplot_data = []
for frac in fracs:
    frac_df = big_c_df.loc[(np.abs(big_c_df["IL-2_Tsec_fraction"] - frac) < 1e-5) & (big_c_df["time_index"] == 0.)]
    cv_list = []
    for rep in frac_df["replicat_index"].unique():
        rep_df = frac_df.loc[frac_df["replicat_index"] == rep]
        cv_list.append(rep_df["IL-2_surf_c"].std())
    cv_boxplot_data.append(np.array(cv_list))



#%%
np.random.seed(0)
# rc_ticks["figure.figsize"] = (0.95 * 1.25, 1.19)
# rc_ticks['figure.figsize'] = (1.67475 * 1.15, 1.386 * 1.15)
rc_ticks['figure.figsize'] = (1.67475 * 1.02, 1.386 * 1.07)
# rc_ticks["xtick.labelsize"] = 7
# rc_ticks["ytick.labelsize"] = 6
# rc_ticks['axes.labelsize'] = 7
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

linewidth = 0.5
alpha = 0.9

color1 = "#aa0000ff"
# color2 = "#e4c137fd"
color2 = "blue"

flierprops = dict(marker='o', markerfacecolor=color1, markersize=0,
                  linestyle='none', markeredgecolor=color1, linewidth=linewidth)
boxprops=dict(linewidth=linewidth)
whiskerprops=dict(linewidth=linewidth)
capsprops=dict(linewidth=linewidth)


# ODE_boxp = plt.boxplot(ODE_cv_cv_boxplot_data, flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops, widths=0.5, patch_artist=True)

two_colours = [[color1] * len(fracs), [color2] * len(fracs)]
from matplotlib.colors import to_rgba
rgba_1 = list(to_rgba(color1))
rgba_1[-1] = 0.3
rgba_2 = list(to_rgba(color2))
rgba_2[-1] = 0.3
rgb_colors = [[rgba_1] * len(fracs),
              [rgba_2] * len(fracs)]
four_colours = [[color1] * len(fracs) * 2, [color2] * len(fracs) * 2]

for b, bdata in enumerate([cv_boxplot_data]):
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


plt.ylabel("surface conc. s.d. (pM)")
# plt.xlabel("%pSTAT$^+$")
plt.xscale("log")
plt.yscale("linear")
# plt.xlim((0.65, 2.35))
plt.ylim(0, 22)
plt.yticks([0, 11, 22])

plt.xlabel(r"secreting cells (%)")
plt.xlim((3, 51))
plt.xticks([3, 10, 40], [3, 10, 40])

fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()