import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = False
import numpy as np
import pandas as pd
import seaborn as sns
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.scripts.paper_models.utilities.plot_helper import EC50_calculation

# rc_ticks['figure.figsize'] = (1.2, 1.1)
rc_ticks['figure.figsize'] = (1.67475, 1.386)

# path = "/extra2/brunner/paper_models/statics/receptor_distribution/"
# path = "/extra2/kiwitz/20220727_paper_models/receptor_distribution_box_240_2d/test_240_fsec_linear_1/"
path = "/extra2/kiwitz/20220727_paper_models/receptor_distribution_box_240_2d/test_240_fsec_linear_5_lognorm/"

fig_path = "/home/brunner/Documents/Current work/2024_10_25/"

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df["SD"] *= 1e3
global_df["surf_c_std"] *= 1e3
cell_df["IL-2_surf_c"] *= 1e3

cell_df.v_v.unique() #receptor distribution
cmap = sns.cm.rocket_r
# cell_df.numeric_ratio.unique()  #ratio between sec and abs


#%%
for plot_df in [cell_df]:
    plot_df = plot_df.loc[plot_df["type_name"] != "sec"]
    plot_df["act_cells"] = None
    plot_df["pSTAT"] = plot_df["IL-2_surf_c"] ** 3 / (
            (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=plot_df["IL-2_R"]) * 1e12) ** 3 + plot_df[
        "IL-2_surf_c"] ** 3).values
    for f, frac in enumerate(plot_df["fractions_sec"].unique()):
        for v, v_v in enumerate(plot_df["v_v"].unique()):
            act_df = plot_df.loc[(plot_df["fractions_sec"] == frac) & (plot_df["v_v"] == v_v)]
            act_cells = len(act_df.loc[act_df["pSTAT"] >= 0.5])/len(act_df)
            plot_df.loc[(plot_df["fractions_sec"] == frac) & (plot_df["v_v"] == v_v), "act_cells"] = act_cells

    plot_df["act_cells"] *= 100
    plot_df["fractions_sec"] *= 100
    pt = pd.pivot_table(plot_df, values="act_cells", columns=["v_v"], index=["fractions_sec"], aggfunc=np.mean)
    pt = pt.iloc[:, ::-1] #reverse columns (x-axis)
    pt = pt[::-1] #and indices (y-axis)

    sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
    fig, ax = plt.subplots()
    from matplotlib.colors import LogNorm
    color = "#aa0000ff"
    ax = sns.heatmap(pt, cmap=sns.light_palette(color, as_cmap=True), cbar_kws={"label": r"%pSTAT5$^+$"}, vmax=100)

    ax.set(xlabel="% receptors on T$_{resp}$", ylabel="secreting cells (%)")
    len_x_ticks = len(plot_df["v_v"].unique())
    plt.xticks([0, len_x_ticks/4, len_x_ticks/2, len_x_ticks * 3/4, len_x_ticks], [0, 25, 50, 75, 100])
    len_y_ticks = len(plot_df["fractions_sec"].unique())
    plt.yticks([0, len_y_ticks/4, len_y_ticks/2, len_y_ticks * 3/4, len_y_ticks],
               np.around([np.max(plot_df["fractions_sec"].unique()), np.max(plot_df["fractions_sec"].unique())* 3/4,
                          np.max(plot_df["fractions_sec"].unique())/2, np.max(plot_df["fractions_sec"].unique()) * 1/4,
                          np.min(plot_df["fractions_sec"].unique())], 1))
    plt.xticks(rotation=0, ha='center')
    plt.xlim((len_x_ticks/2, len_x_ticks))
    plt.ylim((15, 0))
    fig.savefig(fig_path + "Fig1_D_receptor_distribution_surf_c_ACTIVATION.pdf", bbox_inches='tight', transparent=True)
    plt.tight_layout()
    plt.show()


# for plot_df in [cell_df]:
#     plot_df = plot_df.loc[plot_df["type_name"] != "sec"]
#     plot_df["act_cells"] = None
#     plot_df["pSTAT"] = plot_df["IL-2_surf_c"] ** 3 / (
#             (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=plot_df["IL-2_R"]) * 1e12) ** 3 + plot_df[
#         "IL-2_surf_c"] ** 3).values
#     fracs = plot_df["fractions_sec"].unique()
#     fracs_range = fracs[(fracs > 0.08) & (fracs < 0.15)]
#
#     vs = plot_df["v_v"].unique()
#     vs_range = vs[(vs > 0) & (vs <= 0.5)]
#     for f, frac in enumerate(fracs_range):
#         for v, v_v in enumerate(vs_range):
#             act_df = plot_df.loc[(plot_df["fractions_sec"] == frac) & (plot_df["v_v"] == v_v)]
#             act_cells = len(act_df.loc[act_df["pSTAT"] >= 0.5])/len(act_df)
#             plot_df.loc[(plot_df["fractions_sec"] == frac) & (plot_df["v_v"] == v_v), "act_cells"] = act_cells
#
#     plot_df["act_cells"] *= 100
#     pt = pd.pivot_table(plot_df, values="act_cells", columns=["v_v"], index=["fractions_sec"], aggfunc=np.mean)
#     pt = pt.iloc[:, ::-1] #reverse columns (x-axis)
#     pt = pt[::-1] #and indices (y-axis)
#
#     sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
#     fig, ax = plt.subplots()
#     ax = sns.heatmap(pt, cmap=cmap, cbar_kws={"label": "% activated cells"})
#
#     ax.set(xlabel="% receptors on T$_{resp}$", ylabel="fraction of sec. cells")
#     len_x_ticks = len(vs_range)
#     plt.xticks([0, len_x_ticks/2, len_x_ticks], np.around([(1 - vs_range)[-1] * 100, (1 - vs_range)[-int(len_x_ticks/2)] * 100, (1 - vs_range)[0]], 3))
#     len_y_ticks = len(fracs_range)
#     plt.yticks([0, len_y_ticks/2, len_y_ticks], np.around([np.max(fracs_range), np.max(fracs_range)/2, np.min(fracs_range)], 3))
#
#     fig.savefig(fig_path + "Fig1_D_receptor_distribution_surf_c_ACTIVATION_inset.pdf", bbox_inches='tight', transparent=True)
#     plt.tight_layout()
#     plt.show()


# for plot_df in [global_df]:
#     plot_df["Gradient"] *= 1e3
#     pt = pd.pivot_table(plot_df, values="Gradient", columns=["v_v"], index=["fractions_sec"], aggfunc=np.mean)
#     pt = pt.iloc[:, ::-1]  # reverse columns (x-axis)
#     pt = pt[::-1]  # and indices (y-axis)
#
#     sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
#     fig, ax = plt.subplots()
#     ax = sns.heatmap(pt, cmap=cmap, cbar_kws={"label": "gradient (pM/µm)"})
#
#     ax.set(xlabel="% receptors on T$_{resp}$", ylabel="fraction of sec. cells")
#
#     len_x_ticks = len(plot_df["v_v"].unique())
#     plt.xticks([0, len_x_ticks / 2, len_x_ticks], [0, 50, 100])
#     len_y_ticks = len(plot_df["fractions_sec"].unique())
#     plt.yticks([0, len_y_ticks / 2, len_y_ticks],
#                [np.max(plot_df["fractions_sec"].unique()), np.max(plot_df["fractions_sec"].unique()) / 2,
#                 np.min(plot_df["fractions_sec"].unique())])
#
#     fig.savefig(fig_path + "Fig1_D_receptor_distribution_GRAD_CONTOUR.pdf", bbox_inches='tight', transparent=True)
#     plt.tight_layout()
#     plt.show()
fig, ax = plt.subplots()
plt.plot([0], [0])
plt.ylim(0,100)
plt.xlim(50,100)
plt.yticks([0,50,100])
plt.xticks([50,75,100])
plt.ylabel("secreting cells (%)")
plt.xlabel(r"% receptors on T$_{\rm resp}$")
fig.savefig(fig_path + "test.pdf", bbox_inches='tight', transparent=True)
plt.show()