import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = False
import numpy as np
import pandas as pd
import seaborn as sns
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

# rc_ticks['figure.figsize'] = (1.2, 1.1) #(1.7, 1.2)

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

path = "/extra2/brunner/paper_models/statics/receptor_distribution_box_240/20220905_final_quality_3/"
path = "/extra2/brunner/paper_models/statics/receptor_distribution/dataframes__steady_state/"

fig_path = "/home/brunner/Documents/Current work/2024_10_25/"

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")

print(cell_df["IL-2_surf_c"].mean() * 1e3)
print(cell_df.loc[cell_df["scan_index"] == 12, "IL-2_R"].mean())

cell_df["IL-2_surf_c"] *= 1e3
global_df["surf_c"] *= 1e3
global_df["surf_c_std"] *= 1e3
global_df["Gradient"] *= 1e3
global_df["cv"] = global_df["surf_c_std"] / global_df["surf_c"]


print("pSTAT5 values not available, calculating")
cell_df["pSTAT5"] = cell_df["IL-2_surf_c"] ** 3 / (
    (EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=cell_df["IL-2_R"]) * 1e12) ** 3 + cell_df[
            "IL-2_surf_c"] ** 3).values

my_hue = "scan_value"

vs = np.sort(cell_df[my_hue].unique())
vs = vs[~np.isnan(vs)]
# cell_df = cell_df.loc[cell_df["type_name"] == "Th"]
cell_df = cell_df.loc[cell_df["model_name"] == "pde_model"]

pSTAT5 = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))
Th_activated_cells = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))
Th_activated_cells[:] = None
Tsec_activated_cells = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))
Tsec_activated_cells[:] = None
SD = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))
mean = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))
grad = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))
CV = np.zeros((len(cell_df["replicat_index"].unique()), len(vs)))

for r, rep in enumerate(cell_df["replicat_index"].unique()):
    rep_df = cell_df.loc[cell_df["replicat_index"] == rep]
    global_rep_df = global_df.loc[global_df["replicat_index"] == rep]
    for idx,v in enumerate(reversed(vs)):
        scan_df = rep_df.loc[(rep_df[my_hue] == v)]
        act_Th = len(scan_df.loc[(scan_df.type_name == "Th") & (scan_df["pSTAT5"] >= 0.5), "id_id"].unique())
        act_Tsec = len(scan_df.loc[(scan_df.type_name == "Tsec") & (scan_df["pSTAT5"] >= 0.5), "id_id"].unique())
        SD[r][idx] = (scan_df["IL-2_surf_c"].std())
        mean[r][idx] = (scan_df["IL-2_surf_c"].mean())
        # mean.append(global_rep_df.loc[(global_rep_df[my_hue] == v), "surf_c"].mean())
        Th_activated_cells[r][idx] = (act_Th/len(scan_df["id_id"].unique()))
        Tsec_activated_cells[r][idx] = (act_Tsec/len(scan_df["id_id"].unique()))
        pSTAT5[r][idx] = (scan_df["pSTAT5"].mean())
        grad[r][idx] = (global_rep_df.loc[(global_rep_df[my_hue] == v) , "Gradient"].mean()) #/global_rep_df.loc[(global_rep_df[my_hue] == v) , "surf_c"].mean())
        CV[r][idx] = (global_rep_df.loc[(global_rep_df[my_hue] == v) , "cv"].mean())

vs *= 100

#%%
# plotting_var, var_name = (grad, "grad")
plotting_var, var_name = (SD, "SD")
# plotting_var, var_name = (SD/mean, "CV")

plots_ax3 = [
    Th_activated_cells
]
ax3_colours = ["#aa0000ff", "#e5c237f6", "green"]
labels = ["Tresp", "Tsec", "All cells"]

if var_name == "grad" or var_name == "SD":
    # rc_ticks['figure.figsize'] = (1.67475 , 1.386)
    rc_ticks['figure.figsize'] = (1.67475 * 0.57 * 1.3, 1.2 * 1.3)
    pass
else:
    rc_ticks['figure.figsize'] = (1.67475 * 0.57 * 1.1, 1.2 * 1.1)
# rc_ticks['figure.figsize'] = (5, 5)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig, ax = plt.subplots()

CV_colour = "black"
# colour2 = sns.color_palette("dark", 9)[-1]
ax2 = sns.lineplot(x=vs, y=np.mean(plotting_var, axis=0), color=CV_colour)
neg_fill = np.mean(plotting_var, axis=0) - np.std(plotting_var, axis=0)/np.sqrt(len(cell_df.replicat_index.unique()))
pos_fill = np.mean(plotting_var, axis=0) + np.std(plotting_var, axis=0)/np.sqrt(len(cell_df.replicat_index.unique()))

# neg_fill = [x for x in reversed(neg_fill)]
# pos_fill = [x for x in reversed(pos_fill)]
plt.fill_between(vs, neg_fill, pos_fill,color=CV_colour, alpha=0.15, linewidth=0)

ax2.tick_params(axis='y', colors=CV_colour)
ax2.yaxis.label.set_color(CV_colour)
# ax2.set(ylabel="gradient (pM/µm)")


ax3 = ax2.twinx()
for e, entry in enumerate(plots_ax3):
    ax3.plot(vs, np.mean(entry, axis=0) * 100, color=ax3_colours[e], label=labels[e])
    neg_fill = (np.mean(entry, axis=0) * 100) - np.std(entry, axis=0)/np.sqrt(len(cell_df.replicat_index.unique()))
    pos_fill = (np.mean(entry, axis=0) * 100) + np.std(entry, axis=0)/np.sqrt(len(cell_df.replicat_index.unique()))
    plt.fill_between(vs, neg_fill, pos_fill, color=ax3_colours[e], alpha=0.15, linewidth=0)

ax3.set_xlabel("% receptors on T$_{resp}$")
ax3.set_ylabel(r"%pSTAT$^+$ T$_{\rm resp}$ cells")
# ax3.set(ylabel="avg. pSTAT")
ax3.tick_params(axis='y', colors=ax3_colours[0])
ax3.yaxis.label.set_color(ax3_colours[0])

# plt.legend()
if var_name == "grad":
    # ax2.set_yticks([0, 0.1, 0.2, 0.3])
    ax2.set(xlabel="% receptors on T$_{resp}$", ylabel=r"gradient (pM/µm)", ylim=(0, 0.4), yticks=[0, 0.2, 0.4], yticklabels=[0, 0.2, 0.4])
    ax3.set(ylim=(-1,52), yticks=[0, 25, 50], yticklabels=[0, 25, 50])
    fig.savefig(fig_path + "Fig1_D_receptor_distribution_grad.pdf", bbox_inches='tight', transparent=True)
elif var_name == "CV":
    ax2.set_yticks([0, 0.4, 0.8])
    # ax3.set_yticks([0, 15, 30])
    ax2.set_xticks([0, 50, 100])
    ax2.set(xlabel=r"% receptors on T$_{\rm resp}$", ylabel=r"CV")
    fig.savefig(fig_path + "Fig1_D_receptor_distribution_CV.pdf", bbox_inches='tight', transparent=True)
elif var_name == "SD":
    # ax2.set_yticks([0, 0.4, 0.8])
    # ax3.set_yticks([0, 15, 30])
    # ax2.set_xticks([0, 50, 100])
    ax2.set(xlabel=r"% receptors on T$_{\rm resp}$", ylabel=r"surface conc. s.d. (pM)", ylim=(0, 8), yticks=[0, 4, 8], yticklabels=[0, 4, 8])
    ax3.set(ylim=(-1,52), yticks=[0, 25, 50], yticklabels=[0, 25, 50])
    fig.savefig(fig_path + "Fig1_D_receptor_distribution_SD.pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

# cell_df["misc_v"].unique()
# cell_df.loc[(cell_df["misc_v"] == 0.), "IL-2_R"].mean()
