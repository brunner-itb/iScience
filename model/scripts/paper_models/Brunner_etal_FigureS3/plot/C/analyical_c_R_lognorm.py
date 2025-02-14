import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
from thesis.scripts.paper_models.utilities.plot_helper import my_load_df
from thesis.scripts.patrick.ODE.driver import ODEdriver
from thesis.scripts.patrick.ODE.parameters import p
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
# rc_ticks["legend.title_fontsize"] = 10
# rc_ticks["axes.titlesize"] = 7
rc_ticks["figure.figsize"] = (1.7, 1.2)

from thesis.scripts.patrick.ODE.statics.parameters import p
from thesis.scripts.patrick.ODE.statics.driver import ODEdriver

from thesis.scripts.paper_models.utilities.plot_helper import my_interpolation
def I_lin(R_resp, q, f_Tsec, f_Treg, k_on, R_Tsec, R_Treg, N_cells, kd):
    A = q * f_Tsec * N_cells
    B = k_on / dr.molar_to_molecules(1) * R_resp * (1 - (f_Tsec + f_Treg)) * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Tsec * f_Tsec * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Treg * f_Treg * N_cells\
        + kd
    return A/B

def high_density_concentration(p, Rh, N=6): # number of sourrounding cells (Kugelschale)
    enum = 4 * np.pi * p["D"] * p["rho"] * (p["L"] + p["rho"]) + p["k_on"] / ODEdr.molar_to_molecules(1) * p["L"] * N * \
           Rh
    denum = p["k_on"] / ODEdr.molar_to_molecules(1) * p["L"] * p["R_Tsec"] * N * Rh + 4 * np.pi * p["D"] * p[
        "rho"] * (p["L"] + p["rho"]) * (
                    p["R_Tsec"] + N * Rh)
    return p["q"] * p["rho"] / (p["k_on"] / ODEdr.molar_to_molecules(1) * p["rho"]) * (enum / denum)

hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
spatial_path = hdd + "/brunner/paper_models/statics/linear/R_hetero_log10/"
spatial_cell_df, spatial_global_df = my_load_df(spatial_path, offset = 0)
spatial_cell_df["time"] = spatial_cell_df["time"].div(3600)
spatial_cell_df["IL-2_surf_c"] = spatial_cell_df["IL-2_surf_c"].mul(1e3)
spatial_cell_df = spatial_cell_df.loc[spatial_cell_df["type_name"] == "Th"]

ODE_path = hdd + "/brunner/paper_models/ODE/linear/R_lognorm_5/"
ODE_cell_df, ODE_global_df = my_load_df(ODE_path, offset = 0)
ODE_cell_df["time"] = ODE_cell_df["time"].div(3600)
ODE_cell_df["IL-2_surf_c"] = ODE_cell_df["IL-2_surf_c"].mul(1e3)
ODE_cell_df = ODE_cell_df[ODE_cell_df["replicat_index"] < 10]
ODE_cell_df = ODE_cell_df.loc[ODE_cell_df["type_name"] == "Th"]
#%%
ODEdr = ODEdriver(p)
p["N_cells"] = 1000
p["R_Th"] = 1e4
p["q"] = 10
p["L"] = 10 # cell-cell-distance surface
p["f_Tsec"] = 0.1
p["f_Treg"] = 0.0
p["R_Tsec"] = 100
p["kd"] = 0.1/3600
dr = ODEdriver(p)

colours = ["red", "black"]

var_log_c = []
var_norm_c = []
ODE_log_c = []
ODE_norm_c = []
vars = np.logspace(-1,1,2000)
log_avg_R = []
norm_avg_R = []
for v, var in enumerate(vars):
    log_Rs = ODEdr.set_R_lognorm_states(1, p["R_Th"], var * p["R_Th"], length=p["N_cells"], write_draws=False).flatten()
    log_avg_R.append(log_Rs.mean())
    norm_Rs = np.random.normal(p["R_Th"], var * p["R_Th"], p["N_cells"])
    norm_Rs = norm_Rs[norm_Rs > 0]
    norm_avg_R.append(norm_Rs.mean())
    # norm_Rs[np.where(norm_Rs < 0)] = p["R_Th"]
    log_c = ODEdr.molecules_to_molar(high_density_concentration(p, log_Rs, N = 9)) * 1e12
    norm_c = ODEdr.molecules_to_molar(high_density_concentration(p, norm_Rs, N = 9)) * 1e12

    var_log_c.append(np.mean(log_c))
    var_norm_c.append(np.mean(norm_c))

    # ODE_log_c.append(dr.molecules_to_molar(
    #     I_lin(log_Rs.mean(), p["q"], p["f_Tsec"], p["f_Treg"], p["k_on"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
    #           p["kd"])) * 1e12)
    # ODE_norm_c.append(dr.molecules_to_molar(
    #     I_lin(norm_Rs.mean(), p["q"], p["f_Tsec"], p["f_Treg"], p["k_on"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
    #           p["kd"])) * 1e12)

#%%
saving_string =r"/home/brunner/Documents/Current work/2023_11_03/" + "analytical_R_lognorm" + ".pdf"
analytical_df = pd.DataFrame()
analytical_df["IL-2_sigma"] = vars
analytical_df["IL-2_surf_c"] = var_log_c
analytical_df["replicat_index"] = 0

analytical_df["IL-2_sigma"] = pd.cut(analytical_df["IL-2_sigma"], bins=np.linspace(analytical_df["IL-2_sigma"].min(), analytical_df["IL-2_sigma"].max(), 30), include_lowest=True, labels=np.linspace(analytical_df["IL-2_sigma"].min(), analytical_df["IL-2_sigma"].max(), 29))

sns.set_style("ticks")
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig,ax = plt.subplots()
s=0.5
N=80
# for idx, c in enumerate([var_log_c]):
#     # sns.scatterplot(vars, np.array(c), s=s, color=colours[idx])
#     sns.lineplot(vars[N:], np.convolve(c, np.ones(N)/N, mode='same')[N:], color=colours[idx])
#     sns.lineplot(vars[:N], c[:N], color=colours[idx], label="analytical")

# colors = ["#aa0000ff", "black", "green", ]
colors = ["black", "black", "green", ]
labels = ["RD-system", "well-mixed", "analytical approx."]
for c, c_df in enumerate([spatial_cell_df, ODE_cell_df, analytical_df]):
    run_mean = []
    run_ste = []
    for s, sigma in enumerate(c_df["IL-2_sigma"].unique()):
        s_df = c_df.loc[c_df["IL-2_sigma"] == sigma]
        mean = []
        # ste = []
        for r, rep in enumerate(s_df.replicat_index.unique()):
            mean.append(s_df.loc[s_df.replicat_index == rep, "IL-2_surf_c"].mean())
            # ste.append(s_df.loc[s_df.replicat_index == rep, "IL-2_surf_c"].std()/ np.sqrt(len(s_df.replicat_index.unique())))
        run_mean.append(np.mean(mean))
        run_ste.append(np.std(mean))# / np.sqrt(len(s_df.replicat_index.unique())))
    x = c_df["IL-2_sigma"].unique()
    plt.plot(x, run_mean, color=colors[c], label=labels[c], linestyle = "-" if c != 1 else "--")
    plt.fill_between(x, np.array(run_mean) + np.array(run_ste), np.array(run_mean) - np.array(run_ste), color=colors[c], linewidth=0, alpha=0.3)


ax.set(ylabel=r"surface conc. avg. (pM)", xlabel="receptor heterogeneity", yscale="log", xscale="log", ylim=(1, 300), xlim=(None, None)) #, title=None)

ax.set_xlim((0.1, 10))
ax.tick_params(axis='y', colors=colors[0], which="both")
ax.yaxis.label.set_color(colors[0])
# plt.title("linear uptake")
# from matplotlib.lines import Line2D
# handles, labels = ax.get_legend_handles_labels()
# ax.get_legend().remove()
# plt.legend()
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()