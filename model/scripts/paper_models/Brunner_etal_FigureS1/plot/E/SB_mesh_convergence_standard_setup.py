import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.scripts.paper_models.utilities.plot_helper import activation_to_global_df
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
# path = hdd + "/brunner/paper_models/prelim_boxed_static/parameter_scan_sigma/"
# path = "/extra2/brunner/paper_models/boxed_static/mesh_convergence_linear_2/"
frac = 0.1

path = f"/extra2/brunner/paper_models/statics/cells_scan/{str(frac)}_Tsecs_linear_extension/"

fig_path = "/home/brunner/Documents/Current work/2023_11_03/"

global_df = pd.read_hdf(path + 'global_df_combined.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df_combined.h5', mode="r")
global_df = global_df.loc[(global_df["model_name"] == "pde_model")]
cell_df = cell_df.loc[(cell_df["model_name"] == "pde_model")]
cell_df["IL-2_surf_c"] *= 1e3

# path2 = f"/extra2/brunner/paper_models/statics/cells_scan/0.1_Tsecs_linear_extension/"
# global_df_2 = pd.read_hdf(path2 + 'global_df.h5', mode="r")
# cell_df_2 = pd.read_hdf(path2 + 'cell_df.h5', mode="r")
#
# combined_df = pd.concat([cell_df, cell_df_2])
# combined_df.to_hdf(path2 + "cell_df_combined.h5", key="data")
# combined_df = pd.concat([global_df, global_df_2])
# combined_df.to_hdf(path2 + "global_df_combined.h5", key="data")
# global_df = activation_to_global_df(cell_df, global_df)
#%%
print(cell_df["IL-2_q"].unique())
print(cell_df["IL-2_k_on"].unique())
print(cell_df["fractions_Tsec"].unique())
print(cell_df.loc[cell_df.type_name == "Th", "IL-2_R"].mean())
#%%
from thesis.scripts.patrick.ODE.statics.parameters import p
from thesis.scripts.patrick.ODE.statics.driver import ODEdriver
def high_density_concentration(p, Rh, N, N_Tsecs):
    enum = 4 * np.pi * p["D"] * p["rho"] * (p["L"] + p["rho"]) + p["k_on"] / ODEdr.molar_to_molecules(1) * p["L"] * N * \
           Rh
    denum = p["k_on"] / ODEdr.molar_to_molecules(1) * p["L"] * p["R_Tsec"] * N * Rh + 4 * np.pi * p["D"] * p[
        "rho"] * (p["L"] + p["rho"]) * (p["R_Tsec"] + N * Rh)
    return p["q"] * N_Tsecs * p["rho"] / (p["k_on"] / ODEdr.molar_to_molecules(1) * p["rho"]) * (enum / denum)

def I_lin(R_resp, q, f_Tsec, f_Treg, k_on, R_Tsec, R_Treg, N_cells, kd):
    A = q * f_Tsec * N_cells
    B = k_on / dr.molar_to_molecules(1) * R_resp * (1 - (f_Tsec + f_Treg)) * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Tsec * f_Tsec * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Treg * f_Treg * N_cells\
        + kd
    return A/B

def I3(R_resp, K_D, q, f_Tsec, f_Treg, k_endo, R_Tsec, R_Treg):
    a = q * f_Tsec
    b = k_endo * (R_resp * (1 - (f_Tsec + f_Treg)) + R_Tsec * f_Tsec + R_Treg * f_Treg)
    return -(a * K_D)/(a - b)

def I2(R_resp, K_D, q, f_Tsec, f_Treg, k_endo, R_Tsec, R_Treg, N_cells, kd):
    a = q * N_cells * f_Tsec
    b = k_endo  * N_cells * (R_resp * (1 - (f_Tsec + f_Treg)) + R_Tsec * f_Tsec + R_Treg * f_Treg)
    return (np.sqrt(kd ** 2 * K_D ** 2 + 2 * kd * K_D * (a + b) + (a - b) ** 2) - kd * K_D + a - b) / (2 * kd)

ODEdr = ODEdriver(p)

p["N_R"] = 3
# p["k_endo"] = 0.0011
# # p["K_D"] = 58805.37
p["q"] = 10
p["L"] = 10 # cell-cell-distance surface
p["N_cells"] = 1000
p["f_Tsec"] = frac
p["f_Treg"] = 0.0
p["R_Tsec"] = 1
p["R_Treg"] = 5000
p["R_Th"] = 1e4
p["k_endo"] = np.log(2)/(25*60)
p["kd"] = 0.1/3600
p["IL-2_sigma"] = 1

dr = ODEdriver(p)

def ODE(cells):
    p["N_cells"] = 1000
    dr = ODEdriver(p)
    K_D = p["k_off"]/(p["k_on"]/dr.molar_to_molecules(1))
    real_f_Tsec = np.round(p["N_cells"] * p["f_Tsec"], 0)/p["N_cells"]

    varied_Rs = dr.set_R_lognorm(p["R_Th"], p["IL-2_sigma"] * p["R_Th"],
                                 length=int(p["N_cells"] * (1 - p["f_Treg"] - real_f_Tsec)))
    # c = dr.molecules_to_molar(
    #     I3(p["R_Th"], K_D, p["q"], real_f_Tsec, p["f_Treg"], p["k_endo"], 0, p["R_Treg"])
    #     ) * 1e12
    c = dr.molecules_to_molar(
        I_lin(varied_Rs.mean(), p["q"], real_f_Tsec, p["f_Treg"], p["k_on"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
              p["kd"])) * 1e12
    if c <= 0:
        c = np.nan
    return c

def analytical(cells, N=26):
    p["N_cells"] = cells
    real_f_Tsec = np.round(cells * p["f_Tsec"], 0)/cells
    c = ODEdr.molecules_to_molar(
        high_density_concentration(p, p["R_Th"], int(p["N_cells"] * (1 - p["f_Treg"] - p["f_Tsec"])), N_Tsecs = int(p["N_cells"] * p["f_Tsec"]))) * 1e12
    return c


margin = cell_df.geometry_margin.unique()[0]
L = cell_df.geometry_distance.unique()[0]
def cells_from_x_grid(grid, margin, L):
    reduced_grid = grid - margin
    cells = (reduced_grid/L) ** 3
    return cells

cells = cells_from_x_grid(cell_df.geometry_x_grid.unique(), margin, L)[0:]
cell_df["cells"] = cell_df["geometry_x_grid"].map(lambda x: cells_from_x_grid(x, margin, L))

cell_df["normed_cells"] = cell_df["IL-2_surf_c"] * cell_df["cells"]

#%%
cell_df["IL-2_surf_c"] -= 1.3
#%%

sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
alpha = 1
RD = []
RD_std = []

desired_frac = cell_df["fractions_Tsec"].unique()[0]
for c, cell in enumerate(cells):
    c_df = cell_df.loc[cell_df.cells == cell]
    rep_RD = []
    for r, rep in enumerate(c_df.replicat_index.unique()):
        rep_df = c_df.loc[c_df["replicat_index"] == rep]
        real_frac = len(rep_df.loc[rep_df.type_name == "Tsec"])/len(rep_df)
        rep_RD.append(np.mean(rep_df["IL-2_surf_c"]) * desired_frac/real_frac)
    RD.append(np.mean(rep_RD))
    RD_std.append(np.std(rep_RD))

RD = np.array(RD)
RD_std = np.array(RD_std)/len(c_df.replicat_index.unique())
plt.plot(cell_df.cells.unique(), RD, label="RD-system", color="#e4c137fd", alpha=alpha, marker="o", markersize = 2)
plt.fill_between(cell_df.cells.unique(), RD - RD_std, RD + RD_std, color="#e4c137fd", alpha=alpha/4)
plt.plot(cells, [np.mean([ODE(value) for x in range(2000)]) for v,value in enumerate(cells)], color="black", label="ODE",
         alpha=alpha, marker="o", markersize = 2)
plt.plot(cells, [analytical(value, 1) for value in cells], color="green", linestyle="--", label="analytical approx.", marker="o", markersize = 2)

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("system size (# cells)")
plt.ylabel("surface conc. (pM)")

plt.xscale("log")

plt.ylim(0, 5)
plt.minorticks_off()
plt.xticks([N**3 for N in np.arange(3, 11)], [rf"{N}$^3$" for N in np.arange(3, 11)])
plt.xlim(27, 10**3)
# plt.yscale("log")
plt.yticks([0, 2.5, 5])
ax.get_legend().remove()
fig.savefig(fig_path + "FigS1_mesh_convergence_surf_c_standard_setup.pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()

#%%

sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
alpha = 1
y = []
std = []

cells = cells_from_x_grid(cell_df.geometry_x_grid.unique(), margin, L)
analytical_sol = [analytical(value) for value in cells]

RD = []
RD_std = []
desired_frac = cell_df["fractions_Tsec"].unique()[0]
for c, cell in enumerate(cells):
    c_df = cell_df.loc[cell_df.cells == cell]
    rep_RD = []
    for r, rep in enumerate(c_df.replicat_index.unique()):
        rep_df = c_df.loc[c_df["replicat_index"] == rep]
        real_frac = len(rep_df.loc[rep_df.type_name == "Tsec"])/len(rep_df)
        rep_RD.append(np.mean(rep_df["IL-2_surf_c"] * desired_frac/real_frac - analytical_sol[c]))
    RD.append(np.mean(rep_RD))
    RD_std.append(np.std(rep_RD))
RD = np.array(RD)
RD_std = np.array(RD_std)/len(c_df.replicat_index.unique())
# sns.lineplot(data=cell_df, x = "cells", y = "IL-2_surf_c", label="RD-system", color="blue")
# plt.axhline(ODEdr.molecules_to_molar(high_density_concentration(p, p["R_Th"])) * 1e12, 0, 1, color="black", label="analytical")
plt.plot(cells, RD, color="#e4c137fd", label="RD-system", marker="o", markersize = 2)
plt.fill_between(cells, RD - RD_std, RD + RD_std, color="#e4c137fd", alpha=alpha/3)
plt.plot(cells, [np.mean([ODE(value) - analytical_sol[v] for x in range(2000)]) for v,value in enumerate(cells)], 
         color="black", label="ODE", marker="o", markersize = 2)


ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("system size (# cells)")
plt.ylabel(r"difference to" "\n" "analytical solution (pM)")

plt.xscale("log")
# plt.ylim(1.8, 7)
plt.minorticks_off()
plt.xticks([N**3 for N in np.arange(3, 11)], [rf"{N}$^3$" for N in np.arange(3, 11)])
plt.xlim(27, 10**3)
# plt.yscale("log")
plt.yticks([0, 1, 2])
plt.ylim(-0.0, 2)
ax.get_legend().remove()
fig.savefig(fig_path + "FigS1_mesh_convergence_surf_c_diff_to_analytical_standard_setup.pdf", bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()