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
path = "/extra2/brunner/paper_models/statics/cells_scan/single_Tsec/"

fig_path = "/home/brunner/Documents/Current work/2023_11_03/"  + "Fig1C_pSTAT5.pdf"

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
global_df = global_df.loc[(global_df["model_name"] == "pde_model")]
cell_df = cell_df.loc[(cell_df["model_name"] == "pde_model")]
cell_df["IL-2_surf_c"] *= 1e3
print(cell_df.geometry_x_grid.unique())
# global_df = activation_to_global_df(cell_df, global_df)
#%%
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
from thesis.scripts.patrick.ODE.statics.parameters import p
from thesis.scripts.patrick.ODE.statics.driver import ODEdriver
def high_density_concentration(p, Rh, N=26):
    # N = 26  # number of sourrounding cells (Kugelschale)
    enum = 4 * np.pi * p["D"] * p["rho"] * (p["L"] + p["rho"]) + p["k_on"] / ODEdr.molar_to_molecules(1) * p["L"] * N * \
           Rh
    denum = p["k_on"] / ODEdr.molar_to_molecules(1) * p["L"] * p["R_Tsec"] * N * Rh + 4 * np.pi * p["D"] * p[
        "rho"] * (p["L"] + p["rho"]) * (p["R_Tsec"] + N * Rh)
    return p["q"] * p["rho"] / (p["k_on"] / ODEdr.molar_to_molecules(1) * p["rho"]) * (enum / denum)

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
p["f_Tsec"] = 1/1000
p["f_Treg"] = 0.0
p["R_Tsec"] = 1
p["R_Treg"] = 5000
p["R_Th"] = 5e3
p["k_endo"] = 0.00046
p["kd"] = 0.1 / 3600
p["IL-2_sigma"] = 1

dr = ODEdriver(p)

# c = dr.molecules_to_molar(
#     I3(p["R_Th"], dr.molar_to_molecules(40 * 1e-12), p["q"], p["f_Tsec"], p["f_Treg"], p["k_endo"], p["R_Tsec"], p["R_Treg"])) * 1e12

# c = dr.molecules_to_molar(
#     I2(np.mean(varied_Rs), dr.molar_to_molecules(40 * 1e-12), p["q"], p["f_Tsec"], p["f_Treg"], p["k_endo"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
#        p["kd"])) * 1e12

def ODE(cells):
    p["N_cells"] = cells
    p["f_Tsec"] = 1 / cells

    varied_Rs = dr.set_R_lognorm(p["R_Th"], p["IL-2_sigma"] * p["R_Th"],
                                 length=int(p["N_cells"] * (1 - p["f_Treg"] - p["f_Tsec"])))
    c = dr.molecules_to_molar(
        I_lin(varied_Rs.mean(), p["q"], p["f_Tsec"], p["f_Treg"], p["k_on"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
              p["kd"])) * 1e12
    return c

def analytical(cells, N=26):
    p["N_cells"] = cells
    p["f_Tsec"] = 1 / cells
    c = ODEdr.molecules_to_molar(high_density_concentration(p, p["R_Th"] * int(p["N_cells"] * (1 - p["f_Treg"] - p["f_Tsec"])), N)) * 1e12
    return c

margin = cell_df.geometry_margin.unique()[0]
L = cell_df.geometry_distance.unique()[0]
def cells_from_x_grid(grid, margin, L):
    reduced_grid = grid - margin
    cells = (reduced_grid/L) ** 3
    return cells

cells = cells_from_x_grid(cell_df.geometry_x_grid.unique(), margin, L)
cell_df["cells"] = cell_df["geometry_x_grid"].map(lambda x: cells_from_x_grid(x, margin, L))

cell_df["normed_cells"] = cell_df["IL-2_surf_c"] * cell_df["cells"]

sns.lineplot(data=cell_df, x = "cells", y = "IL-2_surf_c", label="RD-system", color="blue")
# plt.axhline(ODEdr.molecules_to_molar(high_density_concentration(p, p["R_Th"])) * 1e12, 0, 1, color="black", label="analytical")
plt.plot(cells, [analytical(value, 1) for value in cells], color="black", label="analytical approx.")
plt.plot(cells, [ODE(value) for value in cells], color="red", label="ODE")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("# responding cells")
plt.ylabel("surface conc. (pM)")
plt.yscale("log")

import matplotlib as matplotlib
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
locmaj = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.yaxis.set_major_locator(locmaj)

plt.ylim(0.01, 10)

fig.savefig(fig_path, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()
