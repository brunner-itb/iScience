import getpass
import os

import numpy as np

path = "/extra2/brunner/paper_models/kinetics/Tsec_scan_25/" # where to save the results of this run
ext_cache = r"../Tsec_scan_25_3D_ext_cache/" # cache location, contains save states, the mesh etc.

"""defines cytokines. Currently their parameters can be changed using the scan sample interface.
"field_quantity" determines which boundary contitions interact with which fields, 
the name is only used in IO/post processing."""

cytokines = [
    {
        "name": "IL-2", # name of the cytokine
        "field_quantity": "il2", # internal name of the cytokine
        "k_on": 111.6,  # receptor binding constant 1/(nM*h),
        "D": 10,  # Diffusion constant mu^2
        "kd": 0.1,  # cytokine decay in medium in 1/h, Η in the paper...
        "k_endo": 0.46e-3, # endocytosis rate in 1/s (oversight, can be fixed in the box_grid.py)
        "k_off": 0.83 # receptor dissociation rate in 1/h
    }
]

"""Sets up cells types. 
The first entry is the default cell type. The "fraction" entry is meaningless.
"""

cell_types_dict = [
    {

     "name": "Tnaive",
     "fraction": 0,
     "il2": {"R": 1e2, "q": 0, "bc_type": "patrick_saturation", "global_q": False},  # [Receptor number per cell, secretion in molecules/s, bc_type defines boundary condition type which is either linear or saturated, global_q makes the secretion systemic, as in independent of number of Tsecs]
     "misc": {"sigma": 1, # Receptor heterogeneity
              "states":[], # saves the cell states in the linear chain trick. Was used to set initial conditions for the linear chain trick, not sure if it still works properly. If empty the initial state of the cell is chosen
              "hill_factor": 3, # hill factor for the receptor feedback calculation
              "Km_pos": 0.5, # The receptors of half saturation in the receptor feedback function during postive fb, K_m in the paper, Eq. 5
              "Km_neg": 0.5, # equivalent for negative fb
              "pSTAT5_signal": True, # whether to calculate the pSTAT5 signal
              "KD": 7.437e-3, # concentration of half saturation for the uptake saturation, K_D in the paper
              # "nu": 1e-3, # receptor decay, deprecated, can be safely removed. Removed from paper in the correction.
              "name": "Tnaive"},
     "internal_solver": "kineticSolver"
     },
    {
     "name": "Tsec",
     "fraction": 0.05, # fraction of secreting cells, f_sec in the paper
     "il2": {"R": 1e2, "q": 10, "bc_type": "patrick_saturation", "global_q": False},# [Receptor number per cell -> R_sec in the paper, secretion in molecules/s, bc_type defines boundary condition type which is either linear or saturated, global_q makes the secretion systemic, as in independent of number of Tsecs]
     "misc": {"sigma": 1, # Receptor heterogeneity, also sigma in the paper.
              "states":[], # saves the cell states in the linear chain trick. Was used to set initial conditions for the linear chain trick, not sure if it still works properly. If empty the initial state of the cell is chosen
              "hill_factor": 3, # hill factor for the receptor feedback calculation, N_pSTAT5 in the paper, Eq. 5
              "Km_pos": 0.5, # The receptors of half saturation in the receptor feedback function during postive fb, K_m in the paper, Eq. 5
              "Km_neg": 0.5, # equivalent for negative fb
              "pSTAT5_signal": True, # whether to calculate the pSTAT5 signal
              "KD": 7.437e-3, # concentration of half saturation for the uptake saturation, K_D in the paper
              "name": "Tsec"},
     "internal_solver": "kineticSolver"
     },
    {
     "name": "Th", # name of the cell in the physical parameter set
     "fraction": 0.95, # fraction of responding cells, f_resp in the paper
     "il2": {"R": 1e3, "q": 0, "bc_type": "patrick_saturation", "global_q": False, "pSTAT5": 0}, # [Receptor number per cell -> R_resp in the paper, secretion in molecules/s, bc_type defines boundary condition type which is either linear or saturated, global_q makes the secretion systemic, as in independent of number of Tsecs]
     "misc": {"sigma": 1.5, # Receptor heterogeneity
              "states":[0,0,0,0,0,0],# saves the cell states in the linear chain trick. Was used to set initial conditions for the linear chain trick, not sure if it still works properly. If empty the initial state of the cell is chosen
              "gamma": 10, # feedback strength
              "pos_half_time": 1, # positive half time parameter in h, is calculated into lambda in the paper Eq. 6 via "np.log(2)/(half_time*3600)"
              "neg_half_time": 1, # negative equivalent
              "hill_factor": 3, # hill factor for the receptor feedback calculation, N_pSTAT5 in the paper, Eq. 5
              "Km_pos": 0.5,# The receptors of half saturation in the receptor feedback function during postive fb, K_m in the paper, Eq. 5
              "Km_neg": 0.5, # equivalent for negative fb
              "R_start_pos": 1.5e3, # initial receptors in the positive fb case, R_resp in the paper. Is only used in the feedback case, else remains at the first given "R"
              "R_start_neg": 1.5e3, # initial receptors in the negative fb case, R_resp in the paper. Is only used in the feedback case, else remains at the first given "R"
              "pSTAT5_signal": True, # whether to calculate the pSTAT5 signal
              "KD": 7.437e-3, # concentration of half saturation for the uptake saturation, K_D in the paper
              "pSTAT5": 0, # container to save the pSTAT5 value in, gets calculated and saved to the cell_df during postprocessing
              "name": "Th", # name of the cell in the misc parameter set, just to be safe. Can be omitted
              "EC50_N": 1.5}, # hill factor for the EC50 calculation, N_EC50 in the paper
     "internal_solver": "kineticSolver"
     },
    {
        "name": "Treg",
        "fraction": 0.0,
        "il2": {"R": 1e3, "q": 0, "bc_type": "patrick_saturation", "global_q": False, "pSTAT5": 0},
        "misc": {"sigma": 0.5,
                 "states": [0, 0, 0, 0, 0, 0],
                 "pos_half_time": 1,
                 "neg_half_time": 1,
                 "hill_factor": 3,
                 "Km_pos": 0.5,
                 "Km_neg": 0.5,
                 "R_start_neg": 1e5,
                 "R_start_pos": 1e4,
                 "pSTAT5_signal": True,
                 "KD": 7.437e-3,
                 "pSTAT5": 0,
                 "name": "Treg"},
        "internal_solver": "kineticSolver"
    },
]

"""defines the variable aspects of the geometry. Unit is micro meters"""
geometry = {
    "margin": 40,  # margin around the cell grid in µm
    "distance": 20,  # distance between cell centers in µm, d_c in the paper
    "rho": 5,  # cell radius in µm, r_c in the paper
    "x_grid": 240,  # dimensions of the cell grid in µm
    "y_grid": 240,
    "z_grid": 240,# comment out for single cell layer
    "norm_area": 4 * np.pi * 5 **2, # the "5" is the radius "rho"
    "randomize": False, # alternatively "bridson", which randomly positions the cells in the system (no longer on grid)
}

boundary = [
    {"name": "box", # outer boundary type
     "expr":"true",
     "il2":{"q":0, "R":0, "bc_type": "linear"}, # initial conditions
     },

]

"""
parameters regarding meshing and fenics. unit_length_exponent is necessary for calculation concentrations. 
-6  means the simulation is in micro meters.
"""
numeric = {
    "linear_solver": "gmres", # if the uptake function is chose as linear, which solver is used
    "preconditioner": "hypre_amg", # preconditioner, check fenics documentations
    "linear": False, # if a linear upake function is chosen, choose the linear solver for better performance (--> True)
    "krylov_atol": 1e-35, # tolerances for the different solvers used, check fenics documentation for details
    "krylov_rtol": 1e-5,
    "newton_atol": 1e-35,
    "newton_rtol": 1e-5,
    "dofs_per_node": 13000, # how many dofs to use per node/thread/cpu core. If you have a mesh with 100000 dofs and set this to 10000, the sytem will use 100000/10000 = 10 cpu cores for this simulation
    "max_mpi_nodes": int(os.cpu_count()/4), # maximum amount of threads available for this run. if set to "os.cpu_count" it can use all of the available threads. If lower than "dofs_per_node" it will use the lower number
    "min_char_length": 0.07,  # mesh, smaller = finer
    "max_char_length": 7,  # mesh, smaller = finer
    "unit_length_exponent": -6  # for concentration conversion
}

IMGPATH = path + "images/" # where to save the images if you use this specific feature. Typically turned off

