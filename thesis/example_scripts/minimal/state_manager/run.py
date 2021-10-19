import logging
import os

from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.main.StateManager import StateManager
from thesis.scenarios.box_grid import setup, assign_fractions

"""Define paths for result folder structure"""
solution_path = "/extra/kiwitz/statemanager_example/test_1/"
os.makedirs(solution_path, exist_ok=True)
logging.basicConfig(
    filename=os.path.join(solution_path, "sim.log"),
    level=logging.INFO,
    filemode="w",
    format='%(levelname)s::%(asctime)s %(message)s',
    datefmt='%I:%M:%S')

scenario = setup(
    cytokines,
    cell_types_dict,
    boundary,
    geometry,
    numeric)

stMan = StateManager(solution_path)
stMan.scenario = scenario

scenario.markers = ["type_name", "IL-2_surf_c", "IL-2_R"]

stMan.T = [0, 1]


def pre_scan(state_manager, scan_index):
    assign_fractions(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan
stMan.run()
