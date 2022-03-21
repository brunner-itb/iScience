import logging

import numpy as np
from scipy.constants import N_A

from thesis.main.ParameterSet import MiscParameter
from thesis.main.my_debug import message

module_logger = logging.getLogger(__name__)


class updateState():

    def __init__(self, replicat_index, t, geometry, parameter_pool, Tsec_distribution_array, Th_distribution_array, Treg_distribution_array = [], offset=0):
        self.Tsec_distribution_array = Tsec_distribution_array
        self.Th_distribution_array = Th_distribution_array
        self.Treg_distribution_array = Treg_distribution_array
        self.geometry = geometry
        self.offset = offset
        self.parameter_pool = parameter_pool
        self.t = t
        self.drawn_fraction = 0
        self.replicat_index = replicat_index

    def get_offset_ids(self,sc):
        no_of_cells = len(sc.entity_list)
        offset = self.offset
        if offset > 0:
            try:
                a = self.geometry["z_grid"]
                dims = 3
                del a
            except KeyError:
                dims = 2
            message("dims = " + str(dims), module_logger)
            if dims == 3:
                cube_size = round(np.cbrt(no_of_cells), 0)
            else:
                cube_size = round(np.sqrt(no_of_cells), 0)
            xr = np.arange(0, cube_size, 1)
            yr = np.arange(0, cube_size, 1)

            if dims == 3:
                zr = np.arange(0, cube_size, 1)
                z_offset = np.array([zr[:offset], zr[-offset:]]).flatten()
            else:
                zr = [None]
                z_offset = []

            x_offset = np.array([xr[:offset], xr[-offset:]]).flatten()
            y_offset = np.array([yr[:offset], yr[-offset:]]).flatten()

            anti_draws = []
            counter = 0
            for x in xr:
                for y in yr:
                    for z in zr:
                        if x in x_offset or y in y_offset or z in z_offset:
                            anti_draws.append(counter)
                        counter += 1
        else:
            anti_draws = []
        return anti_draws


    def get_draws(self, sc, offset_ids, excluded_ids, fraction, at_least_one = False):
        possible_draws = np.setdiff1d(range(len(sc.entity_list)), offset_ids)
        free_draws = np.setdiff1d(possible_draws, excluded_ids)

        if fraction == 0 or len(free_draws) == 0:
            adj_fraction = 0
        else:
            adj_fraction = 1 / (len(free_draws) / (fraction * len(possible_draws)))
        # if np.abs(1.0 - adj_fraction) < 0.009:
        #     adj_fraction = 1

        amount_of_draws = int(round(len(free_draws) * adj_fraction))
        if at_least_one == True:
            if amount_of_draws == 0:
                amount_of_draws = 1
                print("set to 1")
        try:
            draws = np.random.choice(free_draws, amount_of_draws, replace=False)
        except ValueError:
            if amount_of_draws > free_draws:
                amount_of_draws = free_draws
                message("Tried to pick more than available, set to maximum. Initial fraction was " + str(fraction), module_logger)
            draws = np.random.choice(free_draws, amount_of_draws, replace=False)
        return draws

    @staticmethod
    def set_parameters(self, sc, Tsec_draws, Th_draws, Treg_draws = [], q_il2_sum = None): # fallback function, mostly not used
        t_R_start = self.parameter_pool.get_template("R_start")
        t_pSTAT5 = self.parameter_pool.get_template("pSTAT5")
        t_EC50 = self.parameter_pool.get_template("EC50")
        #calculate each Th's receptors to maintain a constant systemic R
        no_of_cells = len(Tsec_draws) + len(Th_draws) + len(Treg_draws)
        systemic_R = int(round(no_of_cells * 0.1)) * 100 + (no_of_cells - int(round(no_of_cells * 0.1))) * 5e3
        Th_R = systemic_R - len(Tsec_draws) * 100

        for i, e in enumerate(sc.entity_list):
            e.p.add_parameter_with_collection(t_pSTAT5(None, in_sim=False))
            e.p.add_parameter_with_collection(t_EC50(0, in_sim=False))
            if len(Tsec_draws) != 0 and e.type_name == "Tsec":
                if q_il2_sum != None:
                    e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum / len(Tsec_draws))
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
            elif e.type_name == "Th":
                e.p.add_parameter_with_collection(t_R_start(Th_R/len(Th_draws), in_sim=False))
                e.p.get_physical_parameter("R", "IL-2").set_in_post_unit(Th_R/len(Th_draws))
                # print(Th_R/len(Th_draws))
            elif e.type_name == "Treg":
                e.p.add_parameter_with_collection(t_R_start(5000, in_sim=False))
            elif e.type_name == "Tnaive":
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
            elif e.type_name == "blank":
                e.p.add_parameter_with_collection(t_R_start(0, in_sim=False))

    @staticmethod
    def set_R_lognorm_parameters(self, sc, Tsec_draws, Th_draws, Treg_draws, q_il2_sum):
        t_R_start = self.parameter_pool.get_template("R_start")
        t_EC50 = self.parameter_pool.get_template("EC50")
        t_pSTAT5 = self.parameter_pool.get_template("pSTAT5")
        no_of_cells = len(Tsec_draws) + len(Th_draws) + len(Treg_draws)
        systemic_R = int(round(no_of_cells * 0.1)) * 100 + (no_of_cells - int(round(no_of_cells * 0.1))) * 1.5e3
        Th_R = systemic_R - len(Tsec_draws) * 100

        for i, e in enumerate(sc.entity_list):
            e.p.add_parameter_with_collection(t_EC50(0, in_sim=False))
            e.p.add_parameter_with_collection(t_pSTAT5(0, in_sim=False))
            try:
                gamma = e.p.get_physical_parameter("gamma", "IL-2").get_in_sim_unit()
            except:
                gamma = float(e.p.get_as_dictionary()["misc_gamma"])
            if len(Tsec_draws) != 0 and e.type_name == "Tsec":
                if q_il2_sum != None:
                    e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum / len(Tsec_draws))
                Tsec_R = e.p.get_physical_parameter("R", "IL-2").get_in_post_unit()
                e.p.add_parameter_with_collection(t_R_start(Tsec_R, in_sim=False))
                E = 0

            elif e.type_name == "Th":
                if gamma < 1:
                    Th_start_R = e.p.get_misc_parameter("R_start_neg", "misc").get_in_post_unit()
                elif gamma > 1:
                    Th_start_R = e.p.get_misc_parameter("R_start_pos", "misc").get_in_post_unit()
                elif gamma == 1:
                    if e.p.get_as_dictionary()["scan_name_scan_name"] == "R_scan":
                        Th_start_R = e.p.get_as_dictionary()["scan_value"]
                    else:
                        Th_start_R = Th_R/len(Th_draws)
                # print(Th_start_R)
                e.p.add_parameter_with_collection(t_R_start(Th_start_R, in_sim=False))
                E = Th_start_R * N_A ** -1 * 1e9

            elif e.type_name == "Treg":
                if gamma < 1:
                    Treg_start_R = e.p.get_misc_parameter("R_start_neg", "misc").get_in_post_unit()
                elif gamma > 1:
                    Treg_start_R = e.p.get_misc_parameter("R_start_pos", "misc").get_in_post_unit()
                elif gamma == 1:
                    if e.p.get_as_dictionary()["scan_name_scan_name"] == "R_scan":
                        Treg_start_R = e.p.get_physical_parameter("R", "IL-2").get_in_post_unit()
                    else:
                        Treg_start_R = 5001
                e.p.add_parameter_with_collection(t_R_start(Treg_start_R, in_sim=False))
                E = Treg_start_R * N_A ** -1 * 1e9

            elif e.type_name == "Tnaive":
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
                E = 0

            elif e.type_name == "blank":
                e.p.add_parameter_with_collection(t_R_start(0, in_sim=False))

            if E != 0:
                if e.p.get_as_dictionary()["scan_name_scan_name"] == "sigma":
                    var = e.p.get_physical_parameter("sigma", "IL-2").get_in_post_unit()
                else:
                    var = e.p.get_misc_parameter("sigma", "misc").get_in_post_unit()
                tmp_sigma = np.sqrt(np.log((var * E) ** 2 / E ** 2 + 1))
                mean = np.log(E) - 1 / 2 * tmp_sigma ** 2
                R_draw = np.random.lognormal(mean, tmp_sigma)

                e.p.get_physical_parameter("R", "IL-2").set_in_sim_unit(R_draw)
                e.p.add_parameter_with_collection(t_R_start(R_draw, in_sim=True))


    def set_cell_types(self,sc):
        offset_ids = self.get_offset_ids(sc)
        varying_Tsec_fraction = False
        try:
            Tsec_fraction = sc.entity_list[0].p.get_physical_parameter("Tsec_fraction", "IL-2").get_in_sim_unit()
            varying_Tsec_fraction = True
        except:
            Tsec_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Tsec"]

        try:
            Treg_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Treg"]
            Treg_fraction = sc.entity_list[0].p.get_physical_parameter("Treg_fraction", "IL-2").get_in_sim_unit()
        except AttributeError:
            Treg_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Treg"]
        except KeyError:
            Treg_fraction = 0


        Tsec_draws = self.get_draws(sc, offset_ids = offset_ids, excluded_ids=[], fraction=Tsec_fraction)


        if len(self.Tsec_distribution_array) != len(Tsec_draws):
            self.Tsec_distribution_array = Tsec_draws
        else:
            Tsec_draws = self.Tsec_distribution_array
        if varying_Tsec_fraction == True or Treg_fraction !=0:
            Th_fraction = 1 - Tsec_fraction - Treg_fraction
        elif Treg_fraction == 0:
            Th_fraction = sc.entity_list[0].p.get_as_dictionary()["fractions_Th"]

        Treg_draws = self.get_draws(sc, offset_ids = offset_ids, excluded_ids=list(Tsec_draws), fraction=Treg_fraction)
        if len(self.Treg_distribution_array) != len(Treg_draws):
            self.Treg_distribution_array = Treg_draws
        else:
            Treg_draws = self.Treg_distribution_array

        Th_draws = self.get_draws(sc, offset_ids = offset_ids, excluded_ids=list(Tsec_draws) + list(Treg_draws), fraction=Th_fraction)
        if len(self.Th_distribution_array) != len(Th_draws):
            self.Th_distribution_array = Th_draws
        else:
            Th_draws = self.Th_distribution_array

        no_of_Th = 0
        no_of_Treg = 0

        for i, e in enumerate(sc.entity_list):
            e.change_type = "Th"
            if i in Tsec_draws:
                e.change_type = "Tsec"
            if i in offset_ids:
                e.change_type = "Th"
                no_of_Th += 1
            if i in Th_draws:
                e.change_type = "Th"
                no_of_Th += 1
            if i in Treg_draws:
                e.change_type = "Treg"
                no_of_Treg += 1
            e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
        sc.apply_type_changes(self.replicat_index)

        message("Number of secreting cells: " + str(len(Tsec_draws)), module_logger)
        message("Number of Ths: " + str(no_of_Th), module_logger)
        message("Number of Tregs: " + str(no_of_Treg), module_logger)

        try:
            global_q = sc.entity_list[Tsec_draws[0]].p.get_physical_parameter("global_q", "IL-2").get_in_sim_unit()
        except:
            global_q = sc.entity_list[Tsec_draws[0]].p.get_as_dictionary()["IL-2_global_q"]
        q = sc.entity_list[Tsec_draws[0]].p.get_physical_parameter("q", "IL-2").get_in_post_unit()
        if global_q == True:
            q_il2_sum = len(sc.entity_list) * 0.1 * q * N_A ** -1 * 1e9
            if len(Tsec_draws) != 0:
                message("Cells q: " + str(q_il2_sum / len(Tsec_draws) / (N_A ** -1 * 1e9)), module_logger)
            else:
                message("Cells q = 0", module_logger)
        else:
            for draw in Tsec_draws:
                half_time = 1
                lamb = np.log(2) / half_time
                beta = 1 / lamb
                try:
                    sc.entity_list[draw].p.get_misc_parameter("sec_start", "misc").set_in_post_unit(np.random.exponential(beta, 1) * 3600)
                except AttributeError:
                    pass
            if len(Tsec_draws) != 0:
                message("Cells q: " + str(q), module_logger)
            else:
                message("Cells q: 0", module_logger)
            q_il2_sum = None
        return Tsec_draws, Th_draws, Treg_draws, q_il2_sum


    def step(self,sc):
        Tsec_draws, Th_draws, Treg_draws, q_il2_sum = self.set_cell_types(sc)
        if sc.entity_list[0].p.get_as_dictionary()["scan_name_scan_name"] == "sigma":
            var = sc.entity_list[0].p.get_as_dictionary()["scan_value"]
        else:
            var = sc.entity_list[0].p.get_misc_parameter("sigma", "misc").get_in_post_unit()
        self.set_R_lognorm_parameters(self, sc, Tsec_draws, Th_draws, Treg_draws, q_il2_sum)


    def reset_draw_arrays(self):
        self.Tsec_distribution_array = np.array([])
        self.Th_distribution_array = np.array([])




    def step_close_to_ODE(self,sc):
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        print("special step, only if f_Tsec >= 0.5 !!!!!")
        Tsec_draws, Th_draws, Treg_draws, q_il2_sum = self.set_cell_types(sc)
        self.set_parameters_close_to_ODE(sc, Tsec_draws, Th_draws, Treg_draws, q_il2_sum)

    def set_parameters_close_to_ODE(self, sc, Tsec_draws, Th_draws, Treg_draws, q_il2_sum):
        t_R_start = self.parameter_pool.get_template("R_start")
        t_EC50 = self.parameter_pool.get_template("EC50")
        no_of_cells = len(Tsec_draws) + len(Th_draws) + len(Treg_draws)
        systemic_R = int(round(no_of_cells * 0.1)) * 100 + (no_of_cells - int(round(no_of_cells * 0.1))) * 1.5e3
        Th_R = (systemic_R - no_of_cells * 0.5 * 100)/(no_of_cells * 0.5)
        systemic_Th_R = Th_R * len(Th_draws)
        Tsec_R = (systemic_R - systemic_Th_R) / len(Tsec_draws)
        message("Tsec_R: " + str(Tsec_R))
        message("Th_R: " + str(Th_R))

        for i, e in enumerate(sc.entity_list):
            if len(Tsec_draws) != 0 and e.type_name == "Tsec":
                if q_il2_sum != None:
                    e.p.get_physical_parameter("q", "IL-2").set_in_sim_unit(q_il2_sum / len(Tsec_draws))
                e.p.add_parameter_with_collection(t_R_start(Tsec_R, in_sim=False))
                E = Tsec_R * N_A ** -1 * 1e9

            elif e.type_name == "Th":

                if e.p.get_as_dictionary()["scan_name_scan_name"] == "R_scan":
                    Th_R = e.p.get_as_dictionary()["scan_value"]
                # print(Th_start_R)
                e.p.add_parameter_with_collection(t_R_start(Th_R, in_sim=False))
                E = Th_R * N_A ** -1 * 1e9

            elif e.type_name == "Tnaive":
                e.p.add_parameter_with_collection(t_R_start(1e2, in_sim=False))
                E = 0

            elif e.type_name == "blank":
                e.p.add_parameter_with_collection(t_R_start(0, in_sim=False))

            if E != 0:
                if e.p.get_as_dictionary()["scan_name_scan_name"] == "sigma":
                    var = e.p.get_physical_parameter("sigma", "IL-2").get_in_post_unit()
                else:
                    var = e.p.get_misc_parameter("sigma", "misc").get_in_post_unit()
                tmp_sigma = np.sqrt(np.log((var * E) ** 2 / E ** 2 + 1))
                mean = np.log(E) - 1 / 2 * tmp_sigma ** 2
                R_draw = np.random.lognormal(mean, tmp_sigma)

                e.p.get_physical_parameter("R", "IL-2").set_in_sim_unit(R_draw)
                e.p.add_parameter_with_collection(t_R_start(R_draw, in_sim=True))

