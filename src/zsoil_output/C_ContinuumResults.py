"""Read and post-process results for two- and three-dimensional continuum elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Continuum_EleResults` reader supports triangular, quadrilateral, brick,
tetrahedral, and wedge formulations. It also provides helpers for transforming stresses and
integrating continuum stresses into equivalent beam resultants.
"""

import os
import numpy as np

from .C_EleResults import *
from .C_Element_cont_Q4 import *
from .C_Element_cont_T3 import *
from .C_Element_cont_B8 import *
from .C_Element_cont_TH4 import *
from .C_Element_cont_W6 import *
from .C_NodalResults import *

#from C_PlasticCodes       import *

from .C_math_utl import *


class Continuum_EleResults (EleResults):
    SCALAR_TYPE = 0
    VECTOR_TYPE = 1
    TENSOR_TYPE = 2

    comps_tns = ['XX', 'YY', 'XY', 'ZZ', 'XZ', 'YZ']
    comps_tns_normal = ['XX', 'YY', 'ZZ']
    comps_vec = ['X', 'Y', 'Z']
    comps_sca = ['']

    # use these keys in functions for rsl_key item
    STRESS_KEY = "STRESSES"
    STRAIN_KEY = "STRAINS"
    PLCODE_KEY = "PLA_CODE"
    STRLEVEL_KEY = "STR_LEVEL"
    SATURATION_KEY = "SATUR"
    FLUID_VEL_KEY = "FLU_VELOC"
    HARD_PC_KEY = "HARD_PC"
    HARD_GAMMA_KEY = "HARD_GAMMA"
    YOUNG_ET_KEY = "YOUNG_ET"
    YOUNG_E0_KEY = "YOUNG_E0"
    DAMAGE_KEY = "DAMAGE"
    UNDR_PRESS_KEY = "UNDR_PRESS"
    EPSVP_EQ_KEY = "EPSVP_EQ"
    HEAT_FLUX_KEY = "HEAT_FLUX"
    HUMI_FLUX_KEY = "HEAT_FLUX"
    FT_KEY = "FT"
    MATURITY_KEY = "MATURITY"
    TEMPERATURE_KEY = "TEMP"
    GRADH_KEY = "GRAD_H"

    # this is an extra key as it combines few results
    TSTRESS_KEY_EXTRA = "TSTRESSES"


# TODO
    rsl_types = {STRESS_KEY: TENSOR_TYPE,
                 STRAIN_KEY: TENSOR_TYPE,
                 PLCODE_KEY: SCALAR_TYPE,
                 STRLEVEL_KEY: SCALAR_TYPE,
                 SATURATION_KEY: SCALAR_TYPE,
                 FLUID_VEL_KEY: VECTOR_TYPE,
                 HARD_PC_KEY: SCALAR_TYPE,
                 HARD_GAMMA_KEY: SCALAR_TYPE,
                 YOUNG_ET_KEY: SCALAR_TYPE,
                 YOUNG_E0_KEY: SCALAR_TYPE,
                 UNDR_PRESS_KEY: SCALAR_TYPE,
                 EPSVP_EQ_KEY: SCALAR_TYPE,
                 HEAT_FLUX_KEY: VECTOR_TYPE,
                 HUMI_FLUX_KEY: VECTOR_TYPE,
                 FT_KEY: SCALAR_TYPE,
                 MATURITY_KEY: SCALAR_TYPE,
                 TEMPERATURE_KEY: SCALAR_TYPE,
                 GRADH_KEY: VECTOR_TYPE
                 }

    def __init__(self, mesh, his, rcf, nodal_results_instance=None):
        """Initialize the continuum element results object for the supplied project data.

        :param mesh: Mesh containing the model nodes, elements, and materials.
        :type mesh: Mesh
        :param his: Execution history that maps solution indices to times and states.
        :type his: HistoryOfExecution
        :param rcf: Result configuration describing stored nodal and element fields.
        :type rcf: RCF_info
        :param nodal_results_instance: Existing nodal-results reader to reuse, if available.
        :type nodal_results_instance: Optional[NodalResults]
        """
        EleResults.__init__(self, mesh, his, rcf)
        self._nodal_results = nodal_results_instance

    def _get_rsl_time_history(self, ele_index_1, solution_indices, rsl_key, comp=None):
        """Return the results history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to select; None requests the complete vector.
        :type comp: Optional[str]
        :return: The results history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        element = self.mesh.get_element(ele_index_1)
        ngaus = len(element.xsiGP)
        if not isinstance(element, Element_cont_Q4) and not isinstance(element, Element_cont_T3) \
           and not isinstance(element, Element_cont_B8) and not isinstance(element, Element_cont_TH4) \
           and not isinstance(element, Element_cont_W6):
            return None

        sel_elements = [self.mesh.get_element(ele_index_1)]
        size_times = len(solution_indices)

        if self.rsl_types[rsl_key] == self.TENSOR_TYPE:
            comps_aux = self.comps_tns
        elif self.rsl_types[rsl_key] == self.VECTOR_TYPE:
            comps_aux = self.comps_vec
        elif self.rsl_types[rsl_key] == self.SCALAR_TYPE:
            comps_aux = [""]

        cols = []
        for comp_ex in comps_aux:
            col, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, comp_ex)
            if col >= 0:
                cols.append(col)

        if comp == None:
            size = len(cols)
            if size == 1:
                out = np.zeros(size_times * ngaus).reshape(size_times, ngaus)
            else:
                out = np.zeros(size_times * ngaus * size).reshape(size_times, ngaus, size)
        else:
            col, ncomps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, comp)
            #cols.append (col)
            # BUG 22.07.2019
            cols = [col]
            out = np.zeros(size_times*ngaus).reshape(size_times, ngaus)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        res_all = self.get_element_results_ex(sel_elements, solution_indices)

        for igaus in range(ngaus):
            res = res_all[igaus]
            offs = 2  # skip element index and gauss point index
            for i in range(size_times):
                for k, col in enumerate(cols):
                    if comp == None and len(cols) > 1:
                        out[i, igaus, k] = res[offs + col]
                    else:
                        out[i, igaus] = res[offs + col]
                offs = offs + gp_size_int
        return out

    def _get_rsl_for_sel_elements(self, ele_indices_1, solution_index, rsl_key, comp=None):
        """Return results for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to select; None requests the complete vector.
        :type comp: Optional[str]
        :return: The results for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        element = self.mesh.get_element(ele_indices_1[0])
        ngaus = len(element.xsiGP)

        solution_indices = [solution_index]
        size = len(ele_indices_1)

        if self.rsl_types[rsl_key] == self.TENSOR_TYPE:
            if comp == None:
                out = np.zeros(size * ngaus * 6).reshape(size, ngaus, 6)
            else:
                out = np.zeros(size * ngaus).reshape(size, ngaus)

        elif self.rsl_types[rsl_key] == self.VECTOR_TYPE:

            if comp == None:
                out = np.zeros(size * ngaus * 3).reshape(size, ngaus, 3)
            else:
                out = np.zeros(size * ngaus).reshape(size, ngaus)

        elif self.rsl_types[rsl_key] == self.SCALAR_TYPE:

            if comp == None:
                out = np.zeros(size * ngaus).reshape(size, ngaus)
            else:
                out = np.zeros(size * ngaus).reshape(size, ngaus)

        for i, ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history(ele_index_1, solution_indices, rsl_key, comp)
            for igaus in range(ret.shape[1]):
                if ret.ndim == 2:
                    out[i][igaus] = ret[0][igaus]
                else:
                    for k in range(ret.shape[2]):  # loop over components
                        out[i][igaus][k] = ret[0][igaus][k]
        return out

    def get_rsl_time_history(self, ele_index_1, solution_indices, rsl_key):
        """Return the results history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :return: NumPy array of shape ``(n_time_instances, ngaus, size)`` for tensor/vector result or NumPy array of
            shape ``(n_time_instances, ngaus)`` for scalar.
        :rtype: Optional[numpy.ndarray]
        """
        rsl_key_ex = rsl_key
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            rsl_key_ex = self.STRESS_KEY
        out = self._get_rsl_time_history(ele_index_1, solution_indices, rsl_key_ex, None)
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            for i, solution_index in enumerate(solution_indices):
                self.make_tot_stress(ele_index_1, solution_index, out[i])
        return out

    def get_rsl_time_history_ex(self, ele_index_1, solution_indices, rsl_key, comp):
        """Return the results history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to return.
        :type comp: str
        :return: NumPy array of shape ``(n_time_instances, ngaus)``.
        :rtype: Optional[numpy.ndarray]
        """
        rsl_key_ex = rsl_key
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            rsl_key_ex = self.STRESS_KEY
        out = self._get_rsl_time_history(ele_index_1, solution_indices, rsl_key_ex, comp)
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            for i, solution_index in enumerate(solution_indices):
                self.make_tot_stress(ele_index_1, solution_index, out[i], comp)
        return out

    def get_rsl_for_sel_elements(self, ele_indices_1, solution_index, rsl_key):
        """Return results for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :return: NumPy array of shape ``(n_elements, ngaus, size)`` for tensor/vector result or NumPy array of shape
            ``(n_elements, ngaus)`` for scalar result.
        :rtype: numpy.ndarray
        """
        rsl_key_ex = rsl_key
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            rsl_key_ex = self.STRESS_KEY
        out = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, rsl_key_ex, None)
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            for i, ele_index_1 in enumerate(ele_indices_1):
                self.make_tot_stress(ele_index_1, solution_index, out[i])
        return out

    def get_rsl_for_sel_elements_ex(self, ele_indices_1, solution_index, rsl_key, comp):
        """Return results for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to return.
        :type comp: str
        :return: NumPy array of shape ``(n_elements, ngaus)``.
        :rtype: numpy.ndarray
        """
        rsl_key_ex = rsl_key
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            rsl_key_ex = self.STRESS_KEY
        out = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, rsl_key_ex, comp)
        if rsl_key == self.TSTRESS_KEY_EXTRA:
            for i, ele_index_1 in enumerate(ele_indices_1):
                self.make_tot_stress(ele_index_1, solution_index, out[i], comp)
        return out

    def make_tot_stress(self, ele_index_1, solution_index, out, comp=None):
        """Convert effective stresses to total stresses for the requested element and solution step.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param out: Out used by this operation.
        :type out: numpy.ndarray
        :param comp: Component label to select; None requests the complete vector.
        :type comp: Optional[str]
        """
        if comp != None and not comp in self.comps_tns_normal:
            return

        if self._nodal_results == None:
            print("NodalResults instance must be given to the constructor of ContinuumEleResults class")
            print("Total stress is not computed")

        # input is [1, ngaus] if comp = 'XX' or any other
        # input is [1, ngaus,size] if comp = None

        # first take nodal pressures if they exist
        element = self.mesh.get_element(ele_index_1)
        p = self._nodal_results.get_nodes_pore_pressures(element.nodes, solution_index)

        # then take saturation
        #output is [1, ngaus]
        ele_indices_1 = [element.index]
        S = self.get_rsl_for_sel_elements(ele_indices_1, solution_index, self.SATURATION_KEY)

        # then take undrained pressure
        # output is [1, ngaus]
        dp = self.get_rsl_for_sel_elements(ele_indices_1, solution_index, self.UNDR_PRESS_KEY)

        # TODO only center so far
        po = 0.0
        for i in range(len(p)):
            po = po + p[i]
        po = po / len(p)

        dpo = dp[0, 0]

        So = S[0, 0]

        if comp == None:
            out[0, 0] = out[0, 0] + So * po + dpo
            out[0, 1] = out[0, 1] + So * po + dpo
            out[0, 3] = out[0, 3] + So * po + dpo
        else:
            pos = self.comps_tns.index(comp)
            out[0] = out[0] + So * po + dpo

    def get_stress_transform_mtrx(self, T_GL):
        """Return the matrix that transforms stress components between coordinate systems.

        :param T_GL: T gl used by this operation.
        :type T_GL: numpy.ndarray
        :return: The matrix that transforms stress components between coordinate systems.
        :rtype: numpy.ndarray
        """
        T = np.zeros(36).reshape(6, 6)

        T[0, 0] = T_GL[0, 0] * T_GL[0, 0]
        T[0, 1] = T_GL[0, 1] * T_GL[0, 1]
        T[0, 2] = 2.0 * T_GL[0, 0] * T_GL[0, 1]
        T[0, 3] = T_GL[0, 2] * T_GL[0, 2]
        T[0, 4] = 2.0 * T_GL[0, 0] * T_GL[0, 2]
        T[0, 5] = 2.0 * T_GL[0, 1] * T_GL[0, 2]

        T[1, 0] = T_GL[1, 0] * T_GL[1, 0]
        T[1, 1] = T_GL[1, 1] * T_GL[1, 1]
        T[1, 2] = 2.0 * T_GL[1, 0] * T_GL[1, 1]
        T[1, 3] = T_GL[1, 2] * T_GL[1, 2]
        T[1, 4] = 2.0 * T_GL[1, 0] * T_GL[1, 2]
        T[1, 5] = 2.0 * T_GL[1, 1] * T_GL[1, 2]

        T[2, 0] = T_GL[0, 0] * T_GL[1, 0]
        T[2, 1] = T_GL[0, 1] * T_GL[1, 1]
        T[2, 2] = T_GL[0, 0] * T_GL[1, 1] + T_GL[0, 1] * T_GL[1, 0]
        T[2, 3] = T_GL[0, 2] * T_GL[1, 2]
        T[2, 4] = T_GL[0, 0] * T_GL[1, 2] + T_GL[0, 2] * T_GL[1, 0]
        T[2, 5] = T_GL[0, 1] * T_GL[1, 2] + T_GL[0, 2] * T_GL[1, 1]

        T[3, 0] = T_GL[2, 0] * T_GL[2, 0]
        T[3, 1] = T_GL[2, 1] * T_GL[2, 1]
        T[3, 2] = 2.0 * T_GL[2, 0] * T_GL[2, 1]
        T[3, 3] = T_GL[2, 2] * T_GL[2, 2]
        T[3, 4] = 2.0 * T_GL[2, 0] * T_GL[2, 2]
        T[3, 5] = 2.0 * T_GL[2, 1] * T_GL[2, 2]

        T[4, 0] = T_GL[0, 0] * T_GL[2, 0]
        T[4, 1] = T_GL[0, 1] * T_GL[2, 1]
        T[4, 2] = T_GL[0, 0] * T_GL[2, 1] + T_GL[0, 1] * T_GL[2, 0]
        T[4, 3] = T_GL[0, 2] * T_GL[2, 2]
        T[4, 4] = T_GL[0, 0] * T_GL[2, 2] + T_GL[0, 2] * T_GL[2, 0]
        T[4, 5] = T_GL[0, 1] * T_GL[2, 2] + T_GL[0, 2] * T_GL[2, 1]

        T[5, 0] = T_GL[1, 0] * T_GL[2, 0]
        T[5, 1] = T_GL[1, 1] * T_GL[2, 1]
        T[5, 2] = T_GL[1, 0] * T_GL[2, 1] + T_GL[1, 1] * T_GL[2, 0]
        T[5, 3] = T_GL[1, 2] * T_GL[2, 2]
        T[5, 4] = T_GL[1, 0] * T_GL[2, 2] + T_GL[1, 2] * T_GL[2, 0]
        T[5, 5] = T_GL[1, 1] * T_GL[2, 2] + T_GL[1, 2] * T_GL[2, 1]

        return T

    def get_MNQ_time_history_for_fict_uniform_beam(self, sel_elements, solution_indices,
                                                   x_axis_pt, x_axis_vec, y_axis_vec):
        """Return the axial, shear, and bending resultants history across the requested solution steps.

        :param sel_elements: Elements to process.
        :type sel_elements: Sequence[Element]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param x_axis_pt: X axis pt used by this operation.
        :type x_axis_pt: Sequence[float]
        :param x_axis_vec: X axis vector used by this operation.
        :type x_axis_vec: Sequence[float]
        :param y_axis_vec: Y axis vector used by this operation.
        :type y_axis_vec: Sequence[float]
        :return: Status flag followed by coordinate, section-size, and force-resultant data; failure paths contain empty
            placeholders.
        :rtype: tuple
        """
        def swap(a, b):
            """Return the two supplied values in reversed order.

            :param a: A used by this operation.
            :type a: array_like
            :param b: B used by this operation.
            :type b: array_like
            """
            tmp = a
            a = b
            b = tmp

        if len(sel_elements) == 0:
            return False, [], [], []

        pt = np.array(x_axis_pt)
        versor = np.array(x_axis_vec)
        length = np.linalg.norm(versor)
        versor = versor / length

        sorted_sel_ele, proj_measures = self.mesh.sort_sel_elements_by_dist_along_dir(sel_elements,
                                                                                      x_axis_pt, x_axis_vec, True)
        indices_of_sorted_sel_ele = self.mesh.get_ele_indices_for_sel_elements(sorted_sel_ele)

        ele_surf_in_dir_X = np.zeros(len(indices_of_sorted_sel_ele))
        ele_size_in_dir_X = np.zeros(len(indices_of_sorted_sel_ele))

        T_GL = np.zeros(9).reshape(3, 3)

        length = np.linalg.norm(x_axis_vec)
        e_x = x_axis_vec / length
        length = np.linalg.norm(y_axis_vec)
        e_y = y_axis_vec / length
        e_z = np.cross(e_x, e_y)

        T_GL[0, :] = e_x
        T_GL[1, :] = e_y
        T_GL[2, :] = e_z

        T_sig = self.get_stress_transform_mtrx(T_GL)

        for i, element in enumerate(sorted_sel_ele):
            ele_surf_in_dir_X[i] = element.get_surface_along_dir(e_x)
            ele_size_in_dir_X[i] = element.get_size_along_dir(e_x)

        # find largest projection measure distance
        if len(proj_measures) == 1:
            n_ele_in_perp_dir = 1
            n_ele_along_axis = 1
        else:
            last_proj_meas = proj_measures[0]
            max_proj_meas = 0.0
            for i in range(len(proj_measures)-1):
                next_proj_meas = proj_measures[i+1]
                max_proj_meas = max(max_proj_meas, abs(next_proj_meas-last_proj_meas))
                last_proj_meas = next_proj_meas

            # detect how many elements are present in direction perpendicular to the axis
            check = [1]
            last_proj_meas = proj_measures[0]

            for i in range(len(proj_measures)-1):
                next_proj_meas = proj_measures[i+1]
                dist = abs(next_proj_meas-last_proj_meas)
                if dist/max_proj_meas < 1.0e-2:
                    check[-1] = check[-1] + 1
                else:
                    check.append(1)
                last_proj_meas = next_proj_meas

            for i in range(len(check)-1):
                if check[i+1] - check[i] != 0:
                    return False, [], [], [], [], []

            n_ele_in_perp_dir = check[0]
            n_ele_along_axis = len(check)

        # generalized forces as in beams [NX,QY,QZ,MX,MY,MZ]
        Forces = np.zeros(n_ele_along_axis * len(solution_indices) *
                          6).reshape(n_ele_along_axis, len(solution_indices), 6)

        x_tab = np.zeros(n_ele_along_axis)
        x_size = np.zeros(n_ele_along_axis)

        sig_aux = np.zeros(6)
        axis_pt = np.zeros(3)

        count = 0
        for j in range(n_ele_along_axis):
            x_tab[j] = proj_measures[count]
            x_size[j] = ele_size_in_dir_X[count]
            count = count + n_ele_in_perp_dir

        for i, solution_index in enumerate(solution_indices):
            # these results are for all elements
            sig = self.get_rsl_for_sel_elements(indices_of_sorted_sel_ele,
                                                solution_index, self.STRESS_KEY)
            ngaus = sig.shape[1]
            # TODO 1 GP so far (we have that in ZSOIL anyway)
            ngaus = 1
            count = 0
            for j in range(n_ele_along_axis):
                axis_pt[:] = pt[:] + versor[:] * proj_measures[count]
                for k in range(n_ele_in_perp_dir):
                    area = ele_surf_in_dir_X[count] / float(ngaus)
                    # TODO this ok for 1GP results
                    ele_center = sorted_sel_ele[count].get_center()
                    dX = ele_center - axis_pt
                    dx = T_GL.dot(dX)

                    for igaus in range(ngaus):
                        sig_aux[:] = sig[count, igaus]
                        sig_tr = T_sig.dot(sig_aux)
                        # axial force
                        jaja = Forces[j, i, 0]
                        Forces[j, i, 0] = Forces[j, i, 0] + sig_tr[0] * area
                        # shear force Qy
                        Forces[j, i, 1] = Forces[j, i, 1] + sig_tr[2] * area
                        # shear force Qz
                        Forces[j, i, 2] = Forces[j, i, 2] + sig_tr[4] * area
                        # torsion moment MX
                        Forces[j, i, 3] = Forces[j, i, 3] + \
                            (sig_tr[2] * dx[2] - sig_tr[4] * dx[1]) * area
                        # bending moment MY
                        Forces[j, i, 4] = Forces[j, i, 4] + sig_tr[0] * area * dx[2]
                        # bending moment MZ
                        Forces[j, i, 5] = Forces[j, i, 5] + sig_tr[0] * area * dx[1]

                    # be careful here because in beam elements Y and Z axes are interchanged

                    count = count + 1

                swap(Forces[j, i, 1], Forces[j, i, 2])  # shear forces
                swap(Forces[j, i, 4], Forces[j, i, 5])  # bending moments

        return True, x_tab, x_size, Forces

# static function


def postprocess_MNQ_time_history_for_fict_uniform_beam(Forces, x_tab, x_size):
    """Convert raw section resultants into a consistent beam time-history layout.

    Prepare the resultants for plotting and correct the bending moments using the shear forces.

    :param Forces: Forces used by this operation.
    :type Forces: numpy.ndarray
    :param x_tab: X tab used by this operation.
    :type x_tab: array_like
    :param x_size: X size used by this operation.
    :type x_size: int
    :return: Processed force-resultant array and corresponding coordinates.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    Forces_ex = np.zeros(Forces.size * 2).reshape(Forces.shape[0] * 2, Forces.shape[1], Forces.shape[2])
    x_tab_ex = np.zeros(x_tab.size * 2)

    count = 0
    for i in range(x_tab.size):
        dx2 = x_size[i] / 2.0
        x_tab_ex[count+0] = x_tab[i] - dx2
        x_tab_ex[count+1] = x_tab[i] + dx2
        count = count + 2

    # loop over time instances
    for i in range(Forces.shape[1]):
        # loop over element patches along beam axis
        count = 0
        for j in range(Forces.shape[0]):
            Forces_ex[count+0, i, 0] = Forces[j, i, 0]
            Forces_ex[count+1, i, 0] = Forces[j, i, 0]

            Forces_ex[count+0, i, 1] = Forces[j, i, 1]
            Forces_ex[count+1, i, 1] = Forces[j, i, 1]

            Forces_ex[count+0, i, 2] = Forces[j, i, 2]
            Forces_ex[count+1, i, 2] = Forces[j, i, 2]

            Forces_ex[count+0, i, 3] = Forces[j, i, 3]
            Forces_ex[count+1, i, 3] = Forces[j, i, 3]

            Forces_ex[count+0, i, 4] = Forces[j, i, 4] - Forces[j, i, 2] * x_size[j] / 2.0
            Forces_ex[count+1, i, 4] = Forces[j, i, 4] + Forces[j, i, 2] * x_size[j] / 2.0

            Forces_ex[count+0, i, 5] = Forces[j, i, 5] - Forces[j, i, 1] * x_size[j] / 2.0
            Forces_ex[count+1, i, 5] = Forces[j, i, 5] + Forces[j, i, 1] * x_size[j] / 2.0

            count = count + 2

    return x_tab_ex, Forces_ex

# static function


def get_MNQ_envelopes_for_fict_uniform_beam(Forces, comp_index_6):
    """Return the minimum and maximum axial, shear, and bending resultant values across the requested solution steps.

    :param Forces: Forces used by this operation.
    :type Forces: numpy.ndarray
    :param comp_index_6: => [NX,QY,QZ,MX,MY,MZ].
    :type comp_index_6: int
    :return: The minimum and maximum axial, shear, and bending resultant values across the requested solution steps.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    out_max = np.zeros(Forces.shape[0] * 6).reshape(Forces.shape[0], 6)
    out_min = np.zeros(Forces.shape[0] * 6).reshape(Forces.shape[0], 6)

    # loop over time instances
    for i in range(Forces.shape[1]):
        # loop over element patches along beam axis
        for j in range(Forces.shape[0]):
            if Forces[j, i, comp_index_6] > out_max[j, comp_index_6]:
                out_max[j, comp_index_6] = Forces[j, i, comp_index_6]
                for k in range(6):
                    if k != comp_index_6:
                        out_max[j, k] = Forces[j, i, k]

            if Forces[j, i, comp_index_6] < out_min[j, comp_index_6]:
                out_min[j, comp_index_6] = Forces[j, i, comp_index_6]
                for k in range(6):
                    if k != comp_index_6:
                        out_min[j, k] = Forces[j, i, k]

    return out_min, out_max


def smooth_jumps_in_stress_resultant(stress_result, comp):
    """Smooth discontinuities in one stress-resultant component.

    :param stress_result: Stress result used by this operation.
    :type stress_result: numpy.ndarray
    :param comp: Component label to return.
    :type comp: str
    """
    count = 1
    while count < stress_result.shape[0] - 1:
        stress_result[count, comp] = \
            0.5 * (stress_result[count, comp] + stress_result[count + 1, comp])
        stress_result[count + 1, comp] = stress_result[count, comp]
        count = count + 2


def smooth_jumps_in_stress_resultant_time_history(stress_result, comp):
    """Smooth discontinuities in a stress-resultant time history.

    :param stress_result: Stress result used by this operation.
    :type stress_result: numpy.ndarray
    :param comp: Component label to return.
    :type comp: str
    """
    for i in range(stress_result.shape[1]):
        smooth_jumps_in_stress_resultant(stress_result[:, i, :], comp)


def combine_stress_resultants(x_tab1, x_size1, Forces1, x_tab2, x_size2, Forces2, tol=None):
    """Combine two sets of stress resultants on a shared coordinate grid.

    :param x_tab1: X tab1 used by this operation.
    :type x_tab1: array_like
    :param x_size1: X size1 used by this operation.
    :type x_size1: int
    :param Forces1: Forces1 used by this operation.
    :type Forces1: Union[float, array_like]
    :param x_tab2: X tab2 used by this operation.
    :type x_tab2: array_like
    :param x_size2: X size2 used by this operation.
    :type x_size2: int
    :param Forces2: Forces2 used by this operation.
    :type Forces2: Union[float, array_like]
    :param tol: Tol used by this operation.
    :type tol: Optional[float]
    :return: Combined coordinates, section sizes, and force-resultant arrays.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    if tol == None:
        dx1 = abs(np.mean(x_tab1))
        dx2 = abs(np.mean(x_tab2))
        tol_ex = 1.0e-4 * min(dx1, dx2)
    else:
        tol_ex = tol

    x_tab, info = math_concatenate_arrays(x_tab1, x_tab2, tol_ex)
    #print "concatanate", x_tab1.shape,x_tab2.shape,x_tab.shape
    #print x_tab1
    #print x_tab2

    size = x_tab.size
    dim1 = Forces1.shape[1]
    dim2 = Forces1.shape[2]
    Forces = np.zeros(size * dim1 * dim2).reshape(size, dim1, dim2)
    x_size = np.zeros(size)

    for i in range(len(x_tab)):
        for pos in info[i]:
            if pos < x_tab1.size:
                Forces[i] = Forces[i] + Forces1[pos]
                x_size[i] = x_size1[pos]
            else:
                pos2 = pos - x_tab1.size
                Forces[i] = Forces[i] + Forces2[pos2]
                x_size[i] = x_size2[pos2]
    return x_tab, x_size, Forces


# debug function to verify results
def plot_figures(y_tab, N, Q, M, out_png, ti=[]):
    """Plot figures using the supplied output settings.

    :param y_tab: Y tab used by this operation.
    :type y_tab: array_like
    :param N: N used by this operation.
    :type N: array_like
    :param Q: Q used by this operation.
    :type Q: array_like
    :param M: M used by this operation.
    :type M: array_like
    :param out_png: Path of the image file to create.
    :type out_png: str
    :param ti: Ti used by this operation.
    :type ti: list
    """
    import pylab

    colors = ['b', 'r', 'c', 'm', 'k']
    fig, (axes_N, axes_Q, axes_M) = pylab.subplots(1, 3, sharey=True, figsize=(16, 8))
    fig.subplots_adjust(hspace=0.0, wspace=0.1, top=0.95, right=0.95)
    # pylab.axis('scaled')

    env_flag = N.ndim == 1
    if env_flag:
        size = 1
    else:
        size = N.shape[1]

    axes_N.grid(True, which="both", ls="-", lw=0.4)
    for i in range(size):
        color_index = i % len(colors)
        color = colors[color_index]
        if env_flag:
            axes_N.plot(N[:], y_tab, color=color, lw=1.75, linestyle='-')
        else:
            axes_N.plot(N[:, i], y_tab, color=color, lw=1.75, linestyle='-', label=str(ti[i]))
    #axes_Nx.plot([0.0, 0.0], ylim, color='black', lw=1.5)
    axes_N.set_xlabel('N [kN]', fontsize=18)
    axes_N.legend(loc='best')

    axes_Q.grid(True, which="both", ls="-", lw=0.4)
    for i in range(size):
        color_index = i % len(colors)
        color = colors[color_index]
        if env_flag:
            axes_Q.plot(Q[:], y_tab, color=color, lw=1.75, linestyle='-')
        else:
            axes_Q.plot(Q[:, i], y_tab, color=color, lw=1.75, linestyle='-', label=str(ti[i]))
    axes_Q.set_xlabel('Q [kN]', fontsize=18)
    axes_Q.legend(loc='best')

    axes_M.grid(True, which="both", ls="-", lw=0.4)
    for i in range(size):
        color_index = i % len(colors)
        color = colors[color_index]
        if env_flag:
            axes_M.plot(M[:], y_tab, color=color, lw=1.75, linestyle='-')
        else:
            axes_M.plot(M[:, i], y_tab, color=color, lw=1.75, linestyle='-', label=str(ti[i]))
    axes_M.set_xlabel('M [kNm]', fontsize=18)
    axes_M.legend(loc='best')

    fig.savefig(out_png)


def main():
    """Run this module's example or command-line workflow."""
    from .C_HistoryOfExecution import HistoryOfExecution
    from .C_Mesh import Mesh
    from .C_Rcf_info import RCF_info
    from .C_ShellResults import Shell_EleResults

    project = "d:/vxx_zsoil/TEMPLATES/tests/rc_pile_coarse"  # TEST-MNT-RECOVER-CONT-2"
    mesh = Mesh(project)
    rcf = RCF_info(project)
    his = HistoryOfExecution(project)

    zoom = [[], [], []]
    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]

    cont_rsl = Continuum_EleResults(mesh, his, rcf)
    mat_filter = [4]
    sel_elements = mesh.get_list_of_elements(Element.GROUP_CONTINUUM,
                                             0.0, False, mat_filter=mat_filter,
                                             ele_class_filter=[], zoom_filter=zoom_filter)

    solution_indices_aux = his.give_converged_time_solutions()
    # take only last one
    solution_indices = [1, solution_indices_aux[-1]]  # solution_indices_aux [:]

    ti = np.zeros(len(solution_indices))
    for i, solution_index in enumerate(solution_indices):
        row = his.data[solution_index]
        ti[i] = row[his.DATA_TIME]

    x_dir = np.array([0.0, -1.0, 0.0])
    pt = np.array([0.0,  0.0, 0.0])
    #pt    = np.array([0.25, 4.0, 0.5])
    y_dir = np.array([1.0, 0.0, 0.0])

    print(len(sel_elements))

    statusOK, x_tab1, x_size1, Forces1 = cont_rsl.get_MNQ_time_history_for_fict_uniform_beam(
        sel_elements, solution_indices,
        pt, x_dir, y_dir)

    #print "continuum elements size in x dir"
    #print x_size1

    shell_rsl = Shell_EleResults(mesh, his, rcf)
    mat_filter = [3]
    sel_elements = mesh.get_list_of_elements(Element.GROUP_SHELL,
                                             0.0, False, mat_filter=mat_filter,
                                             ele_class_filter=[], zoom_filter=zoom_filter)
    if len(sel_elements) > 0:
        #print len(sel_elements)

        # solution_indices = his.give_converged_time_solutions()

        statusOK, x_tab2, x_size2, Forces2 = shell_rsl.get_MNQ_time_history_for_fict_uniform_beam(
            sel_elements, solution_indices,
            pt, x_dir)

        #print "shell elements size in x dir"
        #print x_size2

        x_tab, x_size, Forces = combine_stress_resultants(x_tab1, x_size1, Forces1,
                                                          x_tab2, x_size2, Forces2)
    else:

        x_tab = x_tab1[:]
        Forces = Forces1[:]
        x_size = x_size1[:]

    Forces = Forces * 2.0

    x_tab_ex, Forces_ex = postprocess_MNQ_time_history_for_fict_uniform_beam(Forces, x_tab, x_size)
    x_tab_ex = x_tab_ex * (-1.0)

    smooth_jumps_in_stress_resultant_time_history(Forces_ex, 5)

    plot_figures(x_tab_ex, Forces_ex[:, :, 0], Forces_ex[:, :, 1], Forces_ex[:, :, 5], "Forces-test.png", ti)

    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam(Forces_ex, 0)  # NX envelope
    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam(Forces_ex, 1)  # QY envelope
    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam(Forces_ex, 5)  # MZ envelope

    #plot_figures (x_tab_ex, env_min [:, 0], env_max [:, 1], env_min [:, 5], "env-min-test.png")
    #plot_figures (x_tab_ex, env_max [:, 0], env_max [:, 1], env_max [:, 5], "env-max-test.png")

    pylab.show()


if __name__ == '__main__':
    main()
