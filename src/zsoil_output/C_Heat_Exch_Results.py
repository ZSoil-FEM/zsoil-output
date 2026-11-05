"""Read heat-flux results for heat-exchange elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Heat_Exch_Results` reader returns convective heat-flux histories and
selected-element values from heat-exchange result files.
"""

import os
import numpy as np

from .C_EleResults import *
from .C_Element_heat_exch_L2 import *


class Heat_Exch_Results (EleResults):
    CONVECTIVE_FLUX_KEY = 'HEAT_FLUX'
    SCALAR_TYPE = 0
    VECTOR_TYPE = 1
    TENSOR_TYPE = 2
    rsl_types = {CONVECTIVE_FLUX_KEY: SCALAR_TYPE}

    def __init__(self, mesh, his, rcf):
        """Initialize the heat exch results object for the supplied project data.

        :param mesh: Mesh containing the model nodes, elements, and materials.
        :type mesh: Mesh
        :param his: Execution history that maps solution indices to times and states.
        :type his: HistoryOfExecution
        :param rcf: Result configuration describing stored nodal and element fields.
        :type rcf: RCF_info
        """
        EleResults.__init__(self, mesh, his, rcf)

    def _get_rsl_time_history(self, ele_index_1, solution_indices, rsl_key, comp=None):
        """Return the results history across the requested solution steps.

        .. note::
           This is an internal helper and is not part of the public API.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Result key; currently only ``CONVECTIVE_FLUX_KEY`` is supported.
        :type rsl_key: str
        :param comp: Component selector; use ``None`` for the scalar heat-flux result.
        :type comp: Optional[str]
        :return: NumPy array of shape ``(ntimes, ngaus)`` for scalar result.
        :rtype: Optional[numpy.ndarray]
        """
        element = self.mesh.get_element(ele_index_1)
        ngaus = len(element.xsiGP)
        if not isinstance(element, Element_heat_exch_L2):
            return

        sel_elements = [self.mesh.get_element(ele_index_1)]
        res_all = self.get_element_results_ex(sel_elements, solution_indices)
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
            cols = [col]
            out = np.zeros(size_times*ngaus).reshape(size_times, ngaus)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

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

        .. note::
           This is an internal helper and is not part of the public API.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Result key; currently only ``CONVECTIVE_FLUX_KEY`` is supported.
        :type rsl_key: str
        :param comp: Component selector; use ``None`` for the scalar heat-flux result.
        :type comp: Optional[str]
        :return: NumPy array of shape ``(nele, ngaus)`` for any other scalar result.
        :rtype: numpy.ndarray
        """
        element = self.mesh.get_element(ele_indices_1[0])
        ngaus = len(element.xsiGP)
        solution_indices = [solution_index]
        size = len(ele_indices_1)

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

    def get_rsl_vec_time_history(self, ele_index_1, solution_indices, rsl_key):
        """Return the results vector history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :return: NumPy array of shape ``(n_time_instances, ngaus, 3)`` for stresses/strains. NumPy array of shape
            ``(n_time_instances, ngaus)`` for other.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices,
                                          rsl_key, None)

    def get_rsl_vec_time_history_ex(self, ele_index_1, solution_indices, rsl_key, comp):
        """Return the results vector history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to return.
        :type comp: str
        :return: NumPy array of shape ``(n_time_instances, ngaus, 3)`` for stresses/strains. NumPy array of shape
            ``(n_time_instances, ngaus)`` for other.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices,
                                          rsl_key, comp)

    def get_rsl_vec_for_sel_elements(self, ele_indices_1, solution_index, rsl_key):
        """Return results vector for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :return: NumPy array of shape ``(n_elements, ngaus, ncomps)`` for strains/stresses. NumPy array of shape
            ``(n_elements, ngaus)`` for other.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements_ex(ele_indices_1, solution_index,
                                                 rsl_key, None)

    def get_rsl_vec_for_sel_elements(self, ele_indices_1, solution_index, rsl_key, comp):
        """Return results vector for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to return.
        :type comp: str
        :return: NumPy array of shape ``(n_elements, ngaus)`` with given stress/strain component. NumPy array of shape
            ``(n_elements, ngaus)`` with given other result.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index,
                                              rsl_key, comp)
