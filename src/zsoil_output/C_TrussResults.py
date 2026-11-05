"""Read axial-force and state results for truss elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Truss_EleResults` reader returns force, stress, strain, plasticity-code, and
stress-level histories, selected-element values, and force envelopes.
"""

import os
import numpy as np

from .C_EleResults import *
from .C_Element_truss_L2 import *
from .C_PlasticCodes import *


class Truss_EleResults (EleResults):
    def __init__(self, mesh, his, rcf):
        """Initialize the truss element results object for the supplied project data.

        :param mesh: Mesh containing the model nodes, elements, and materials.
        :type mesh: Mesh
        :param his: Execution history that maps solution indices to times and states.
        :type his: HistoryOfExecution
        :param rcf: Result configuration describing stored nodal and element fields.
        :type rcf: RCF_info
        """
        EleResults.__init__(self, mesh, his, rcf)

    def _get_rsl_time_history(self, ele_index_1, solution_indices, rsl_key, comp):
        """Return the results history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to return.
        :type comp: str
        :return: The results history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        element = self.mesh.get_element(ele_index_1)
        if not isinstance(element, Element_truss_L2):
            return None
        sel_elements = [self.mesh.get_element(ele_index_1)]
        res_all = self.get_element_results_ex(sel_elements, solution_indices)
        res = res_all[0]
        col, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(),
                                                     rsl_key, comp)
        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        size = len(solution_indices)
        out = np.zeros(size)
        offs = 2
        for i in range(size):
            if col != -1:
                out[i] = res[offs+col]
            offs = offs + gp_size_int
        return out

    def _get_rsl_for_sel_elements(self, ele_indices_1, solution_index, rsl_key, comp):
        """Return results for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param rsl_key: Name of the result item to read.
        :type rsl_key: str
        :param comp: Component label to return.
        :type comp: str
        :return: The results for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        solution_indices = [solution_index]
        size = len(ele_indices_1)
        out = np.zeros(size)
        for i, ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history(ele_index_1, solution_indices, rsl_key, comp)
            out[i] = ret[0]

        return out

    # public functions
    def get_force_time_history(self, ele_index_1, solution_indices):
        """Return the force history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The force history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'FORCE', 'X')

    def get_force_for_sel_elements(self, ele_indices_1, solution_index):
        """Return force for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The force for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', 'X')

    def get_force_envelope(self, ele_indices_1, solution_indices):
        """Return the minimum and maximum force values across the requested solution steps.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The minimum and maximum force values across the requested solution steps.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """
        env_min = np.zeros(len(ele_indices_1))
        env_max = np.zeros(len(ele_indices_1))

        for solution_index in solution_indices:
            out = self.get_force_for_sel_elements(ele_indices_1, solution_index)
            for i in range(out.shape[0]):
                env_min[i] = min(env_min[i], out[i])
                env_max[i] = max(env_max[i], out[i])
        return env_min, env_max

    def get_stress_time_history(self, ele_index_1, solution_indices):
        """Return the stress history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The stress history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'STRESSES', 'XX')

    def get_stress_for_sel_elements(self, ele_indices_1, solution_index):
        """Return stress for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The stress for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'STRESSES', 'XX')

    def get_strain_time_history(self, ele_index_1, solution_indices):
        """Return the strain history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The strain history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'STRAINS', 'XX')

    def get_strain_for_sel_elements(self, ele_indices_1, solution_index):
        """Return strain for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The strain for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'STRAINS', 'XX')

    def get_plast_code_time_history(self, ele_index_1, solution_indices):
        """Return the plasticity code history across the requested solution steps.

        .. seealso::
           :class:`PlasticCodes` for the meaning of individual bit flags.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The plasticity code history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'PLA_CODE', '')

    def get_plast_code_for_sel_elements(self, ele_indices_1, solution_index):
        """Return plasticity code for the selected elements at one solution step.

        .. seealso::
           :class:`PlasticCodes` for the meaning of individual bit flags.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The plasticity code for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'PLA_CODE', '')

    def get_stress_level_time_history(self, ele_index_1, solution_indices):
        """Return the stress level history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The stress level history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'STR_LEVEL', '')

    def get_stress_level_for_sel_elements(self, ele_indices_1, solution_index):
        """Return stress level for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The stress level for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'STR_LEVEL', '')
