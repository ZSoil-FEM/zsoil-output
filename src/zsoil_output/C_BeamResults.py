"""Read force and moment results for standard beam elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Beam_EleResults` reader provides histories, selected-element values, nodal
values, and minimum/maximum envelopes for generalized beam forces and moments.
"""

import os
import numpy as np

from .C_EleResults import *
from .C_Element_beam_L2 import *
from .C_PlasticCodes import *

class Beam_EleResults (EleResults):
    def __init__(self, mesh, his, rcf):
        """Initialize the beam element results object for the supplied project data.

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
        if not isinstance(element, Element_beam_L2):
            return None
        sel_elements = [self.mesh.get_element(ele_index_1)]
        res_all = self.get_element_results_ex(sel_elements, solution_indices)
        res = res_all[0]
        size = len(solution_indices)

        if comp == None:
            col_1, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'X')
            col_2, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'Y')
            col_3, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'Z')
            out = np.zeros(size*3).reshape(size, 3)
        else:
            col_1, ncomps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, comp)
            col_2 = -1
            col_3 = -1
            out = np.zeros(size)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        offs = 2  # skip element index and gauss point index
        for i in range(size):
            if comp == None:
                if col_1 != -1:
                    out[i][0] = res[offs+col_1]
                if col_2 != -1:
                    out[i][1] = res[offs+col_2]
                if col_3 != -1:
                    out[i][2] = res[offs+col_3]
            else:
                if col_1 != -1:
                    out[i] = res[offs+col_1]
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
        solution_indices = [solution_index]
        size = len(ele_indices_1)
        if comp == None:
            out = np.zeros(size*3).reshape(size, 3)
        else:
            out = np.zeros(size)

        for i, ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history(ele_index_1, solution_indices, rsl_key, comp)
            out[i] = ret[0]

        return out

    # public functions
    def get_forces_vec_time_history(self, ele_index_1, solution_indices):
        """Return the forces vector history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: NumPy array of shape ``(n_time_instances, 3)`` containing Nx,Qy,Qz forces.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'FORCE', None)

    def get_force_time_history(self, ele_index_1, solution_indices, comp):
        """Return the force history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param comp: One of 'X', 'Y', 'Z'.
        :type comp: str
        :return: NumPy array of shape ``(n_time_instances,)`` containing given force component.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'FORCE', comp)

    def get_forces_vec_for_sel_elements(self, ele_indices_1, solution_index):
        """Return forces vector for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: NumPy array of shape ``(n_elements, 3)``.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', None)

    def get_force_for_sel_elements(self, ele_indices_1, solution_index, comp):
        """Return force for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param comp: Component label to return.
        :type comp: str
        :return: NumPy array of shape ``(n_elements,)``.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', comp)

    # here moments - same explanation as for forces (see above)

    def get_moments_vec_time_history(self, ele_index_1, solution_indices):
        """Return the moments vector history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The moments vector history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'MOMENT', None)

    def get_moment_time_history(self, ele_index_1, solution_indices, comp):
        """Return the moment history across the requested solution steps.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param comp: Component label to return.
        :type comp: str
        :return: The moment history across the requested solution steps.
        :rtype: Optional[numpy.ndarray]
        """
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'MOMENT', comp)

    def get_moments_vec_for_sel_elements(self, ele_indices_1, solution_index):
        """Return moments vector for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The moments vector for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'MOMENT', None)

    def get_moment_for_sel_elements(self, ele_indices_1, solution_index, comp):
        """Return moment for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param comp: Component label to return.
        :type comp: str
        :return: The moment for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'MOMENT', comp)

    def get_force_at_nodes_for_sel_elements(self, ele_indices_1, solution_index, comp):
        """Return force at nodes for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param comp: Component label to return.
        :type comp: str
        :return: The force at nodes for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        F = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', comp)
        out = np.zeros(F.shape[0] * 2)
        offs = 0
        for i in range(F.shape[0]):
            out[offs+0] = F[i]
            out[offs+1] = F[i]
            offs = offs + 2
        return out

    def get_forces_vec_at_nodes_for_sel_elements(self, ele_indices_1, solution_index):
        """Return forces vector at nodes for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The forces vector at nodes for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        F = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', None)
        out = np.zeros(F.shape[0] * 2 * 3).reshape(F.shape[0] * 2, 3)
        offs = 0
        for i in range(F.shape[0]):
            out[offs+0] = F[i]
            out[offs+1] = F[i]
            offs = offs + 2
        return out

    def get_moment_at_nodes_for_sel_elements(self, ele_indices_1, solution_index, comp):
        """Return moment at nodes for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :param comp: Component label to return.
        :type comp: str
        :return: The moment at nodes for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        F = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', None)
        M = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'MOMENT', comp)
        out = np.zeros(M.shape[0] * 2)
        offs = 0
        for i in range(M.shape[0]):
            if comp == 'X':
                out[offs+0] = M[i]
                out[offs+1] = M[i]
            elif comp == 'Y':
                element = self.mesh.get_element(ele_indices_1[i])
                L = element.get_ele_dim()
                Qz = F[i][2]
                out[offs+0] = M[i] - Qz * L / 2.0
                out[offs+1] = M[i] + Qz * L / 2.0
            elif comp == 'Z':
                element = self.mesh.get_element(ele_indices_1[i])
                L = element.get_ele_dim()
                Qy = F[i][1]
                out[offs+0] = M[i] - Qy * L / 2.0
                out[offs+1] = M[i] + Qy * L / 2.0
            offs = offs + 2
        return out

    def get_moments_vec_at_nodes_for_sel_elements(self, ele_indices_1, solution_index):
        """Return moments vector at nodes for the selected elements at one solution step.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_index: Solution index defined by :class:`HistoryOfExecution`.
        :type solution_index: int
        :return: The moments vector at nodes for the selected elements at one solution step.
        :rtype: numpy.ndarray
        """
        F = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'FORCE', None)
        M = self._get_rsl_for_sel_elements(ele_indices_1, solution_index, 'MOMENT', None)
        out = np.zeros(M.shape[0] * 2 * 3).reshape(M.shape[0] * 2, 3)
        offs = 0
        for i in range(M.shape[0]):
            element = self.mesh.get_element(ele_indices_1[i])
            L = element.get_ele_dim()
            out[offs+0][0] = M[i][0]
            out[offs+1][0] = M[i][0]
            Qz = F[i][2]
            out[offs+0][1] = M[i][1] - Qz * L / 2.0
            out[offs+1][1] = M[i][1] + Qz * L / 2.0
            Qy = F[i][1]
            out[offs+0][2] = M[i][2] - Qy * L / 2.0
            out[offs+1][2] = M[i][2] + Qy * L / 2.0

            offs = offs + 2
        return out

    def _get_envelope_at_nodes_for_sel_elements(self, ele_indices_1, solution_indices, comp_index_6):
        """Return the minimum and maximum result values across the requested solution steps.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param comp_index_6: Component index 6 used by this operation.
        :type comp_index_6: int
        :return: A ``(minimum, maximum)`` tuple of NumPy arrays in which each row consists of
            ``[NX, QY, QZ, MX, MY, MZ]`` generalized force
            components.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """
        out_max = np.zeros(len(ele_indices_1) * 2 * 6).reshape(len(ele_indices_1)*2, 6)
        out_min = np.zeros(len(ele_indices_1) * 2 * 6).reshape(len(ele_indices_1)*2, 6)
        aux = np.zeros(len(ele_indices_1) * 2 * 6).reshape(len(ele_indices_1)*2, 6)
        vecComps = ['X', 'Y', 'Z']

        for solution_index in solution_indices:
            F = self.get_forces_vec_at_nodes_for_sel_elements(ele_indices_1, solution_index)
            M = self.get_moments_vec_at_nodes_for_sel_elements(ele_indices_1, solution_index)
            offs = 0
            for i in range(F.shape[0]):
                for k in range(3):
                    aux[i, k+0] = F[i, k]
                    aux[i, k+3] = M[i, k]

            offs = 0
            for i in range(aux.shape[0]):
                if aux[i, comp_index_6] > out_max[i, comp_index_6]:
                    out_max[i, comp_index_6] = aux[i, comp_index_6]
                    for k in range(6):
                        if k != comp_index_6:
                            out_max[i, k] = aux[i, k]

                if aux[i, comp_index_6] < out_min[i, comp_index_6]:
                    out_min[i, comp_index_6] = aux[i, comp_index_6]
                    for k in range(6):
                        if k != comp_index_6:
                            out_min[i, k] = aux[i, k]

        return out_min, out_max

    def get_force_envelope_at_nodes_for_sel_elements(self, ele_indices_1, solution_indices, comp):
        """Return the minimum and maximum force values across the requested solution steps.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param comp: Component label to return.
        :type comp: str
        :return: A ``(minimum, maximum)`` tuple of NumPy arrays of shape ``(N_ele, 2, 6)`` in which each row consists of
            ``[NX, QY, QZ, MX, MY, MZ]``
            generalized force components.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """
        vecComps = ['X', 'Y', 'Z']
        comp_index = vecComps.index(comp)
        return self._get_envelope_at_nodes_for_sel_elements(ele_indices_1, solution_indices, comp_index)

    def get_moment_envelope_at_nodes_for_sel_elements(self, ele_indices_1, solution_indices, comp):
        """Return the minimum and maximum moment values across the requested solution steps.

        :param ele_indices_1: 1-based element indices.
        :type ele_indices_1: Sequence[int]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param comp: Component label to return.
        :type comp: str
        :return: A ``(minimum, maximum)`` tuple of NumPy arrays of shape ``(N_ele, 2, 6)`` in which each row consists of
            ``[NX, QY, QZ, MX, MY, MZ]``
            generalized force components.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """
        vecComps = ['X', 'Y', 'Z']
        comp_index = vecComps.index(comp)+3
        return self._get_envelope_at_nodes_for_sel_elements(ele_indices_1, solution_indices, comp_index)
