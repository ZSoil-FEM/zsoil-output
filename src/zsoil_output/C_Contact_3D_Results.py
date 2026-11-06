"""Read and post-process ZSoil FEM output results for three-dimensional contact elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Contact_3D_Results` reader handles stresses, strains, plasticity codes, and
stress levels for quadrilateral, triangular, bonded-line, and contact-tip formulations.
Tensor results use the local interface components ``XZ``, ``YZ``, and ``ZZ``.

Public element indices are 1-based, while solution indices follow
:class:`HistoryOfExecution`. Result histories are returned as NumPy arrays; their shapes
depend on whether all tensor components or a single component are requested.
"""

import os
import numpy as np

from .C_EleResults import *
from .C_Element_cnt_Q4 import *
from .C_Element_cnt_T3 import *
from .C_Element_cnt_bond_L2 import *
#from C_PlasticCodes       import *


class Contact_3D_Results (EleResults):
    SCALAR_TYPE = 0
    VECTOR_TYPE = 1
    TENSOR_TYPE = 2

    comps_tns = ['XZ', 'YZ', 'ZZ']
    comps_vec = ['X', 'Y', 'Z']
    stress_key = 'STRESSES'
    t_stress_key = 'T_STRESSES'
    strains_key = 'STRAINS'
    pl_code_key = 'PLA_CODE'
    str_lev_key = 'STR_LEVEL'

    def __init__(self, mesh, his, rcf):
        """Initialize the 3D contact-results reader for the supplied project data.

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
        :param rsl_key: One of ``stress_key``, ``t_stress_key``, ``strains_key``, ``pl_code_key`` or
            ``str_lev_key``.
        :type rsl_key: str
        :param comp: One of ``'XZ'``, ``'YZ'`` or ``'ZZ'``; ``None`` returns every tensor component.
        :type comp: Optional[str]
        :return: NumPy array of shape ``(ntimes, ngaus, 3)`` for stresses/strains when comp = None. NumPy array of shape
            ``(ntimes, ngaus)`` for selected stress/strain component. NumPy array of shape ``(ntimes, ngaus)`` for
            any other scalar result.
        :rtype: Optional[numpy.ndarray]
        """
        element = self.mesh.get_element(ele_index_1)
        ngaus = len(element.xsiGP)
        if (not isinstance(element, Element_cnt_Q4) and
            not isinstance(element, Element_cnt_T3) and
            not isinstance(element, Element_cnt_bond_L2) and
                not isinstance(element, Element_cnt_tip)):
            return

        # T_TNS = element.get_T_TNS ()
        T_VEC = element.get_T_VEC()

        comp_index = None

        sel_elements = [self.mesh.get_element(ele_index_1)]
        res_all = self.get_element_results_ex(sel_elements, solution_indices)
        size = len(solution_indices)

        if rsl_key == Contact_3D_Results.stress_key or rsl_key == Contact_3D_Results.t_stress_key or \
                rsl_key == Contact_3D_Results.strains_key:
            if comp == None:
                col_1_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XZ')
                col_2_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YZ')
                col_3_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'ZZ')
                out = np.zeros(size*ngaus*3).reshape(size, ngaus, 3)
            else:
                col_1_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XZ')
                col_2_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YZ')
                col_3_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'ZZ')
                out = np.zeros(size*ngaus).reshape(size, ngaus)
        else:
            col, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, '')
            out = np.zeros(size * ngaus).reshape(size, ngaus)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        for igaus in range(ngaus):
            res = res_all[igaus]
            offs = 2  # skip element index and gauss point index
            for i in range(size):
                if rsl_key == Contact_3D_Results.stress_key or rsl_key == Contact_3D_Results.t_stress_key or \
                        rsl_key == Contact_3D_Results.strains_key:
                    sig = np.zeros(3)
                    sig[0] = res[offs+col_1_stress]
                    sig[1] = res[offs+col_2_stress]
                    sig[2] = res[offs+col_3_stress]
                    # transform N to the current local interface system (as in shells)
                    # if T_TNS is None:
                    if T_VEC is None:
                        sig_tr = np.zeros(3)
                        sig_tr[2] = sig[2]
                    else:
                        # sig_tr =  T_TNS.dot (sig)
                        sig_tr = T_VEC.dot(sig)

                    sig_loc = sig_tr[:]

                    if comp == None:
                        out[i][igaus][0] = sig_loc[0]
                        out[i][igaus][1] = sig_loc[1]
                        if sig_loc.shape[0] > 2:
                            out[i][igaus][2] = sig_loc[2]
                    else:
                        out[i][igaus] = sig_loc[self.comps_tns.index(comp)]

                else:
                    out[i][igaus] = res[offs+col]

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
        :param rsl_key: One of ``stress_key``, ``t_stress_key``, ``strains_key``, ``pl_code_key`` or
            ``str_lev_key``.
        :type rsl_key: str
        :param comp: One of ``'XZ'``, ``'YZ'`` or ``'ZZ'``; ``None`` returns every tensor component.
        :type comp: Optional[str]
        :return: NumPy array of shape ``(nele, ngaus, 3)`` for stresses/strains when comp = None. NumPy array of shape
            ``(nele, ngaus)`` for selected stress/strain component. NumPy array of shape ``(nele, ngaus)`` for any
            other scalar result.
        :rtype: numpy.ndarray
        """
        element = self.mesh.get_element(ele_indices_1[0])
        ngaus = len(element.xsiGP)
        solution_indices = [solution_index]
        size = len(ele_indices_1)

        if rsl_key == Contact_3D_Results.stress_key or \
           rsl_key == Contact_3D_Results.t_stress_key or \
           rsl_key == Contact_3D_Results.strains_key:
            if comp == None:
                out = np.zeros(size*ngaus*3).reshape(size, ngaus, 3)
            else:
                out = np.zeros(size*ngaus).reshape(size, ngaus)
        else:
            out = np.zeros(size * ngaus).reshape(size, ngaus)

        for i, ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history(ele_index_1, [solution_index], rsl_key, comp)
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
        return self._get_rsl_for_sel_elements(ele_indices_1, solution_index,
                                              rsl_key, None)

    def get_rsl_vec_for_sel_elements_ex(self, ele_indices_1, solution_index, rsl_key, comp):
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


def main():
    """Run this module's example or command-line workflow."""
    from .C_HistoryOfExecution import HistoryOfExecution
    from .C_Mesh import Mesh
    from .C_Rcf_info import RCF_info

    pass
    #project = "d:/vxx_zsoil/TEMPLATES/tests/TEST-MNT-RECOVER-SHELLS"
    # project = "d:/vxx_zsoil/manual/Templates/ZSsoldierWalls/tutorials/TUT-SOLDIERWALL-2"  # TEST-MNT-RECOVER-CONT-2"
    # mesh = Mesh (project)
    # rcf  = RCF_info (project)
    # his  = HistoryOfExecution (project)

    # zoom = [[], [], []]
    # zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]
    # cnt_rsl = Contact_3D_Results (mesh,his,rcf)

    # interface_mat_labels = ['#Contact_soldier_pile', '#Contact_soldier_beam', '#Contact_soldier_lagging']
    # out = mesh.find_materials_by_labels(interface_mat_labels)
    # mat_filter = []
    # for material in out:
    # mat_filter.append(material.index_in_inp)

    # sel_elements = mesh.get_list_of_elements(Element.GROUP_CONTACT, \
    # 0.0, False, mat_filter=mat_filter, \
    # ele_class_filter=[], zoom_filter=zoom_filter)
    # solution_indices_aux = his.give_converged_time_solutions()
    # solution_indices = [solution_indices_aux [-1]]

    # ti = np.zeros(len(solution_indices))
    # for i, solution_index in enumerate(solution_indices):
    # row = his.data[solution_index]
    # ti[i] = row[his.DATA_TIME]

    # # set local interface base with one in plane axis along vertical one
    # dir = np.array ([0.0, -1.0, 0.0])
    # ele_indices_1 = []
    # for element in sel_elements:
    # if isinstance (element,Element_cnt_bond_L2):
    # pass
    # element.set_transf_matrices(dir)
    # ele_indices_1.append (element.index)

    # ret = cnt_rsl.get_rsl_vec_for_sel_elements_ex (ele_indices_1, solution_indices [-1],\
    # Contact_3D_Results.strains_key, 'ZZ')
    # print(ret)
