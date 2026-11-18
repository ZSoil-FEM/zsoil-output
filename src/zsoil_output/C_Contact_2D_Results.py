import os
import numpy as np

from .C_EleResults         import *
from .C_Element_cnt_L2   import *
from .C_Element_cnt_bond_L2 import *
#from C_PlasticCodes       import *

#debug this class

#class to get results from contact interface elements
#=====================================================
class Contact_2D_Results (EleResults):
#=====================================================
    SCALAR_TYPE = 0
    VECTOR_TYPE = 1
    TENSOR_TYPE = 2

    comps_tns    = ['XY', 'YY']
    comps_vec    = ['X', 'Y']
    stress_key   = 'STRESSES'
    t_stress_key = 'T_STRESSES'
    strains_key  = 'STRAINS'
    pl_code_key  = 'PLA_CODE'
    str_lev_key  = 'STR_LEVEL'


    # =====================================================
    def __init__(self,mesh,his,rcf ):
    # =====================================================
        EleResults.__init__(self,mesh,his,rcf)

    # =====================================================
    def _get_rsl_time_history (self,ele_index_1,solution_indices,rsl_key,comp=None):
    # =====================================================
        # this function should not be used by the user
        # it is treated as a private one
        # this function returns time history of one of the stored results identified
        # in list of key strings =
        # {stress_key, t_stress_key, strains_key, pl_code_key, str_lev_key}
        # ele_index_1 - is the element index starting from 1..
        # solution_indices - list of indices of time instances to be analyzed
        # rsl_key     - is the result key string
        # comp        - component is None or one from set comps_tns = ['XZ', 'YZ', 'ZZ']
        # it returns numpy array of size [ntimes,ngaus,2] for stresses/strains when comp = None
        # it returns numpy array of size [ntimes,ngaus] for selected stress/strain component
        # it returns numpy array of size [ntimes,ngaus] for any other scalar result

        element = self.mesh.get_element (ele_index_1)
        ngaus   = len(element.xsiGP)
        if not (isinstance(element,Element_cnt_L2) or isinstance(element,Element_cnt_bond_L2)) :
            return

        comp_index = None

        sel_elements = [self.mesh.get_element (ele_index_1)]
        res_all = self.get_element_results_ex (sel_elements,solution_indices)
        size = len(solution_indices)

        if rsl_key == Contact_2D_Results.stress_key or rsl_key == Contact_2D_Results.t_stress_key or \
            rsl_key == Contact_2D_Results.strains_key:
            if comp == None:
                col_1_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XY')
                col_2_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YY')
                out = np.zeros(size*ngaus*2).reshape (size,ngaus,2)
            else:
                col_1_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XY')
                col_2_stress, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YY')
                out = np.zeros(size*ngaus).reshape (size,ngaus)
        else:
            col, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, "")
            out = np.zeros(size * ngaus).reshape(size, ngaus)

        gp_size_int = self.rcf.give_one_gp_size (sel_elements[0].get_group())

        for igaus in range (ngaus):
            res = res_all [igaus]
            offs = 2 #skip element index and gauss point index
            for i in range (size):
                if rsl_key == Contact_2D_Results.stress_key or rsl_key == Contact_2D_Results.t_stress_key or \
                        rsl_key == Contact_2D_Results.strains_key:
                    sig = np.zeros(2)
                    sig [0] = res [offs+col_1_stress]
                    sig [1] = res [offs+col_2_stress]
                    if comp == None:
                        out [i][igaus][0] = sig [0]
                        out [i][igaus][1] = sig [1]
                    else:
                        out [i][igaus] = sig [Contact_2D_Results.comps_tns.index(comp)]
                else:
                    out [i][igaus] = res [offs+col]

                offs = offs + gp_size_int
        return out

    # =====================================================
    def _get_rsl_for_sel_elements (self, ele_indices_1, solution_index, rsl_key,\
                                   comp=None):
    # =====================================================
        # this function should not be used by the user
        # it is treated as a private one
        # this function returns required result identified by the rsl_key
        # from list of key strings =
        # {stress_key, t_stress_key, strains_key, pl_code_key, str_lev_key}
        # ele_indices_1 - list of element indices starting from 1
        # solution_index - index of time instance to be analyzed
        # rsl_key     - is the result key string
        # comp        - component is None or one from set comps_tns = ['XY', 'YY']
        # it returns numpy array of size [nele,ngaus,2] for stresses/strains when comp = None
        # it returns numpy array of size [nele,ngaus] for selected stress/strain component
        # it returns numpy array of size [nele,ngaus] for any other scalar result

        element = self.mesh.get_element(ele_indices_1 [0])
        ngaus   = len(element.xsiGP)
        solution_indices = [solution_index]
        size = len(ele_indices_1)

        if rsl_key == Contact_2D_Results.stress_key or \
           rsl_key == Contact_2D_Results.t_stress_key or \
           rsl_key == Contact_2D_Results.strains_key:
            if comp == None:
                out = np.zeros (size*ngaus*2).reshape (size,ngaus,2)
            else:
                out = np.zeros (size*ngaus).reshape(size,ngaus)
        else:
            out = np.zeros(size * ngaus).reshape(size, ngaus)

        for i,ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history (ele_index_1,solution_indices,rsl_key,comp)
            for igaus in range (ret.shape[1]):
                if ret.ndim == 2:
                    out[i][igaus] = ret[0][igaus]
                else:
                    for k in range (ret.shape [2]):#loop over components
                        out [i][igaus][k] = ret [0][igaus][k]
        return out

    # --------------------------------------
    # public functions
    # these function can be used by the user
    # --------------------------------------
    # =====================================================
    def  get_rsl_vec_time_history (self, ele_index_1, solution_indices, rsl_key):
    # =====================================================
        #ele_index_1 : element index (starts from 1)
        #solution_indices : indices of time instances (see HistoryOfExecution class)
        #return numpy array of size (n_time_instances,ngaus,3) for stresses/strains
        # return numpy array of size (n_time_instances,ngaus) for other
        return self._get_rsl_time_history (ele_index_1, solution_indices, \
                                           rsl_key, None)


    # =====================================================
    def  get_rsl_vec_time_history (self, ele_index_1, solution_indices, rsl_key, comp):
    # =====================================================
        # ele_index_1 : element index (starts from 1)
        # solution_indices : indices of time instances (see HistoryOfExecution class)
        # returns numpy array of size (n_time_instances,ngaus,3) for stresses/strains
        # returns numpy array of size (n_time_instances,ngaus) for other
        return self._get_rsl_time_history (ele_index_1, solution_indices, \
                                           rsl_key, comp)

    # =====================================================
    def  get_rsl_vec_for_sel_elements (self, ele_indices_1, solution_index, rsl_key):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # returns numpy array of size (n_elements,ngaus,ncomps) for strains/stresses
        # returns numpy array of size (n_elements,ngaus) for other
        return self._get_rsl_for_sel_elements (ele_indices_1, solution_index, \
                                               rsl_key, None)


    # =====================================================
    def  get_rsl_vec_for_sel_elements_ex (self, ele_indices_1, solution_index, rsl_key, comp):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # returns numpy array of size (n_elements,ngaus)  with given stress/strain component
        # returns numpy array of size (n_elements,ngaus)  with given other result
        return self._get_rsl_for_sel_elements (ele_indices_1, solution_index, \
                                               rsl_key, comp)

#=====================================================
def main():
#=====================================================
    from .C_HistoryOfExecution import HistoryOfExecution
    from .C_Mesh import Mesh
    from .C_Rcf_info import RCF_info
    project = "d:/vxx_zsoil/manual/Templates/ZSsoldierWalls/tutorials/TUT-SOLDIERWALL-2"  # TEST-MNT-RECOVER-CONT-2"
    mesh = Mesh (project)
    rcf  = RCF_info (project)
    his  = HistoryOfExecution (project)

    zoom = [[], [], []]
    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]
    cnt_rsl = Contact_2D_Results (mesh,his,rcf)
    mat_filter = []
    sel_elements = mesh.get_list_of_elements(Element.GROUP_CONTACT, \
                                              0.0, False, mat_filter=mat_filter, \
                                              ele_class_filter=["C_L2"], zoom_filter=zoom_filter)
    solution_indices_aux = his.give_converged_time_solutions()
    solution_indices = [solution_indices_aux [-1]]

    ti = np.zeros(len(solution_indices))
    for i, solution_index in enumerate(solution_indices):
        row = his.data[solution_index]
        ti[i] = row[his.DATA_TIME]

    ele_indices_1 = []
    for element in sel_elements:
        ele_indices_1.append (element.index)

    ret = cnt_rsl.get_rsl_vec_for_sel_elements (ele_indices_1, solution_indices [-1], \
                                                Contact_2D_Results.strains_key, 'YY')
    print(ret)


if __name__ == '__main__':
    main()
