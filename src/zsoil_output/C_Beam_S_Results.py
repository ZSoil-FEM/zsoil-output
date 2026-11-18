import os
import os
import numpy as np

from .C_EleResults     import *
from .C_Element_s_beam_L2 import *
from .C_PlasticCodes   import *

#class to get results for Spacone's beam with 5 integration points

#=====================================================
class Beam_S_EleResults (EleResults):
#=====================================================

    # =====================================================
    def __init__(self,mesh,his,rcf ):
    # =====================================================
        EleResults.__init__(self,mesh,his,rcf)

    # =====================================================
    def _get_rsl_time_history (self,ele_index_1,solution_indices,rsl_key,comp=None):
    # =====================================================
        element = self.mesh.get_element (ele_index_1)
        ngaus   = len(element.xsiGP)
        if not isinstance(element,Element_s_beam_L2):
            return None
        sel_elements = [self.mesh.get_element (ele_index_1)]
        res_all = self.get_element_results_ex (sel_elements,solution_indices)
        size = len(solution_indices)

        if comp == None:
            col_1, n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,'X')
            col_2, n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,'Y')
            col_3, n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,'Z')
            out = np.zeros(size*ngaus*3).reshape (size,ngaus,3)
        else:
            col_1, ncomps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,comp)
            col_2 = -1
            col_3 = -1
            out = np.zeros(size*ngaus).reshape (size,ngaus)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        for igaus in range (ngaus):
            res = res_all [igaus]
            offs = 2 #skip element index and gauss point index
            for i in range (size):
                if comp == None:
                    if col_1 != -1:
                        out [i][igaus][0] = res [offs+col_1]
                    if col_2 != -1:
                        out [i][igaus][1] = res [offs+col_2]
                    if col_3 != -1:
                        out [i][igaus][2] = res [offs+col_3]
                else:
                    if col_1 != -1:
                        out [i][igaus] = res [offs+col_1]
                offs = offs + gp_size_int
        return out

    # =====================================================
    def _get_rsl_for_sel_elements (self,ele_indices_1,solution_index,rsl_key,comp=None):
    # =====================================================
        element = self.mesh.get_element(ele_indices_1 [0])
        ngaus = len(element.xsiGP)

        solution_indices = [solution_index]
        size = len(ele_indices_1)
        if comp == None:
            out = np.zeros (size*ngaus*3).reshape (size,ngaus,3)
        else:
            out = np.zeros (size*ngaus).reshape(size,ngaus)

        for i,ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history(ele_index_1, solution_indices, rsl_key, comp)
            for igaus in range (ngaus):
                out [i][igaus] = ret [0][igaus]
        return out

    #public functions
    # =====================================================
    def  get_forces_vec_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        #ele_index_1 : element index (starts from 1)
        #solution_indices : indices of time instances (see HistoryOfExecution class)
        #return numpy array of size (n_time_instances,ngaus,3) containing Nx,Qy,Qz forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,'FORCE',None)

    # =====================================================
    def  get_force_time_history (self,ele_index_1,solution_indices,comp):
    # =====================================================
        # ele_index_1 : element index (starts from 1)
        # solution_indices : indices of time instances (see HistoryOfExecution class)
        # comp : one of 'X', 'Y', 'Z'
        # return numpy array of size (n_time_instances,ngaus) containing given force component
        return self._get_rsl_time_history (ele_index_1,solution_indices,'FORCE',comp)

    # =====================================================
    def  get_forces_vec_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # return numpy array of size (n_elements,ngaus,3)
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'FORCE',None)


    # =====================================================
    def  get_force_for_sel_elements (self,ele_indices_1,solution_index,comp):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # return numpy array of size (n_elements,ngaus)
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'FORCE',comp)


    #here moments - same explanation as for forces (see above)

    # =====================================================
    def  get_moments_vec_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        return self._get_rsl_time_history (ele_index_1,solution_indices,'MOMENT',None)

    # =====================================================
    def  get_moment_time_history (self,ele_index_1,solution_indices,comp):
    # =====================================================
        return self._get_rsl_time_history (ele_index_1,solution_indices,'MOMENT',comp)

    # =====================================================
    def  get_moments_vec_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'MOMENT',None)

    # =====================================================
    def  get_moment_for_sel_elements (self,ele_indices_1,solution_index,comp):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'MOMENT',comp)

    # =====================================================
    def _get_envelope_for_sel_elements(self, ele_indices_1, solution_indices, comp_index_6):
    # =====================================================
        # returns 2 matrices in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        element = self.mesh.get_element(ele_indices_1[0])
        ngaus = len(element.xsiGP)

        out_max = np.zeros(len(ele_indices_1) * ngaus * 6).reshape(len(ele_indices_1),ngaus, 6)
        out_min = np.zeros(len(ele_indices_1) * ngaus * 6).reshape(len(ele_indices_1),ngaus, 6)
        aux = np.zeros(len(ele_indices_1) * ngaus * 6).reshape(len(ele_indices_1),ngaus, 6)
        vecComps = ['X', 'Y', 'Z']

        for solution_index in solution_indices:
            F = self.get_forces_vec_for_sel_elements  (ele_indices_1, solution_index)
            M = self.get_moments_vec_for_sel_elements (ele_indices_1, solution_index)
            for i in range(F.shape[0]):
                for j in range (F.shape[1]):
                    for k in range (F.shape[2]):
                        aux[i][j][k+0] = F[i][j][k]
                        aux[i][j][k+3] = M[i][j][k]

            for i in range(aux.shape[0]):
                for j in range (aux.shape[1]):
                    if aux[i][j][comp_index_6] > out_max[i][j][comp_index_6]:
                        out_max[i][j][comp_index_6] = aux[i][j][comp_index_6]
                        for k in range(aux.shape[2]):
                            if k != comp_index_6:
                                out_max[i][j][k] = aux[i][j][k]

                    if aux[i][j][comp_index_6] < out_min[i][j][comp_index_6]:
                        out_min [i][j][comp_index_6] = aux[i][j][comp_index_6]
                        for k in range(aux.shape[2]):
                            if k != comp_index_6:
                                out_min [i][j][k] = aux[i][j][k]

        return out_min, out_max


    # =====================================================
    def  get_force_envelope_at_nodes_for_sel_elements (self,ele_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,2,6) in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        # it keeps envelope for the given force component and corresponding forces/moments
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z']
        comp_index = vecComps.index (comp)
        return self._get_envelope_for_sel_elements( ele_indices_1, solution_indices, comp_index)


    # =====================================================
    def  get_moment_envelope_at_nodes_for_sel_elements (self,ele_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,2,6) in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        # it keeps envelope for the given force component and corresponding forces/moments
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z']
        comp_index = vecComps.index(comp)+3
        return self._get_envelope_for_sel_elements(ele_indices_1, solution_indices, comp_index)
