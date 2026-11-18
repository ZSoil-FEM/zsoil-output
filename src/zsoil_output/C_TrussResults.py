import os
import numpy as np

from .C_EleResults     import *
from .C_Element_truss_L2 import *
from .C_PlasticCodes   import *



#=====================================================
class Truss_EleResults (EleResults):
#=====================================================

    # =====================================================
    def __init__(self,mesh,his,rcf ):
    # =====================================================
        EleResults.__init__(self,mesh,his,rcf)

    # =====================================================
    def _get_rsl_time_history (self,ele_index_1,solution_indices,rsl_key,comp):
    # =====================================================
        element = self.mesh.get_element (ele_index_1)
        if not isinstance(element,Element_truss_L2):
            return None
        sel_elements = [self.mesh.get_element (ele_index_1)]
        res_all = self.get_element_results_ex (sel_elements,solution_indices)
        res = res_all [0]
        col,n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(), \
                                                     rsl_key,comp)
        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        size = len(solution_indices)
        out  = np.zeros (size)
        offs = 2
        for i in range (size):
            if col != -1:
                out [i] = res [offs+col]
            offs = offs + gp_size_int
        return out

    # =====================================================
    def _get_rsl_for_sel_elements (self,ele_indices_1,solution_index,rsl_key,comp):
    # =====================================================
        solution_indices = [solution_index]
        size = len(ele_indices_1)
        out = np.zeros (size)
        for i,ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history (ele_index_1,solution_indices,rsl_key,comp)
            out [i] = ret [0]

        return out

    #public functions
    # =====================================================
    def  get_force_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        return self._get_rsl_time_history (ele_index_1,solution_indices,'FORCE','X')

    # =====================================================
    def  get_force_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'FORCE','X')

    # =====================================================
    def  get_force_envelope (self,ele_indices_1,solution_indices):
    # =====================================================
        #for given list of elements and time instances create min/max force envelopes
        env_min = np.zeros (len(ele_indices_1))
        env_max = np.zeros (len(ele_indices_1))

        for solution_index in solution_indices:
            out = self.get_force_for_sel_elements (ele_indices_1,solution_index)
            for i in range (out.shape[0]):
                env_min [i] = min (env_min [i],out [i])
                env_max [i] = max (env_max [i],out [i])
        return env_min,env_max

    # =====================================================
    def  get_stress_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'STRESSES', 'XX')

    # =====================================================
    def  get_stress_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'STRESSES','XX')

    # =====================================================
    def  get_strain_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        return self._get_rsl_time_history(ele_index_1, solution_indices, 'STRAINS', 'XX')

    # =====================================================
    def  get_strain_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'STRAINS','XX')

    # =====================================================
    def  get_plast_code_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        #see PLasticCodes class in C_PlasticCodes.py for the meaning
        return self._get_rsl_time_history (ele_index_1,solution_indices,'PLA_CODE','')

    # =====================================================
    def  get_plast_code_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        # see PLasticCodes class in C_PlasticCodes.py for the meaning
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'PLA_CODE','')

    # =====================================================
    def  get_stress_level_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        return self._get_rsl_time_history (ele_index_1,solution_indices,'STR_LEVEL','')

    # =====================================================
    def  get_stress_level_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,'STR_LEVEL','')



