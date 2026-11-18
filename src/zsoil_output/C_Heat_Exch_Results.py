import os
import numpy as np

from .C_EleResults         import *
from .C_Element_heat_exch_L2   import *
#debug this class

#class to get results from contact interface elements
#=====================================================
class Heat_Exch_Results (EleResults):
#=====================================================
    CONVECTIVE_FLUX_KEY = 'HEAT_FLUX'
    SCALAR_TYPE = 0
    VECTOR_TYPE = 1
    TENSOR_TYPE = 2
    rsl_types  = {CONVECTIVE_FLUX_KEY    : SCALAR_TYPE}

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
        # in list of key strings = {power_key}
        # ele_index_1 - is the element index starting from 1..
        # solution_indices - list of indices of time instances to be analyzed
        # rsl_key     - is the result key string
        # comp        - component is None here
        # it returns numpy array of size [ntimes,ngaus] for scalar result

        element = self.mesh.get_element (ele_index_1)
        ngaus   = len(element.xsiGP)
        if not isinstance(element,Element_heat_exch_L2) :
            return

        sel_elements = [self.mesh.get_element (ele_index_1)]
        res_all = self.get_element_results_ex (sel_elements,solution_indices)
        size_times = len(solution_indices)


        if self.rsl_types[rsl_key] == self.TENSOR_TYPE:
            comps_aux = self.comps_tns
        elif self.rsl_types[rsl_key] == self.VECTOR_TYPE:
            comps_aux = self.comps_vec
        elif self.rsl_types[rsl_key] == self.SCALAR_TYPE:
            comps_aux = [""]

        cols = []
        for comp_ex in comps_aux:
            col, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key,comp_ex)
            if col >= 0:
                cols.append (col)

        if comp == None:
            size = len(cols)
            if size == 1:
                out = np.zeros (size_times * ngaus).reshape (size_times, ngaus)
            else:
                out = np.zeros (size_times * ngaus * size).reshape (size_times,ngaus,size)
        else:
            col, ncomps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,comp)
            cols = [col]
            out = np.zeros (size_times*ngaus).reshape(size_times,ngaus)

        gp_size_int = self.rcf.give_one_gp_size (sel_elements[0].get_group())


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
        # comp        - is None here
        # it returns numpy array of size [nele,ngaus] for any other scalar result

        element = self.mesh.get_element(ele_indices_1 [0])
        ngaus   = len(element.xsiGP)
        solution_indices = [solution_index]
        size = len(ele_indices_1)

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
    def  get_rsl_vec_time_history_ex (self, ele_index_1, solution_indices, rsl_key, comp):
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
        return self._get_rsl_for_sel_elements_ex (ele_indices_1, solution_index, \
                                               rsl_key, None)


    # =====================================================
    def  get_rsl_vec_for_sel_elements (self, ele_indices_1, solution_index, rsl_key, comp):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # returns numpy array of size (n_elements,ngaus)  with given stress/strain component
        # returns numpy array of size (n_elements,ngaus)  with given other result
        return self._get_rsl_for_sel_elements (ele_indices_1, solution_index, \
                                               rsl_key, comp)

