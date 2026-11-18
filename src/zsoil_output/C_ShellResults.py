import os
import numpy as np

from .C_EleResults         import *
from .C_Element_shell_Q4   import *
from .C_Element_shell_T3   import *
from .C_Element_shell_SHQ4 import *
#from C_PlasticCodes       import *

#debug this class


#class to get results for shell elements

#=====================================================
class Shell_EleResults (EleResults):
#=====================================================

    comps_tns  = ['XX', 'YY', 'XY']
    comps_vec  = ['X', 'Y']
    comps_all  = ['NXX','NYY','NXY','MXX','MYY','MXY','QX','QY']
    comps_all_AW = ['NXX','NYY','MXX','MYY','QX', 'QY']
    sNforce_key= 'SMFORCE'
    sQforce_key= 'SQFORCE'
    sMoment_key= 'SMOMENT'
    thick_key  = 'THICK'

    # =====================================================
    def __init__(self,mesh,his,rcf ):
    # =====================================================
        EleResults.__init__(self,mesh,his,rcf)

    # =====================================================
    def Armer_Wood_MN (self,Mxx,Myy,Mxy,fiber_string):
    # =====================================================
        #this function computes Armer Wood moments (top(+)/bottom(-) and membrane forces (+/-)
        #fiber string is '+' (TOP) or '-' (BOT)
        Mxx_AW = 0.0
        Myy_AW = 0.0
        EPS    = 1.0e-10

        if fiber_string == '+' or fiber_string == 'T':
            #----------------------------------------
            Mxx_AW = Mxx + abs (Mxy)
            if Myy + abs(Mxy) < -EPS:
                Mxx_AW = Mxx + Mxy * Mxy / abs(Myy)
            if Mxx_AW < 0.0:
                Mxx_AW = 0.0

            Myy_AW = Myy + abs(Mxy)
            if Mxx + abs(Mxy) < -EPS:
                Myy_AW = Myy + Mxy * Mxy / abs(Mxx)
            if Myy_AW < 0.0:
                Myy_AW = 0.0

        elif fiber_string == '-' or fiber_string == 'B':
            #----------------------------------------
            Mxx_AW = Mxx - abs(Mxy)
            if Myy - abs(Mxy) > EPS:
                Mxx_AW = Mxx - Mxy * Mxy / abs(Myy)
            if Mxx_AW > 0.0:
                Mxx_AW = 0.0

            Myy_AW = Myy - abs(Mxy)
            if Mxx - abs(Mxy) > EPS:
                Myy_AW = Myy - Mxy * Mxy / abs(Mxx)
            if Myy_AW > 0.0:
                Myy_AW = 0.0

        out = np.zeros (2)
        out [0] = Mxx_AW
        out [1] = Myy_AW
        return out

    # =====================================================
    def _get_rsl_time_history (self,ele_index_1,solution_indices,rsl_key,comp=None,\
                               Armer_Wood_flag=False, fiber_flag = ''):
    # =====================================================
        element = self.mesh.get_element (ele_index_1)
        ngaus   = len(element.xsiGP)
        if not isinstance(element,Element_shell_Q4) and not isinstance(element,Element_shell_SHQ4) \
                and not isinstance(element,Element_shell_T3):
            return None

        T_TNS = element.get_T_TNS ()
        T_VEC = element.get_T_VEC ()

        comps_MN = Shell_EleResults.comps_tns
        comps_Q  = Shell_EleResults.comps_vec
        comp_index = None

        sel_elements = [self.mesh.get_element (ele_index_1)]
        res_all = self.get_element_results_ex (sel_elements,solution_indices)
        size = len(solution_indices)

        if rsl_key == Shell_EleResults.sNforce_key or rsl_key == Shell_EleResults.sMoment_key:
            col_1N, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), Shell_EleResults.sNforce_key, 'XX')
            col_2N, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), Shell_EleResults.sNforce_key, 'YY')
            col_3N, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), Shell_EleResults.sNforce_key, 'XY')
            col_1M, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), Shell_EleResults.sMoment_key, 'XX')
            col_2M, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), Shell_EleResults.sMoment_key, 'YY')
            col_3M, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), Shell_EleResults.sMoment_key, 'XY')
            if comp == None:
                col_1, n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,'XX')
                col_2, n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,'YY')
                col_3, n_comps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,'XY')
                if Armer_Wood_flag:
                    out = np.zeros(size*ngaus*2).reshape (size,ngaus,2)
                else:
                    out = np.zeros(size*ngaus*3).reshape (size,ngaus,3)
            else:
                col_1, ncomps = self.rcf.give_col_for_ele_rsl (sel_elements[0].get_group(),rsl_key,comp)
                col_2 = -1
                col_3 = -1
                out = np.zeros(size*ngaus).reshape (size,ngaus)
        else:
            col_1Q, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'X')
            col_2Q, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'Y')
            col_1Q, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'X')
            col_2Q, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'Y')

            if comp == None:
                col_1, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'X')
                col_2, n_comps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'Y')
                col_3 = -1
                out = np.zeros(size * ngaus * 2).reshape(size, ngaus, 2)
            else:
                col_1, ncomps = self.rcf.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, comp)
                col_2 = -1
                col_3 = -1
                out = np.zeros(size * ngaus).reshape(size, ngaus)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        for igaus in range (ngaus):
            res = res_all [igaus]
            offs = 2 #skip element index and gauss point index
            for i in range (size):

                if rsl_key == Shell_EleResults.sNforce_key:
                    N = np.zeros(3)
                    N [0] = res [offs+col_1N]
                    N [1] = res [offs+col_2N]
                    N [2] = res [offs+col_3N]
                    #transform N to the current local shell system
                    #(no transform set via set_transf_matrices() -> fall back
                    #to the element's own untransformed local axes rather than
                    #discarding the real, correctly-read result as zero):
                    if T_TNS is None:
                        N_tr = N
                    else:
                        N_tr =  T_TNS.dot (N)
                    if Armer_Wood_flag:
                        N_loc = self.Armer_Wood_MN (N_tr [0], N_tr [1], N_tr [2], \
                                                    fiber_flag)
                    else:
                        N_loc = N_tr

                    if comp == None:
                        out [i][igaus][0] = N_loc [0]
                        out [i][igaus][1] = N_loc [1]
                        if N_loc.shape [0] > 2:
                            out [i][igaus][2] = N_loc [2]
                    else:
                        out [i][igaus] = N_loc [comps_MN.index(comp)]

                elif rsl_key == Shell_EleResults.sMoment_key:
                    M = np.zeros(3)
                    M [0] = res [offs + col_1M]
                    M [1] = res [offs + col_2M]
                    M [2] = res [offs + col_3M]
                    #transform M to the current local shell system
                    #(same fallback as for N above - see comment there):
                    if T_TNS is None:
                        M_tr = M
                    else:
                        M_tr =  T_TNS.dot (M)
                    if Armer_Wood_flag:
                        M_loc = self.Armer_Wood_MN (M_tr[0], M_tr[1], M_tr[2],\
                                                    fiber_flag)
                    else:
                        M_loc = M_tr
                    if comp == None:
                        out [i][igaus][0] = M_loc [0]
                        out [i][igaus][1] = M_loc [1]
                        if M_loc.shape[0] > 2:
                            out [i][igaus][2] = M_loc [2]
                    else:
                        out [i][igaus] = M_loc [comps_MN.index(comp)]

                else:

                    Q = np.zeros(3)
                    Q [0] = res [offs + col_1Q]
                    Q [1] = res [offs + col_2Q]
                    Q [2] = 0.0
                    #transform Q to the current local shell system
                    #(same fallback as for N above - see comment there):
                    if T_VEC is None:
                        Q_loc = Q
                    else:
                        Q_loc = T_VEC.dot (Q)
                    if comp == None:
                        out [i][igaus][0] = Q_loc [0]
                        out [i][igaus][1] = Q_loc [1]
                    else:
                        out [i][igaus] = Q_loc [comps_Q.index(comp)]

                offs = offs + gp_size_int
        return out

    # =====================================================
    def _get_rsl_for_sel_elements (self,ele_indices_1,solution_index,rsl_key,comp=None, \
                                   Armer_Wood_flag=False, fiber_flag=''):
    # =====================================================
        element = self.mesh.get_element(ele_indices_1 [0])
        ngaus = len(element.xsiGP)

        solution_indices = [solution_index]
        size = len(ele_indices_1)

        if rsl_key == Shell_EleResults.sNforce_key or rsl_key == Shell_EleResults.sMoment_key:
            if comp == None:
                if Armer_Wood_flag:
                    out = np.zeros (size*ngaus*2).reshape (size,ngaus,2)
                else:
                    out = np.zeros (size*ngaus*3).reshape (size,ngaus,3)
            else:
                out = np.zeros (size*ngaus).reshape(size,ngaus)
        else:
            if comp == None:
                out = np.zeros(size * ngaus * 2).reshape(size, ngaus, 2)
            else:
                out = np.zeros(size * ngaus).reshape(size, ngaus)

        for i,ele_index_1 in enumerate(ele_indices_1):
            ret = self._get_rsl_time_history (ele_index_1,solution_indices,rsl_key,comp,\
                                              Armer_Wood_flag,fiber_flag)
            for igaus in range (ret.shape[1]):
                if ret.ndim == 2:
                    out[i][igaus] = ret[0][igaus]
                else:
                    for k in range (ret.shape [2]):#loop over components
                        out [i][igaus][k] = ret [0][igaus][k]
        return out

    #public functions
    # =====================================================
    def get_thickness (self,ele_index_1,solution_index):
    # =====================================================
        element = self.mesh.get_element(ele_index_1)
        ngaus   = len(element.xsiGP)
        if not isinstance(element, Element_shell_Q4) and not isinstance(element, Element_shell_SHQ4) \
                and not isinstance(element, Element_shell_T3):
            return None

        sel_elements     = [self.mesh.get_element(ele_index_1)]
        solution_indices = [solution_index]
        res_all = self.get_element_results_ex (sel_elements, solution_indices)

        rsl_key = Shell_EleResults.thick_key
        col, n_comps = self.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, '')

        out = np.zeros (ngaus)

        gp_size_int = self.rcf.give_one_gp_size(sel_elements[0].get_group())

        for igaus in range(ngaus):
            res = res_all[igaus]
            offs = 2  # skip element index and gauss point index
            out.append (res [col])

        return out

    # =====================================================
    def  get_N_forces_vec_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        #ele_index_1 : element index (starts from 1)
        #solution_indices : indices of time instances (see HistoryOfExecution class)
        #return numpy array of size (n_time_instances,ngaus,3) containing Nxx,Nyy,Nxy forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sNforce_key,None)

    # =====================================================
    def  get_N_AW_forces_vec_time_history (self,ele_index_1,solution_indices,fiber_flag):
    # =====================================================
        #ele_index_1 : element index (starts from 1)
        #solution_indices : indices of time instances (see HistoryOfExecution class)
        #return numpy array of size (n_time_instances,ngaus,2) containing Nxx,Nyy forces
        #fiber flag is '+' or '-'
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sNforce_key,None,True,fiber_flag)

    # =====================================================
    def  get_Q_forces_vec_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        #ele_index_1 : element index (starts from 1)
        #solution_indices : indices of time instances (see HistoryOfExecution class)
        #return numpy array of size (n_time_instances,ngaus,2) containing Qx,Qy forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sQforce_key,None)

    # =====================================================
    def  get_N_force_time_history (self,ele_index_1,solution_indices,comp):
    # =====================================================
        # ele_index_1 : element index (starts from 1)
        # solution_indices : indices of time instances (see HistoryOfExecution class)
        # comp : one of 'XX', 'YY', 'XY'
        # return numpy array of size (n_time_instances,ngaus) containing given force component
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sNforce_key,comp)

    # =====================================================
    def  get_N_AW_force_time_history (self,ele_index_1,solution_indices,comp,fiber_flag):
    # =====================================================
        # ele_index_1 : element index (starts from 1)
        # solution_indices : indices of time instances (see HistoryOfExecution class)
        # comp : one of 'XX', 'YY'
        # return numpy array of size (n_time_instances,ngaus) containing given force component
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sNforce_key,comp,True,fiber_flag)

    # =====================================================
    def  get_Q_force_time_history (self,ele_index_1,solution_indices,comp):
    # =====================================================
        # ele_index_1 : element index (starts from 1)
        # solution_indices : indices of time instances (see HistoryOfExecution class)
        # comp : one of 'X', 'Y'
        # return numpy array of size (n_time_instances,ngaus) containing given force component
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sQforce_key,comp)

    # =====================================================
    def  get_N_forces_vec_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # return numpy array of size (n_elements,ngaus,3) NXX,NYY,NXY
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sNforce_key,None)

    # =====================================================
    def  get_N_AW_forces_vec_for_sel_elements (self,ele_indices_1,solution_index,fiber_flag):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # return numpy array of size (n_elements,ngaus,2) NXX,NYY
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sNforce_key,None,\
                                               True,fiber_flag)

    # =====================================================
    def  get_Q_forces_vec_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # return numpy array of size (n_elements,ngaus,2)
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sQforce_key,None)


    # =====================================================
    def  get_N_force_for_sel_elements (self,ele_indices_1,solution_index,comp):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        #comp is one of 'XX','YY','XY
        # return numpy array of size (n_elements,ngaus)
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sNforce_key,comp)

    # =====================================================
    def  get_N_AW_force_for_sel_elements (self,ele_indices_1,solution_index,comp,fiber_flag):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # comp is one of 'XX','YY','XY
        # return numpy array of size (n_elements,ngaus)
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sNforce_key,comp,\
                                               True,fiber_flag)

    # =====================================================
    def  get_Q_force_for_sel_elements (self,ele_indices_1,solution_index,comp):
    # =====================================================
        # ele_indices_1 : element indices (start from 1)
        # solution_index: index of time instance (see HistoryOfExecution class)
        # comp is one of 'X','Y'
        # return numpy array of size (n_elements,ngaus)
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sQforce_key,comp)


    #here moments - same explanation as for forces (see above)

    # =====================================================
    def  get_moments_vec_time_history (self,ele_index_1,solution_indices):
    # =====================================================
        #as for N forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sNforce_key,None)

    # =====================================================
    def  get_moments_AW_vec_time_history (self,ele_index_1,solution_indices,fiber_flag):
    # =====================================================
        #as for N forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sNforce_key,None,\
                                           True,fiber_flag)

    # =====================================================
    def  get_moment_time_history (self,ele_index_1,solution_indices,comp):
    # =====================================================
        #as for N forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sMoment_key,comp)

    # =====================================================
    def  get_moment_AW_time_history (self,ele_index_1,solution_indices,comp,fiber_flag):
    # =====================================================
        # as for N forces
        return self._get_rsl_time_history (ele_index_1,solution_indices,Shell_EleResults.sMoment_key,comp,\
                                           True,fiber_flag)

    # =====================================================
    def  get_moments_vec_for_sel_elements (self,ele_indices_1,solution_index):
    # =====================================================
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sMoment_key,None)

    # =====================================================
    def  get_moments_AW_vec_for_sel_elements (self,ele_indices_1,solution_index,fiber_flag):
    # =====================================================
        # as for N forces
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sMoment_key,None,\
                                               True,fiber_flag)

    # =====================================================
    def  get_moment_for_sel_elements (self,ele_indices_1,solution_index,comp):
    # =====================================================
        # as for N forces
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sMoment_key,comp)

    # =====================================================
    def  get_moment_AW_for_sel_elements (self,ele_indices_1,solution_index,comp,fiber_flag):
    # =====================================================
        # as for N forces
        return self._get_rsl_for_sel_elements (ele_indices_1,solution_index,Shell_EleResults.sMoment_key,comp,\
                                               True,fiber_flag)


    # =====================================================
    def _get_envelope_for_sel_elements(self, ele_indices_1, solution_indices, comp_index_ex, \
                                       Armer_Wood_flag=False, fiber_flag=''):
    # =====================================================
        # returns 2 matrices in which each row consists of [NXX,NYY,NXY,MXX,MYY,MXY,QX,QY] generalized force components
        # or 2 matrices in which each row consists of [NXX,NYY,MXX,MYY,QX,QY] generalized force components for Armer Wood
        element = self.mesh.get_element(ele_indices_1[0])
        ngaus = len(element.xsiGP)

        if Armer_Wood_flag:
            ncols = len(Shell_EleResults.comps_all_AW)
        else:
            ncols = len(Shell_EleResults.comps_all)

        out_max = np.zeros(len(ele_indices_1) * ngaus * ncols).reshape(len(ele_indices_1),ngaus, ncols)
        out_min = np.zeros(len(ele_indices_1) * ngaus * ncols).reshape(len(ele_indices_1),ngaus, ncols)
        aux = np.zeros(len(ele_indices_1) * ngaus * ncols).reshape(len(ele_indices_1),ngaus, ncols)
        vecComps = Shell_EleResults.comps_vec
        tnsCompos= Shell_EleResults.comps_tns

        for solution_index in solution_indices:
            if Armer_Wood_flag:
                N = self.get_N_AW_forces_vec_for_sel_elements  (ele_indices_1, solution_index,fiber_flag)
                M = self.get_moments_AW_vec_for_sel_elements   (ele_indices_1, solution_index,fiber_flag)
                size_out = M.shape [2]
            else:
                N = self.get_N_forces_vec_for_sel_elements  (ele_indices_1, solution_index)
                M = self.get_moments_vec_for_sel_elements   (ele_indices_1, solution_index)
                size_out = M.shape [2]
            Q = self.get_Q_forces_vec_for_sel_elements  (ele_indices_1, solution_index)

            offs = 0
            for i in range(N.shape[0]):#element
                for j in range (N.shape [1]):#gauss point
                    for k in range (N.shape [2]):#component
                        aux[i][j][k+offs] = N[i][j][k]

            offs = offs + size_out
            for i in range(M.shape[0]):
                for j in range (M.shape [1]):
                    for k in range (M.shape [2]):
                        aux[i][j][k+offs] = M[i][j][k]

            offs = offs + size_out
            for i in range(Q.shape[0]):
                for j in range (Q.shape [1]):
                    for k in range (Q.shape [2]):
                        aux[i][j][k+offs] = Q[i][j][k]


            for i in range(aux.shape [0]):
                for j in range (aux.shape [1]):
                    if aux[i][j][comp_index_ex] > out_max[i][j][comp_index_ex]:
                        out_max[i][j][comp_index_ex] = aux[i][j][comp_index_ex]
                        for k in range(aux.shape [2]):
                            if k != comp_index_ex:
                                out_max[i][j][k] = aux[i][j][k]

                    if aux[i][j][comp_index_ex] < out_min[i][j][comp_index_ex]:
                        out_min [i][j][comp_index_ex] = aux[i][j][comp_index_ex]
                        for k in range(aux.shape [2]):
                            if k != comp_index_ex:
                                out_min [i][j][k] = aux[i][j][k]

        return out_min, out_max


    # =====================================================
    def  get_N_force_envelope_for_sel_elements (self,ele_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,ngaus,8) in which each row consists of [NXX,NYY,NXY,MXX,MYY,MXY,QX,QY] generalized force components
        # it keeps envelope for the given force component and corresponding other forces/moments
        # first matrix is for min envelope and second for the max one
        #tnsComps = Shell_EleResults.comps_tns
        comp_index = Shell_EleResults.comps_all.index ('N'+comp) #tnsComps.index (comp)
        return self._get_envelope_for_sel_elements( ele_indices_1, solution_indices, comp_index)

    # =====================================================
    def  get_N_force_AW_envelope_for_sel_elements (self,ele_indices_1,solution_indices,comp,fiber_flag):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,ngaus,6) in which each row consists of [NXX,NYY,MXX,MYY,QX,QY] generalized force components
        # it keeps envelope for the given force component and corresponding other forces/moments
        # first matrix is for min envelope and second for the max one
        # fiber_flag is one of '+','-'
        #tnsComps = Shell_EleResults.comps_tns
        #comp_index = tnsComps.index (comp)
        comp_index = Shell_EleResults.comps_all_AW.index('N' + comp)
        return self._get_envelope_for_sel_elements( ele_indices_1, solution_indices, comp_index,True,fiber_flag)


    # =====================================================
    def  get_moment_envelope_for_sel_elements (self,ele_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,ngaus,8) in which each row consists of [NXX,NYY,NXY,MXX,MYY,MXY,QX,QY] generalized force components
        # it keeps envelope for the given force component and corresponding other forces/moments
        # first matrix is for min envelope and second for the max one
        #tnsComps = Shell_EleResults.comps_tns
        #comp_index = tnsComps.index (comp) + 3
        comp_index = Shell_EleResults.comps_all.index('M' + comp)
        return self._get_envelope_for_sel_elements(ele_indices_1, solution_indices, comp_index)

    # =====================================================
    def  get_moment_AW_envelope_for_sel_elements (self,ele_indices_1,solution_indices,comp,fiber_flag):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,ngaus,6) in which each row consists of [NXX,NYY,MXX,MYY,QX,QY] generalized force components
        # it keeps envelope for the given force component and corresponding other forces/moments
        # first matrix is for min envelope and second for the max one
        # fiber_flag is one of '+','-'
        #tnsComps = Shell_EleResults.comps_tns
        #comp_index = tnsComps.index (comp) + 3
        comp_index = Shell_EleResults.comps_all_AW.index('M' + comp)
        return self._get_envelope_for_sel_elements(ele_indices_1, solution_indices, comp_index, True,fiber_flag)


    # =====================================================
    def  get_Q_force_envelope_for_sel_elements (self,ele_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,ngaus,8) in which each row consists of [NXX,NYY,NXY,MXX,MYY,MXY,QX,QY] generalized force components
        # it keeps envelope for the given force component and corresponding forces/moments
        # first matrix is for min envelope and second for the max one
        #vecComps = Shell_EleResults.comps_vec
        #comp_index = vecComps.index (comp)
        comp_index = Shell_EleResults.comps_all.index('Q' + comp)
        return self._get_envelope_for_sel_elements( ele_indices_1, solution_indices, comp_index)

    # =====================================================
    def  get_Q_force_AW_envelope_for_sel_elements (self,ele_indices_1,solution_indices,comp,fiber_flag):
    # =====================================================
        # returns 2 matrices of dimension (N_ele,ngaus,6) in which each row consists of  [NXX,NYY,MXX,MYY,QX,QY] generalized force components
        # it keeps envelope for the given force component and corresponding forces/moments
        # first matrix is for min envelope and second for the max one
        comp_index = Shell_EleResults.comps_all_AW.index('Q' + comp)
        return self._get_envelope_for_sel_elements(ele_indices_1, solution_indices, comp_index, \
                                                   True, fiber_flag)


    # =====================================================
    def get_MNQ_time_history_for_fict_uniform_beam (self,sel_elements,solution_indices,\
                                                    x_axis_pt,x_axis_vec):
    # =====================================================
        versor = np.array(x_axis_vec)
        length = np.linalg.norm(versor)
        versor = versor / length
        pt = np.array(x_axis_pt)

        #first set local bases in shell elements
        for element in sel_elements:
            element.set_transf_matrices (x_axis_vec)
        sorted_sel_ele, proj_measures = self.mesh.sort_sel_elements_by_dist_along_dir (sel_elements,\
                                                                      x_axis_pt,x_axis_vec,True)
        #print proj_measures
        indices_of_sorted_sel_ele = self.mesh.get_ele_indices_for_sel_elements (sorted_sel_ele)

        ele_size_in_dir_Z = np.zeros (len(indices_of_sorted_sel_ele))
        ele_size_in_dir_X = np.zeros (len(indices_of_sorted_sel_ele))
        for i,element in enumerate(sorted_sel_ele):
            #set element size in perpendicular direction
            if i==3:
                jaja=1
            T_GL = element.get_T_GL ()
            ez = T_GL [2,:]
            ey = np.cross (ez,versor)
            #in beam context Z (shell) is swapped with Y (beam)
            ele_size_in_dir_Z [i] = element.get_size_along_dir (ey)
            ele_size_in_dir_X [i] = element.get_size_along_dir (versor)
            #print i,ele_size_in_dir_X [i],ele_size_in_dir_Z [i]

        #find largest projection measure distance
        if len(proj_measures) == 1:
            n_ele_in_perp_dir = 1
            n_ele_along_axis  = 1
        else:
            last_proj_meas = proj_measures [0]
            max_proj_meas  = 0.0
            for i in range (len(proj_measures)-1):
                next_proj_meas = proj_measures [i+1]
                max_proj_meas  = max (max_proj_meas,abs(next_proj_meas-last_proj_meas))
                last_proj_meas = next_proj_meas

            #detect how many elements are present in direction perpendicular to the axis
            check = [1]
            last_proj_meas = proj_measures [0]

            for i in range (len(proj_measures)-1):
                next_proj_meas = proj_measures [i+1]
                dist = abs (next_proj_meas-last_proj_meas)
                if dist/max_proj_meas < 1.0e-2:
                    check [-1] = check [-1] + 1
                else:
                    check.append (1)
                last_proj_meas = next_proj_meas

            for i in range (len(check)-1):
                if check [i+1] - check [i] != 0:
                    return False, [], [], [], [], []

            n_ele_in_perp_dir = check [0]
            n_ele_along_axis  = len(check)

        Forces = np.zeros(n_ele_along_axis * len(solution_indices) *6).reshape(n_ele_along_axis, len(solution_indices),6)
        x_tab = np.zeros(n_ele_along_axis)
        x_size= np.zeros(n_ele_along_axis)
        axis_pt = np.zeros(3)

        count = 0
        for j in range(n_ele_along_axis):
            x_tab   [j] = proj_measures     [count]
            x_size  [j] = ele_size_in_dir_X [count]
            count = count + n_ele_in_perp_dir

        for i, solution_index in enumerate (solution_indices):
            #these results are for all elements
            N = self.get_N_forces_vec_for_sel_elements (indices_of_sorted_sel_ele,solution_index)
            #print "N",N
            Q = self.get_Q_forces_vec_for_sel_elements (indices_of_sorted_sel_ele,solution_index)
            M = self.get_moments_vec_for_sel_elements  (indices_of_sorted_sel_ele,solution_index)
            ngaus = N.shape [1]
            count = 0
            for j in range (n_ele_along_axis):
                axis_pt[:] = pt[:] + versor[:] * proj_measures[count]
                for k in range (n_ele_in_perp_dir):
                    dL = ele_size_in_dir_Z [count] / float(ngaus)
                    ele_center = sorted_sel_ele [count].get_center()
                    dX = ele_center - axis_pt
                    dx = T_GL.dot (dX)
                    for igaus in range (ngaus):
                        #NX
                        Forces [j,i,0] = Forces [j,i,0] + N [count,igaus,0] * dL # n_xx dL
                        #QY
                        Forces [j,i,1] = Forces [j,i,1] + Q [count,igaus,0] * dL # q_x dL
                        #QZ
                        Forces [j,i,2] = Forces [j,i,2] + N [count,igaus,2] * dL #n_xy dL
                        #MX
                        Forces [j,i,3] = Forces [j,i,3] + M [count,igaus,2] * dL + \
                                                          Q [count,igaus,0] * dx [1] * dL
                        #MY
                        Forces [j,i,4] = Forces [j,i,5] + N [count,igaus,0] * dx [1] * dL
                        #MZ
                        Forces [j,i,5] = Forces [j,i,5] + M [count,igaus,0] * dL # m_xx dL
                    count = count + 1

        return True, x_tab, x_size, Forces



#=====================================================
def main():
#=====================================================
    #project = "d:/vxx_zsoil/TEMPLATES/tests/TEST-MNT-RECOVER-SHELLS"
    from .C_ContinuumResults import plot_figures
    from .C_HistoryOfExecution import HistoryOfExecution
    from .C_Mesh import Mesh
    from .C_Rcf_info import RCF_info
    project = "d:/vxx_zsoil/TEMPLATES/tests/test-solw-profile-plus-pile-case-1"  # TEST-MNT-RECOVER-CONT-2"
    mesh = Mesh (project)
    rcf  = RCF_info (project)
    his  = HistoryOfExecution (project)

    zoom = [[], [], []]
    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]

    shell_rsl = Shell_EleResults(mesh,his,rcf)
    mat_filter = [3]
    sel_elements = mesh.get_list_of_elements(Element.GROUP_SHELL, \
                                              0.0, False, mat_filter=mat_filter, \
                                              ele_class_filter=[], zoom_filter=zoom_filter)


    solution_indices_aux = his.give_converged_time_solutions()
    solution_indices = [solution_indices_aux [-1]]


    ti = np.zeros(len(solution_indices))
    for i, solution_index in enumerate(solution_indices):
        row = his.data[solution_index]
        ti[i] = row[his.DATA_TIME]


    dir = np.array ([0.0, -1.0, 0.0])
    #pt  = np.array ([0.0,  4.0, 0.0])
    pt  = np.array ([0.0,  0.0, 0.0])

    statusOK, x_tab, x_size, Forces = shell_rsl.get_MNQ_time_history_for_fict_uniform_beam( \
                                     sel_elements, solution_indices, \
                                     pt, dir)
    Forces = Forces * 2.0
    #for plot purposes we want to improve bending moments and plot constant Q,N forces
    x_tab_ex, Forces_ex = postprocess_MNQ_time_history_for_fict_uniform_beam(Forces, x_tab, x_size)
    x_tab_ex = x_tab_ex * (-1.0)

    smooth_jumps_in_stress_resultant_time_history (Forces_ex,5)

    plot_figures(x_tab_ex, Forces_ex[:, :, 0], Forces_ex[:, :, 1], Forces_ex[:, :, 5], "Forces-test.png", ti)
    pylab.show ()
    #print "z_tab",x_tab
    #print "z_size",x_size
    #print "NX",Forces [:,:,0]
    #print "QY",Forces [:,:,1]
    #print "MZ",Forces [:,:,5]

    #function get_MNQ_envelopes_for_fict_uniform_beam is statiuc and kept in C_ContinuumEleResults
    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam (Forces_ex, 0)#NX envelope
    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam (Forces_ex, 1)#QY envelope
    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam (Forces_ex, 5)#MZ envelope

if __name__ == '__main__':
    main()
