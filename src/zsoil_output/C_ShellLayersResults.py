import os
import numpy as np

from .C_EleResults         import *
from .C_Element_shell_Q4   import *
from .C_Element_shell_T3   import *
from .C_Element_shell_SHQ4 import *
#from C_PlasticCodes       import *



#class to get results for shell elements

#=====================================================
class ShellLayers_Results ():
#=====================================================

    COMPS_TNS = ['XX', 'YY', 'XY', 'ZZ', 'XZ', 'YZ']
    COMPS_VEC = ['X', 'Y', 'Z']
    STRESS_KEY  = 'STRESSES'
    STRAIN_KEY  = 'STRAINS'
    POS_KEY     = 'POSITION'
    PL_CODE_KEY = 'PLA_CODE'
    STR_LEV_KEY = 'STR_LEVEL'
    DAMAGE_KEY  = 'DAMAGE'


    # =====================================================
    def __init__(self,mesh,his,rcf, rcf_lay ):
    # =====================================================
        self.mesh = mesh
        self.his = his
        self.rcf = rcf
        self.rcf_lay = rcf_lay
        self.comp_dict = {'XX':0, 'YY':1, 'XY':2,'ZZ':3,'XZ':4,'YZ':5}
        self.stress_map_dict = {'XX': (0,0), 'YY': (1,1), 'XY': (0,1), 'ZZ': (2,2), 'XZ': (0,2), 'YZ': (1,2)}
        # self.stress_map = np.zeros (9,dtype=int).reshape(3,3)
        # self.stress_map [0,0] = self.comp_dict ('XX')
        # self.stress_map [0,1] = self.comp_dict ('XY')
        # self.stress_map [0,2] = self.comp_dict ('XZ')
        # self.stress_map [1,0] = self.comp_dict('XY')
        # self.stress_map [1,1] = self.comp_dict('YY')
        # self.stress_map [1,2] = self.comp_dict('YZ')
        # self.stress_map [2,0] = self.comp_dict('XZ')
        # self.stress_map [2,1] = self.comp_dict('YZ')
        # self.stress_map [2,2] = self.comp_dict('ZZ')





    #======================================
    def _get_layers_results (self, sel_elements, solution_indices, layer_index_1 = None):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf
        rcf_lay = self.rcf_lay

        out = []
        if len(sel_elements) <= 0:
            return out

        ext = sel_elements [0].rsl_lay_ext
        if ext == None:
            return out

        f      = open (mesh.my_project+ext,"rb")
        sum_of_layers_in_rsl_file = mesh.get_sum_of_layers (Element.GROUP_SHELL)
        lay_size_int = rcf_lay.give_one_lay_size (Element.GROUP_SHELL)
        lay_size= lay_size_int * 4 #in bytes


        for e in sel_elements:
            e_single = []
            ngaus = len(e.xsiGP)
            nlayers_all, nlayers_reinf = e.get_nlayers()
            gp_all = []
            for igaus in range (ngaus):
                gp_single = []
                if layer_index_1 == None:
                    size_layers = nlayers_all
                else:
                    size_layers = 1

                for layer_ex in range (size_layers):

                    if layer_index_1 == None:
                        layer = layer_ex
                    else:
                        layer = layer_index_1 - 1

                    single_lay_rsl = [] #lay_rsl = [e.index,igaus+1,layer+1]
                    for sol_index in solution_indices:
                        seek_pos = sol_index * sum_of_layers_in_rsl_file * lay_size
                        seek_pos = seek_pos + e.seek_in_lay_rsl
                        seek_pos = seek_pos + igaus * nlayers_all * lay_size
                        seek_pos = seek_pos + layer * lay_size
                        f.seek    (seek_pos,os.SEEK_SET)  # seek
                        arr = numpy.fromfile (f,dtype=numpy.float32,count=lay_size_int,sep="")  # read the data into numpy
                        #arr = arr.reshape ( (len(mesh.elems),gp_size),order="C")
                        # for i in range(len(arr)):
                        single_lay_rsl.append ( arr )
                    gp_single.append (single_lay_rsl)
                gp_all.append (gp_single)
            e_single.append (gp_all)
            out.append (e_single)

        f.close ()
        return out


    # =====================================================
    def get_rsl_time_history (self,element, layer_index_1, solution_indices, rsl_key, comp=None):
    # =====================================================
        #returns numpy array
        # if comp == None:
        #     out = np.zeros(size(solution_indices) * ngaus * 6).reshape(size, ngaus, 6)
        # else:
        #     out = np.zeros(size(solution_indices) * ngaus).reshape(size, ngaus)

        ngaus   = len(element.xsiGP)
        if not isinstance(element,Element_shell_Q4) and not isinstance(element,Element_shell_SHQ4) \
                and not isinstance(element,Element_shell_T3):
            return None

        # T_TNS = element.get_T_TNS ()
        T_VEC = element.get_T_VEC ()
        T_VEC_t = np.transpose (T_VEC)

        sel_elements = [element]
        res_all = self._get_layers_results (sel_elements, solution_indices, layer_index_1)
        size = len(solution_indices)

        if rsl_key == self.STRESS_KEY:
            col_xx, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XX')
            col_yy, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YY')
            col_xy, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XY')
            col_xz, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XZ')
            col_yz, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YZ')

            if comp == None:
                out = np.zeros(size*ngaus*6).reshape (size,ngaus,6)
            else:
                out = np.zeros(size*ngaus).reshape (size,ngaus)

            gp_rsl_all = res_all [0]
            for igaus in range (ngaus):
                gp_rsl_single = gp_rsl_all [igaus]
                lay_rsl_single = gp_rsl_single[0]
                for sol_index in range (size):
                    whole_rsl_record = lay_rsl_single [sol_index]
                    stress = np.zeros (9).reshape(3,3)
                    stress [0,0] = whole_rsl_record [col_xx]
                    stress [0,1] = whole_rsl_record [col_xy]
                    stress [0,2] = whole_rsl_record [col_xz]
                    stress [1,0] = stress [0,1]
                    stress [1,1] = whole_rsl_record [col_yy]
                    stress [1,2] = whole_rsl_record [col_yz]
                    stress [2,0] = stress [0,2]
                    stress [2,1] = stress [1,2]
                    stress [2,2] = 0.0

                    stress_transformed = T_VEC.dot(stress.dot (T_VEC_t))
                    if comp == None:
                        for i, comp_key in enumerate(self.COMPS_TNS):
                            ij_tuple = self.stress_map_dict [comp_key]
                            out [sol_index,igaus,i] =  stress_transformed [ij_tuple [0],ij_tuple[1]]
                    else:
                        ij_tuple = self.stress_map_dict [comp]
                        out[sol_index, igaus] = stress_transformed [ij_tuple [0],ij_tuple[1]]

        elif rsl_key == self.STRAIN_KEY:

            col_xx, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XX')
            col_yy, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YY')
            col_xy, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XY')
            col_xz, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'XZ')
            col_yz, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, 'YZ')

            if comp == None:
                out = np.zeros(size * ngaus * 6).reshape(size, ngaus, 6)
            else:
                out = np.zeros(size * ngaus * 1).reshape(size, ngaus, 1)

            gp_rsl_all = res_all[0]
            for igaus in range(ngaus):
                gp_rsl_single = gp_rsl_all[igaus]
                lay_rsl_single = gp_rsl_single[0]
                for sol_index in range(size):
                    whole_rsl_record = lay_rsl_single[sol_index]
                    strain = np.zeros(9).reshape(3, 3)
                    strain[0, 0] = whole_rsl_record[col_xx]
                    strain[0, 1] = whole_rsl_record[col_xy] / 2.0
                    strain[0, 2] = whole_rsl_record[col_xz] / 2.0
                    strain[1, 0] = strain[0, 1]
                    strain[1, 1] = whole_rsl_record[col_yy]
                    strain[1, 2] = whole_rsl_record[col_yz] / 2.0
                    strain[2, 0] = strain[0, 2]
                    strain[2, 1] = strain[1, 2]
                    strain[2, 2] = 0.0

                    strain_transformed = T_VEC.dot(strain.dot(T_VEC_t))
                    for k in range (3):
                        for r in range (3):
                            if k != r:
                                strain_transformed [k,r] = strain_transformed [k,r] * 2.0

                    if comp == None:
                        for i, comp_key in enumerate(self.COMPS_TNS):
                            ij_tuple = self.stress_map_dict[comp_key]
                            out[sol_index, igaus, i] = strain_transformed[ij_tuple[0], ij_tuple[1]]
                    else:
                        ij_tuple = self.stress_map_dict[comp]
                        out[sol_index, igaus, 0] = strain_transformed[ij_tuple[0], ij_tuple[1]]

        else: #scalar results

            col, n_comps = self.rcf_lay.give_col_for_ele_rsl(sel_elements[0].get_group(), rsl_key, None)
            out = np.zeros(size * ngaus * 1).reshape(size, ngaus, 1)

            gp_rsl_all = res_all[0]
            for igaus in range(ngaus):
                gp_rsl_single = gp_rsl_all[igaus]
                lay_rsl_single = gp_rsl_single[0]
                for sol_index in range(size):
                    whole_rsl_record = lay_rsl_single[sol_index]
                    out[sol_index, igaus, 0] = whole_rsl_record [col]

        return out

    # =====================================================
    def get_rsl_for_sel_elements (self, sel_elements, layer_index_1, solution_index, rsl_key, comp=None):
    # =====================================================
    #returns list of numpy arrays of size [ngaus,ncomp]
        solution_indices = [solution_index]
        table = []
        for element in sel_elements:
            out = self.get_rsl_time_history(element, layer_index_1, solution_indices, rsl_key, comp)
            out.reshape (out.shape[1],out.shape[2])
            table.append (out)
        return table


