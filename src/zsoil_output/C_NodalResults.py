import os
import numpy
import numpy as np
from   math import *


#=====================================================
class NodalResults ():
#=====================================================

    #====================================================
    def __init__ (self,mesh,his,rcf):
    #=====================================================
        self.mesh = mesh
        self.his  = his
        self.rcf  = rcf

        self.dict_exts          = {"NODAL_RSL":".s00","NODAL_RES":".r00","NODAL_VEL":".v00","NODAL_ACC":".a00"}
        self.dict_node_size_int = {"NODAL_RSL":rcf.give_one_node_size(),"NODAL_RES":rcf.give_one_node_size(),"NODAL_VEL":rcf.give_one_node_size_v(),"NODAL_ACC":rcf.give_one_node_size_a()}

    #======================================
    def get_nodal_results (self,selected_nodes,solution_indices,dict_key):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        if len(selected_nodes) <= 0:
            return out

        rsl_file_ext = self.dict_exts [dict_key]
        node_size_int= self.dict_node_size_int[dict_key]

        f      = open (mesh.my_project+rsl_file_ext,"rb")
        node_size     = node_size_int * 4 #in bytes

        for node in selected_nodes:
            node_rsl = [node.index]
            for sol_index in solution_indices:
                seek_pos = sol_index * len(mesh.nodes) * node_size + (node.index-1) * node_size
                f.seek    (seek_pos,os.SEEK_SET)  # seek
                arr = numpy.fromfile (f,dtype=numpy.float32,count=node_size_int,sep="")  # read the data into numpy
                for i in range(len(arr)):
                    node_rsl.append ( arr [i] )
            out.append (node_rsl)
        f.close ()
        return out



    #======================================
    def _get_node_vectorial_rsl_ths (self,node_index_1,solution_indices,dict_key,rsl_key,comp=None):
    #======================================
        node = self.mesh.get_node(node_index_1)
        selected_nodes = [node]
        res_all = self.get_nodal_results(selected_nodes, solution_indices, dict_key)
        vec = res_all[0]  # only one node is considered here
        size = len(solution_indices)

        offs = 1 # as the first position is occupied by the node number in the vector
        node_size_int = self.dict_node_size_int [dict_key]

        if comp == None:
            col_1, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, 'X')
            col_2, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, 'Y')
            col_3, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, 'Z')
            out = np.zeros(size * 3).reshape (size,3)
        else:
            col_1, ncomps = self.rcf.give_col_for_nod_rsl(rsl_key, comp)
            col_2 = -1
            col_3 = -1
            out = np.zeros(size)

        for i in range(len(solution_indices)):
            if comp == None:
                if col_1 != -1:
                    out[i][0] = vec[col_1+offs]
                if col_2 != -1:
                    out[i][1] = vec[col_2+offs]
                if col_3 != -1:
                    out[i][2] = vec[col_3+offs]
            else:
                if col_1 != -1:
                    out[i] = vec[col_1+offs]

            offs = offs + node_size_int

        return out


    # ======================================
    def _get_node_scalar_rsl_ths(self, node_index_1, solution_indices, dict_key, rsl_key):
    # ======================================
        node = self.mesh.get_node(node_index_1)
        selected_nodes = [node]
        res_all = self.get_nodal_results(selected_nodes, solution_indices, dict_key )
        vec = res_all[0]  # only one node is considered here
        size = len(solution_indices)


        offs = 1  # as the first position is occupied by the node number in the vector
        node_size_int = self.dict_node_size_int [dict_key]

        col, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, '')
        out = np.zeros(size)

        for i in range(len(solution_indices)):
            if col != -1:
                out[i] = vec[col+offs]
            offs = offs + node_size_int

        return out


    # ======================================
    def _get_nodes_vectorial_rsl (self, node_indices_1, solution_index, dict_key, rsl_key, comp=None):
    # ======================================
        # get displacements for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # comp is one of 'X', 'Y', 'Z'; if None take full vector
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        solution_indices = [solution_index]
        if comp == None:
            ret = np.zeros(len(node_indices_1) * 3).reshape (len(node_indices_1),3)
        else:
            ret = np.zeros(len(node_indices_1))

        row = 0
        for node_index_1 in node_indices_1:
            out = self._get_node_vectorial_rsl_ths(node_index_1, solution_indices, dict_key, rsl_key, comp)
            ret[row] = out[0]
            row = row + 1

        return ret


    # ======================================
    def _get_nodes_scalar_rsl (self, node_indices_1, solution_index,dict_key, rsl_key):
    # ======================================
        # get pore pressures for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        solution_indices = [solution_index]
        ret = np.zeros(len(node_indices_1))

        row = 0
        for node_index_1 in node_indices_1:
            out = self._get_node_scalar_rsl_ths (node_index_1, solution_indices,dict_key, rsl_key)
            ret[row] = out[0]
            row = row + 1

        return ret

    # ---------------------------------
    # public  functions
    #functions to be used by the user
    # ---------------------------------

    #======================================
    def get_node_displacements_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get displacement for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_RSL','DISP_TRA',comp)

    #======================================
    def get_node_reaction_forces_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get displacement for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_RES','DISP_TRA',comp)

    #======================================
    def get_node_veloc_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get velocity for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_VEL','DISP_TRA',comp)

    #======================================
    def get_node_accel_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get accelerations for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_ACC','DISP_TRA',comp)

    #======================================
    def get_nodes_displacements (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get displacements for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_RSL','DISP_TRA',comp)

    #======================================
    def get_nodes_reaction_forces (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get displacements for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_RES','DISP_TRA',comp)

    #======================================
    def get_nodes_veloc (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get velocities for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_VEL','DISP_TRA',comp)


    #======================================
    def get_nodes_accel (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get accelerations for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_ACC','DISP_TRA',comp)



    # ======================================
    def get_node_rotations_time_history (self, node_index_1, solution_indices, comp=None):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','DISP_ROT', comp)

    # ======================================
    def get_node_reaction_moments_time_history (self, node_index_1, solution_indices, comp=None):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths(node_index_1, solution_indices,'NODAL_RES','DISP_ROT', comp)


    # ======================================
    def get_nodes_rotations (self, node_indices_1, solution_index, comp=None):
    # ======================================
        # get rotations for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # comp is one of 'X', 'Y', 'Z'; if None take full vector
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl(node_indices_1, solution_index,'NODAL_RSL','DISP_ROT', comp)

    # ======================================
    def get_nodes_reaction_moments (self, node_indices_1, solution_index, comp=None):
    # ======================================
        # get rotations for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # comp is one of 'X', 'Y', 'Z'; if None take full vector
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl(node_indices_1, solution_index,'NODAL_RES','DISP_ROT', comp)



    # ======================================
    def get_node_pore_pressures_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','PPRESS')


    # ======================================
    def get_nodes_pore_pressures(self, node_indices_1, solution_index):
    # ======================================
        # get pore pressures for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return  self._get_nodes_scalar_rsl (node_indices_1, solution_index,'NODAL_RSL','PPRESS')

    # ======================================
    def get_node_pressure_heads_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','PRES_HEAD')


    # ======================================
    def get_nodes_pressure_heads(self, node_indices_1, solution_index):
    # ======================================
        # get pressure heads for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_nodes_scalar_rsl(node_indices_1, solution_index,'NODAL_RSL','PRES_HEAD')

    # ======================================
    def get_node_temperatures_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','TEMP')


    # ======================================
    def get_nodes_temperatures(self, node_indices_1, solution_index):
    # ======================================
        # get temperatures for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_nodes_scalar_rsl(node_indices_1, solution_index,'NODAL_RSL','TEMP')

    # ======================================
    def get_node_humidities_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','HUMID')


    # ======================================
    def get_nodes_humidities(self, node_indices_1, solution_index):
    # ======================================
        # get humidities for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_nodes_scalar_rsl (node_indices_1, solution_index,'NODAL_RSL','HUMID')


    # =====================================================
    def _get_envelope_vec_for_sel_nodes (self, node_indices_1, solution_indices, dict_key, rsl_key, comp_index_4):
    # =====================================================
        # returns 2 matrices in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        out_max = np.zeros(len(node_indices_1) * 4).reshape(len(node_indices_1),4)
        out_min = np.zeros(len(node_indices_1) * 4).reshape(len(node_indices_1),4)
        aux     = np.zeros(len(node_indices_1) * 4).reshape(len(node_indices_1),4)
        vecComps = ['X', 'Y', 'Z','ABS']

        vec = np.zeros (4)
        for solution_index in solution_indices:
            vecx = self._get_nodes_vectorial_rsl(node_indices_1, solution_index, dict_key, rsl_key)
            for i in range(len(node_indices_1)):
                tmp = 0.0
                for k in range (3):
                    vec [k] = vecx [i][k]
                    tmp = tmp + vec [k] * vec [k]
                tmp = sqrt(tmp)
                vec [3] = tmp
                aux [i] = vec

            for i in range(len(node_indices_1)):
                if aux[i][comp_index_4] > out_max[i][comp_index_4]:
                    out_max[i][comp_index_4] = aux[i][comp_index_4]
                    for k in range(4):
                        if k != comp_index_4:
                            out_max[i][k] = aux[i][k]

                if aux[i][comp_index_4] < out_min[i][comp_index_4]:
                    out_min [i][comp_index_4] = aux[i][comp_index_4]
                    for k in range(4):
                        if k != comp_index_4:
                            out_min [i][k] = aux[i][k]

        return out_min, out_max

    # =====================================================
    def _get_envelope_scalar_for_sel_nodes (self, node_indices_1, solution_indices, dict_key, rsl_key):
    # =====================================================
        # returns 2 matrices in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        out_max = np.zeros(len(node_indices_1) )
        out_min = np.zeros(len(node_indices_1) )

        for solution_index in solution_indices:
            aux = self._get_nodes_scalar_rsl(node_indices_1, solution_index, dict_key, rsl_key)

            for i in range(len(node_indices_1)):
                if aux[i] > out_max[i]:
                    out_max [i] = aux[i]
                if aux[i] < out_min[i]:
                    out_min [i] = aux[i]
        return out_min, out_max


    # =====================================================
    def  get_displ_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RSL',"DISP_TRA",comp_index)

    # =====================================================
    def  get_reaction_forces_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RES',"DISP_TRA",comp_index)

    # =====================================================
    def  get_rot_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [RX,RY,RZ,R-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RSL',"DISP_ROT",comp_index)

    # =====================================================
    def  get_reaction_moments_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [RX,RY,RZ,R-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RES',"DISP_ROT",comp_index)

    # =====================================================
    def  get_pressure_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','PPRESS')

    # =====================================================
    def  get_pressure_head_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','PRES_HEAD')

    # =====================================================
    def  get_temperature_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','TEMP')

    # =====================================================
    def  get_humidity_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','HUMID')



    # =====================================================
    def  get_veloc_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_VEL',"DISP_TRA",comp_index)

    # =====================================================
    def  get_accel_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_ACC',"DISP_TRA",comp_index)

