import os
import numpy
import numpy as np
from .C_Element import *


#=====================================================
class EleResults ():
#=====================================================

    BUFFER_F_HANDLE    = 0
    BUFFER_SUM_ALL_GP  = 1
    BUFFER_GP_SIZE_INT = 2
    BUFFER_MAX         = 3

    #====================================================
    def __init__ (self,mesh,his,rcf):
    #=====================================================
        self.mesh = mesh
        self.his  = his
        self.rcf  = rcf

        self.group_buffers = []
        for i in range (Element.GROUP_MAX):
            buffer = []
            for j in range (self.BUFFER_MAX):
                buffer.append (None)
            self.group_buffers.append (buffer)


    #======================================
    def get_element_results_ex (self,sel_elements,solution_indices):
    #======================================
	# all elements must belong to the same group

        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        if len(sel_elements) <= 0:
            return out

        ext    = sel_elements [0].rsl_ext
        group  = sel_elements [0].group
        group_index = sel_elements[0].group_dict_inv [group]

        if self.group_buffers [group_index][self.BUFFER_F_HANDLE] == None:
            f = open (mesh.my_project+ext,"rb")
            sum_of_all_gp_in_rsl_file = mesh.get_sum_of_gauss_points (group)
            gp_size_int = rcf.give_one_gp_size(group)
            self.group_buffers [group_index][self.BUFFER_F_HANDLE   ] = f
            self.group_buffers [group_index][self.BUFFER_SUM_ALL_GP ] = sum_of_all_gp_in_rsl_file
            self.group_buffers [group_index][self.BUFFER_GP_SIZE_INT] = gp_size_int
        else:
            f                         = self.group_buffers [group_index][self.BUFFER_F_HANDLE   ]
            sum_of_all_gp_in_rsl_file = self.group_buffers [group_index][self.BUFFER_SUM_ALL_GP ]
            gp_size_int               = self.group_buffers [group_index][self.BUFFER_GP_SIZE_INT]

        gp_size = gp_size_int * 4 #in bytes

        for e in sel_elements:
            ngaus = len(e.xsiGP)
            for igaus in range (ngaus):
                gauss_rsl = [e.index,igaus+1]
                for sol_index in solution_indices:
                    offs = sol_index * sum_of_all_gp_in_rsl_file * gp_size
                    seek_pos = offs + e.seek_in_rsl + igaus * gp_size # + col_rsl * 4
                    f.seek    (seek_pos,os.SEEK_SET)  # seek
                    arr = numpy.fromfile (f,dtype=numpy.float32,count=gp_size_int,sep="")  # read the data into numpy
                    for i in range(len(arr)):
                        gauss_rsl.append ( arr [i] )
                out.append (gauss_rsl)
        #do not close the file as it takes lot of time
        #f.close ()
        return out

    #======================================
    def get_element_results_ex (self,sel_elements,solution_indices):
    #======================================
	# all elements must belong to the same group

        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        if len(sel_elements) <= 0:
            return out

        ext    = sel_elements [0].rsl_ext
        group  = sel_elements [0].group
        group_index = sel_elements[0].group_dict_inv [group]

        if self.group_buffers [group_index][self.BUFFER_F_HANDLE] == None:
            f = open (mesh.my_project+ext,"rb")
            sum_of_all_gp_in_rsl_file = mesh.get_sum_of_gauss_points (group)
            gp_size_int = rcf.give_one_gp_size(group)
            self.group_buffers [group_index][self.BUFFER_F_HANDLE   ] = f
            self.group_buffers [group_index][self.BUFFER_SUM_ALL_GP ] = sum_of_all_gp_in_rsl_file
            self.group_buffers [group_index][self.BUFFER_GP_SIZE_INT] = gp_size_int
        else:
            f                         = self.group_buffers [group_index][self.BUFFER_F_HANDLE   ]
            sum_of_all_gp_in_rsl_file = self.group_buffers [group_index][self.BUFFER_SUM_ALL_GP ]
            gp_size_int               = self.group_buffers [group_index][self.BUFFER_GP_SIZE_INT]

        gp_size = gp_size_int * 4 #in bytes

        for e in sel_elements:
            ngaus = len(e.xsiGP)
            for igaus in range (ngaus):
                gauss_rsl = [e.index,igaus+1]
                for sol_index in solution_indices:
                    offs = sol_index * sum_of_all_gp_in_rsl_file * gp_size
                    seek_pos = offs + e.seek_in_rsl + igaus * gp_size # + col_rsl * 4
                    f.seek    (seek_pos,os.SEEK_SET)  # seek
                    arr = numpy.fromfile (f,dtype=numpy.float32,count=gp_size_int,sep="")  # read the data into numpy
                    for i in range(len(arr)):
                        gauss_rsl.append ( arr [i] )
                out.append (gauss_rsl)
        #do not close the file as it takes lot of time
        #f.close ()
        return out

    #======================================
    def get_element_results (self,group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]]):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        sel_elements = mesh.get_list_of_elements (group,0.0,False,mat_filter=mat_filter,\
                                                  ele_class_filter=[],zoom_filter=zoom_filter)
        print(len(sel_elements))
        if len(sel_elements) <= 0:
            return out
        out = self.get_element_results_ex (sel_elements,solution_indices)
        return out
