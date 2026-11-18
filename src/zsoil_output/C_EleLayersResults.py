import os
import numpy


#=====================================================
class EleLayersResults ():
#=====================================================

    #====================================================
    def __init__ (self,mesh,his,rcf,rcf_lay):
    #=====================================================
        self.mesh = mesh
        self.his  = his
        self.rcf  = rcf
        self.rcf_lay = rcf_lay

    #======================================
    def get_layers_results (self,group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]]):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf
        rcf_lay = self.rcf_lay

        out = []

        sel_elements = mesh.get_list_of_elements (group,0.0,False,mat_filter,[],zoom_filter)
        if len(sel_elements) <= 0:
            return out

        ext    = sel_elements [0].rsl_lay_ext
        if ext == None:
            return out

        f      = open (mesh.my_project+ext,"rb")
        sum_of_layers_in_rsl_file = mesh.get_sum_of_layers (group)
        lay_size_int = rcf_lay.give_one_lay_size (group)
        lay_size= lay_size_int * 4 #in bytes

        for e in sel_elements:
            ngaus = len(e.xsiGP)
            nlayers_all, nlayers_reinf = e.get_nlayers()
            for igaus in range (ngaus):
                for layer in range (nlayers_all):
                    lay_rsl = [e.index,igaus+1,layer+1]
                    for sol_index in solution_indices:
                        seek_pos = sol_index * sum_of_layers_in_rsl_file * lay_size
                        seek_pos = seek_pos + e.seek_in_lay_rsl
                        seek_pos = seek_pos + igaus * nlayers_all * lay_size
                        seek_pos = seek_pos + layer * lay_size
                        f.seek    (seek_pos,os.SEEK_SET)  # seek
                        arr = numpy.fromfile (f,dtype=numpy.float32,count=lay_size_int,sep="")  # read the data into numpy
                        #arr = arr.reshape ( (len(mesh.elems),gp_size),order="C")
                        for i in range(len(arr)):
                            lay_rsl.append ( arr [i] )
                    out.append (lay_rsl)
        f.close ()
        return out


    #======================================
    def _get_layers_results_ex (self, group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]], core_layers_flag=False):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf
        rcf_lay = self.rcf_lay

        out = []

        sel_elements = mesh.get_list_of_elements (group,0.0,False,mat_filter,[],zoom_filter)
        if len(sel_elements) <= 0:
            return out

        ext    = sel_elements [0].rsl_lay_ext
        if ext == None:
            return out

        f      = open (mesh.my_project+ext,"rb")
        sum_of_layers_in_rsl_file = mesh.get_sum_of_layers (group)
        lay_size_int = rcf_lay.give_one_lay_size (group)
        lay_size= lay_size_int * 4 #in bytes

        for e in sel_elements:
            ngaus = len(e.xsiGP)
            nlayers_all, nlayers_reinf = e.get_nlayers()
            nlayers_core = nlayers_all - nlayers_reinf
            if core_layers_flag:
                layer_beg = 0
                layer_end = nlayers_core
            else:
                layer_beg = nlayers_core
                layer_end = nlayers_all

            for igaus in range (ngaus):
                for layer in range (layer_beg,layer_end):
                    lay_rsl = [e.index,igaus+1,layer+1]
                    for sol_index in solution_indices:
                        seek_pos = sol_index * sum_of_layers_in_rsl_file * lay_size
                        seek_pos = seek_pos + e.seek_in_lay_rsl
                        seek_pos = seek_pos + igaus * nlayers_all * lay_size
                        seek_pos = seek_pos + layer * lay_size
                        f.seek    (seek_pos,os.SEEK_SET)  # seek
                        arr = numpy.fromfile (f,dtype=numpy.float32,count=lay_size_int,sep="")  # read the data into numpy
                        #arr = arr.reshape ( (len(mesh.elems),gp_size),order="C")
                        for i in range(len(arr)):
                            lay_rsl.append ( arr [i] )
                    out.append (lay_rsl)
        f.close ()
        return out



    #======================================
    def get_core_layers_results (self,group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]]):
    #======================================
        out = self._get_layers_results_ex (group,solution_indices,mat_filter,zoom_filter,core_layers_flag=True)
        return out

    #======================================
    def get_reinforcement_layers_results (self,group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]]):
    #======================================
        out = self._get_layers_results_ex (group,solution_indices,mat_filter,zoom_filter,core_layers_flag=False)
        return out


