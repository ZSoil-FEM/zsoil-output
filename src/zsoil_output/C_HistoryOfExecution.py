#-------------------------------------------------------------------------------
# Name:        HistoryOfExecution
# Purpose:
#
# Author:      Andrzej Truty
#
# Created:     19-08-2015
# Copyright:   (c) ZACE 2015
# Licence:     <ZACE licence>
#-------------------------------------------------------------------------------

#old format
##     &'(3I5,e12.5,i5,e14.6,1x,i5,i5,i5,e14.6,e14.6,1x,a21,i5,1x,e14.6,1
##     &x,e14.6,1x,i5,1x,i5,1x,e14.6,6(1x,e14.6),4(1x,e14.6))'

#new format
     #  character*118,parameter:: HisFileFormat_1607_HIS =
     # &'(3I5,e12.5,i5,e20.12,1x,i5,i5,i5,e14.6,e14.6,1x,a21,i5,1x,e14.6,1
     # &x,e14.6,1x,i5,1x,i5,1x,e14.6,6(1x,e14.6),4(1x,e14.6))'

#  character*118,parameter:: HisFileFormat_2409_HIS =
# &'(3I5,e12.5,i5,e20.12,1x,i10,i5,i5,e14.6,e14.6,1x,a21,i5,1x,e14.6,1
# &x,e14.6,1x,i5,1x,i5,1x,e14.6,6(1x,e14.6),4(1x,e14.6))'

import os
from .ZSoil_version import *

#=====================================================
class HistoryOfExecution ():
#=====================================================

    #=====================================================
    def __init__ (self,project):
    #=====================================================
        self.my_project = project

        zv = ZSoil_version (project)
        self.his_version = 0
        self.dat_version = zv.get_dat_version ()

        self.DATA_DRIVER_TYPE = 0 #driver type index
        self.DATA_NITER       = 1 #number of iterations
        self.DATA_CONV_STATUS = 2 #convergence status(-1/0)
        self.DATA_SF          = 3 #current safety factor
        self.DATA_PLOT_FLAG   = 4 #plot storage flag (0/1)
        self.DATA_TIME        = 5 #current time
        self.DATA_STEP_COUNT  = 6 #step counter
        self.DATA_PUSH_FLAG   = 7 #pushover_flag
        self.DATA_NR_EIGEN_MODES = 8 #NoEigenModes
        self.DATA_PUSH_LAMBDA = 9 #pushover load factor
        self.DATA_PUSH_CTRL_DISPL = 10 #pushover control disp
        self.DATA_PUSH_LABEL  = 11#user label for pushover data
        self.DATA_ARC_LENGTH_FLAG = 12 #arc_length_flag
        self.DATA_ARC_LENGTH_U= 13#U norm
        self.DATA_ARC_LENGTH_LAMBDA = 14#load factor
        self.DATA_ARC_LENGTH_NL_SOLVER=15#nl_solver
        self.DATA_MAX_ABS_ITER= 16#max_abs_iter
        self.DATA_AMPLF_CONV = 17#ampl_cnv_factor
        self.DATA_SCALAR_MAX = 18
        self.DATA_CONV_RHS_NORMS = 18 # 6 norms F/M/Qf/Qth/Qhum/Theta
        self.DATA_CONV_ENE_NORMS = 19 # 4 norms Esolid/Efluid/Etherm/Ehumi
        self.DATA_NR_NSTD_LRECS  = 20
        self.DATA_MAX = 21

        self.dtypes = []
        for i in range (self.DATA_SCALAR_MAX):
            self.dtypes.append (type(1))

        self.dtypes[self.DATA_DRIVER_TYPE    ] = type (1)
        self.dtypes[self.DATA_NITER          ] = type (1)
        self.dtypes[self.DATA_CONV_STATUS    ] = type (1)
        self.dtypes[self.DATA_SF             ] = type (1.0)
        self.dtypes[self.DATA_PLOT_FLAG      ] = type (1)
        self.dtypes[self.DATA_TIME           ] = type (1.0)
        self.dtypes[self.DATA_STEP_COUNT     ] = type (1)
        self.dtypes[self.DATA_PUSH_FLAG      ] = type (1)
        self.dtypes[self.DATA_NR_EIGEN_MODES ] = type (1)
        self.dtypes[self.DATA_PUSH_LAMBDA    ] = type (1.0)
        self.dtypes[self.DATA_PUSH_CTRL_DISPL] = type (1.0)
        self.dtypes[self.DATA_PUSH_LABEL     ] = type("")
        self.dtypes[self.DATA_ARC_LENGTH_FLAG] = type (1)
        self.dtypes[self.DATA_ARC_LENGTH_U   ] = type (1.0)
        self.dtypes[self.DATA_ARC_LENGTH_LAMBDA] = type (1.0)
        self.dtypes[self.DATA_ARC_LENGTH_NL_SOLVER] = type (1)
        self.dtypes[self.DATA_MAX_ABS_ITER   ] = type (1)
        self.dtypes[self.DATA_AMPLF_CONV     ] = type (1.0)

        self.CONVERGED_STATUS = -1

        self.debug = False

        self.data = []


        if self.dat_version == None or self.dat_version < 0:
            return

        if self.debug:
            print("")
            print("reading history of execution->"+ project)
        self.instanciate (project+".his")
        if self.debug:
            print("done..")


    #=====================================================
    def instanciate (self,filename):
    #=====================================================

        try:
            f = open (filename,"rt")
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            return
        except:
            print("I/O error({0}): {1}")
            raise
            return

#AT his version since 18.03

        line = f.readline()
        aux  = line.split()
        self.his_version = int(aux[0])

        # 

        # new_format = True
##OK
        # if self.dat_version <= 1607 or (self.dat_version > 1690 and self.dat_version <= 1696):
        #     new_format = False

        for line in f:

            # #AT 3.07.2018 - in case on developer versions
            # if len(line) > 325:
            #     new_format = True

            if len(line) <= 0:
                f.close()
                return

            aux = []

            if self.dat_version <= 1607 or (self.dat_version > 1690 and self.dat_version <= 1696):
                line1 = line [0  : 90]
                line2 = line [91 :111]
                line3 = line [112:   ]
            elif self.dat_version < 2409:
                line1 = line [0  : 96]
                line2 = line [102 :117]
                line3 = line [118:   ]
            else: # extended step counter to 10 since 2409 added nr of nstd lrecs (needed in v25)
                line1 = line [0  : 101]
                line2 = line [102 :122]
                line3 = line [123:   ]

            # if new_format:
            #     line1 = line [0  : 96]
            #     line2 = line [97 :117]
            #     line3 = line [118:   ]
            # else:
            #     line1 = line [0  : 90]
            #     line2 = line [91 :111]
            #     line3 = line [112:   ]

            texts = line1.split()
            for i in range (self.DATA_PUSH_CTRL_DISPL+1):
                if self.dtypes [i] == type(1):
                    aux.append (int(texts [i]))
                elif self.dtypes [i] == type (1.0):
                    aux.append (float(texts [i]))
                else:
                    aux.append (texts [i])

            aux.append (line2.strip())

            texts = line3.split()
            n_read= 0
            for i in range (self.DATA_SCALAR_MAX-self.DATA_PUSH_LABEL-1):
                k = self.DATA_PUSH_LABEL+1+i
                if self.dtypes [k] == type(1):
                    aux.append (int(texts [i]))
                elif self.dtypes [k] == type (1.0):
                    aux.append (float(texts [i]))
                else:
                    aux.append (texts [i])
                n_read = n_read + 1

            n_left = len(texts) - n_read
            offs = n_read
            aux_rhs_norm = []
            for i in range (6):
                aux_rhs_norm.append (float(texts[offs+i]))
            offs = offs + 6
            n_left = n_left - 6
            aux_ene_norms = []
            for i in range (4):
                aux_ene_norms.append (float(texts[offs+i]))
            offs = offs + 4
            n_left = n_left - 4
            aux.append (aux_rhs_norm )
            aux.append (aux_ene_norms)
            if n_left > 0:
                aux.append (int(texts[offs+0]))
            else:
                aux.append(0)
            #AT 28.11.2017 BUG omitt steps for which results are not stored
            if aux [self.DATA_PLOT_FLAG] == 1:
                self.data.append (aux)

        f.close ()

    #=====================================================
    def last_step_converged (self):
    #=====================================================
        size = len(self.data)
        row = self.data [size-1]
        ret = row [self.DATA_CONV_STATUS]
        if ret == self.CONVERGED_STATUS:
            return True
        else:
            return False


    #=====================================================
    def give_converged_solutions (self,time_filter=[-1.0,1.0e99]):
    #=====================================================
        cs    = []
        for i in range(len(self.data)):
            row = self.data [i]
            if row [self.DATA_CONV_STATUS] == self.CONVERGED_STATUS and row [self.DATA_PLOT_FLAG] == 1: ##converged step
                if row [self.DATA_TIME] >= time_filter [0] and row [self.DATA_TIME] <= time_filter [1]:
                    cs.append (i)

        return cs


    #=====================================================
    def give_solutions_with_plot_status (self,time_filter=[-1.0,1.0e99]):
    #=====================================================
        cs    = []
        for i in range(len(self.data)):
            row = self.data [i]
            if row [self.DATA_PLOT_FLAG] == 1: ##converged step
                if row [self.DATA_TIME] >= time_filter [0] and row [self.DATA_TIME] <= time_filter [1]:
                    cs.append (i)
        return cs


    #=====================================================
    def give_converged_time_solutions (self,time_filter = [-1.0,1.0e99]):
    #=====================================================
        cs    = []
        for i in range(len(self.data)):
            row = self.data [i]
			
            if row [self.DATA_CONV_STATUS] == self.CONVERGED_STATUS and row [self.DATA_PLOT_FLAG] == 1: ##converged step
                if abs(row [self.DATA_SF]) < 1.0e-12:
                    if row [self.DATA_PUSH_LABEL] == "":
                        if row [self.DATA_ARC_LENGTH_FLAG] == 0:
                            if row [self.DATA_TIME] >= time_filter [0] and row [self.DATA_TIME] <= time_filter [1]:
                                cs.append (i)

        return cs

    #=====================================================
    def give_sols_for_time (self,time):
    #=====================================================
        good_sols = self.give_converged_time_solutions()
        for i in range (len(good_sols)-1):
            sol_t1_index = good_sols [i]
            sol_t2_index = good_sols [i+1]
            row1 = self.data [sol_t1_index]
            row2 = self.data [sol_t2_index]
            t1 = row1 [self.DATA_TIME]
            t2 = row2 [self.DATA_TIME]
            if time >= t1 and time <= t2:
                aux = [sol_t1_index,sol_t2_index]
                return aux
        return []

    #=====================================================
    def give_closest_sol_index_for_time (self,time):
    #=====================================================
        good_sols = self.give_converged_time_solutions()
        closest   = -1
        dt_min    = 1.0e38

        for i in range (len(good_sols)):
            row = self.data [good_sols [i]]
            t = row [self.DATA_TIME]
            if abs(t-time) < dt_min:
                dt_min = abs(t-time)
                closest = good_sols [i]
        return closest

    #=====================================================
    def give_converged_sol_index_for_time (self,time,good_sols_in):
    #=====================================================
        if len(good_sols_in)>0:
            good_sols = good_sols_in
        else:
            good_sols = self.give_converged_time_solutions()

        closest   = -1
        dt_min    = 1.0e38

        for i in range (len(good_sols)):
            row = self.data [good_sols [i]]
            t = row [self.DATA_TIME]
            if abs(t-time) < dt_min:
                dt_min = abs(t-time)
                closest = good_sols [i]

        if dt_min<0.001:
            return closest
        return -1

    #=====================================================
    def give_converged_SF_for_time (self,time):
    #=====================================================
        maxSF = 0.0
        time_filter = [time-0.0001,time+0.0001]
        for i in range(len(self.data)):
            row = self.data [i]
			
            if row [self.DATA_CONV_STATUS] == self.CONVERGED_STATUS and row [self.DATA_PLOT_FLAG] == 1: ##converged step
                if row [self.DATA_TIME] >= time_filter [0] and row [self.DATA_TIME] <= time_filter [1]:
                    curSF = abs(row [self.DATA_SF])
                    if  curSF > maxSF:
                        maxSF = curSF

        return maxSF
        
    #=====================================================
    def give_converged_SF_solutions_for_time (self,time):
    #=====================================================
        out = []
        out_SF = []
        time_filter = [time-0.0001,time+0.0001]
        for i in range(len(self.data)):
            row = self.data [i]			
            if row [self.DATA_CONV_STATUS] == self.CONVERGED_STATUS and row [self.DATA_PLOT_FLAG] == 1: ##converged step
                if row [self.DATA_TIME] >= time_filter [0] and row [self.DATA_TIME] <= time_filter [1]:
                    if abs(row [self.DATA_SF]) > 1.0e-12:
                        currSF = abs(row [self.DATA_SF])
                        out.append (i)
                        out_SF.append (currSF)
        return out, out_SF
        
        
    #==================================================================
    def get_data_value_for_sol_instances (self,list_of_time_instances,par_index):
    #==================================================================
        out = []
        for solution_index in list_of_time_instances:
            row = self.data [solution_index]
            out.append (row [par_index])
        return out

#=====================================================
def main():
#=====================================================
    project = 'd:\\ZSOIL_work\\v24\\#inp\\Benchmarks\\SSH-ADVECTION-BENCH-DOMENICO-HE-4M-A'
    his  = HistoryOfExecution (project)
    ret = his.give_converged_solutions ()
    print(his.dat_version)
    print(ret)
    
if __name__ == '__main__':
    main()
