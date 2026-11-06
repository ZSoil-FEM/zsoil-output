"""Provide the legacy execution-history reader retained for compatibility investigations.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

This variant of :class:`HistoryOfExecution` is kept for reproducing and diagnosing
historical parsing behaviour. New code should normally use :mod:`C_HistoryOfExecution`.
"""

# Name:        HistoryOfExecution
# Purpose:
#
# Author:      Andrzej Truty
#
# Created:     19-08-2015
# Copyright:   (c) ZACE 2015
# Licence:     <ZACE licence>
# &'(3I5,e12.5,i5,e14.6,1x,i5,i5,i5,e14.6,e14.6,1x,a21,i5,1x,e14.6,1
# &x,e14.6,1x,i5,1x,i5,1x,e14.6,6(1x,e14.6),4(1x,e14.6))'


class HistoryOfExecution ():
    def __init__(self, project):
        """Initialize the history of execution object for the supplied project data.

        :param project: Path or filename prefix of the ZSoil project.
        :type project: str
        """
        self.my_project = project
        self.DATA_DRIVER_TYPE = 0  # driver type index
        self.DATA_NITER = 1  # number of iterations
        self.DATA_CONV_STATUS = 2  # convergence status(-1/0)
        self.DATA_SF = 3  # current safety factor
        self.DATA_PLOT_FLAG = 4  # plot storage flag (0/1)
        self.DATA_TIME = 5  # current time
        self.DATA_STEP_COUNT = 6  # step counter
        self.DATA_PUSH_FLAG = 7  # pushover_flag
        self.DATA_NR_EIGEN_MODES = 8  # NoEigenModes
        self.DATA_PUSH_LAMBDA = 9  # pushover load factor
        self.DATA_PUSH_CTRL_DISPL = 10  # pushover control disp
        self.DATA_PUSH_LABEL = 11  # user label for pushover data
        self.DATA_ARC_LENGTH_FLAG = 12  # arc_length_flag
        self.DATA_ARC_LENGTH_U = 13  # U norm
        self.DATA_ARC_LENGTH_LAMBDA = 14  # load factor
        self.DATA_ARC_LENGTH_NL_SOLVER = 15  # nl_solver
        self.DATA_MAX_ABS_ITER = 16  # max_abs_iter
        self.DATA_AMPLF_CONV = 17  # ampl_cnv_factor
        self.DATA_SCALAR_MAX = 18
        self.DATA_CONV_RHS_NORMS = 18
        self.DATA_CONV_ENE_NORMS = 19
        self.DATA_MAX = 20

        self.dtypes = []
        for i in range(self.DATA_SCALAR_MAX):
            self.dtypes.append(type(1))

        self.dtypes[self.DATA_DRIVER_TYPE] = type(1)
        self.dtypes[self.DATA_NITER] = type(1)
        self.dtypes[self.DATA_CONV_STATUS] = type(1)
        self.dtypes[self.DATA_SF] = type(1.0)
        self.dtypes[self.DATA_PLOT_FLAG] = type(1)
        self.dtypes[self.DATA_TIME] = type(1.0)
        self.dtypes[self.DATA_STEP_COUNT] = type(1)
        self.dtypes[self.DATA_PUSH_FLAG] = type(1)
        self.dtypes[self.DATA_NR_EIGEN_MODES] = type(1)
        self.dtypes[self.DATA_PUSH_LAMBDA] = type(1.0)
        self.dtypes[self.DATA_PUSH_CTRL_DISPL] = type(1.0)
        self.dtypes[self.DATA_PUSH_LABEL] = type("")
        self.dtypes[self.DATA_ARC_LENGTH_FLAG] = type(1)
        self.dtypes[self.DATA_ARC_LENGTH_U] = type(1.0)
        self.dtypes[self.DATA_ARC_LENGTH_LAMBDA] = type(1.0)
        self.dtypes[self.DATA_ARC_LENGTH_NL_SOLVER] = type(1)
        self.dtypes[self.DATA_MAX_ABS_ITER] = type(1)
        self.dtypes[self.DATA_AMPLF_CONV] = type(1.0)

        self.CONVERGED_STATUS = -1

        self.debug = False

        self.data = []

        if self.debug:
            print("")
            print("reading history of execution->" + project)
        self.instanciate(project+".his")
        if self.debug:
            print("done..")

    def instanciate(self, filename):
        """Load this object's data from the associated ZSoil project input.

        :param filename: Filename used by this operation.
        :type filename: str
        """
        try:
            f = open(filename, "rt")
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            return
        except:
            print("I/O error({0}): {1}")
            raise
            return

        f.readline()

        for line in f:
            if len(line) <= 0:
                f.close()
                return

            aux = []
            line1 = line[0: 90]
            line2 = line[91:111]
            line3 = line[112:]
            texts = line1.split()
            for i in range(self.DATA_PUSH_CTRL_DISPL+1):
                if self.dtypes[i] == type(1):
                    aux.append(int(texts[i]))
                elif self.dtypes[i] == type(1.0):
                    aux.append(float(texts[i]))
                else:
                    aux.append(texts[i])

            aux.append(line2.strip())

            texts = line3.split()
            n_read = 0
            for i in range(self.DATA_SCALAR_MAX-self.DATA_PUSH_LABEL-1):
                k = self.DATA_PUSH_LABEL+1+i
                if self.dtypes[k] == type(1):
                    aux.append(int(texts[i]))
                elif self.dtypes[k] == type(1.0):
                    aux.append(float(texts[i]))
                else:
                    aux.append(texts[i])
                n_read = n_read + 1

            n_left = len(texts) - n_read
            offs = n_read
            aux_rhs_norm = []
            for i in range(6):
                aux_rhs_norm.append(float(texts[offs+i]))
            offs = offs + 6
            n_left = n_left - 6
            aux_ene_norms = []
            for i in range(n_left):
                aux_ene_norms.append(float(texts[offs+i]))
            aux.append(aux_rhs_norm)
            aux.append(aux_ene_norms)

            self.data.append(aux)

        f.close()

    def last_step_converged(self):
        """Return whether the most recently recorded solution step converged.

        :return: Whether the most recently recorded solution step converged.
        :rtype: bool
        """
        size = len(self.data)
        row = self.data[size-1]
        ret = row[self.DATA_CONV_STATUS]
        if ret == self.CONVERGED_STATUS:
            return True
        else:
            return False

    def give_converged_solutions(self, time_filter=[-1.0, 1.0e99], give_cs_calc_flag=False):
        """Return indices of converged, plottable solutions within the requested time range.

        :param time_filter: Inclusive minimum and maximum physical times.
        :type time_filter: Sequence[float]
        :param give_cs_calc_flag: Whether to enable give cs calc.
        :type give_cs_calc_flag: bool
        :return: The indices of converged, plottable solutions within the requested time range.
        :rtype: Union[list, tuple[list, list]]
        """
        cs = []
        cs_calc = []
        count = 0
        for i in range(len(self.data)):
            row = self.data[i]
            if row[self.DATA_CONV_STATUS] == self.CONVERGED_STATUS and row[self.DATA_PLOT_FLAG] == 1:  # converged step
                if row[self.DATA_TIME] >= time_filter[0] and row[self.DATA_TIME] <= time_filter[1]:
                    cs.append(count)
                    cs_calc.append(i)
            if row[self.DATA_PLOT_FLAG] == 1:  # converged step
                count = count + 1

        if give_cs_calc_flag:
            return cs, cs_calc
        else:
            return cs

    def give_converged_time_solutions(self, time_filter=[-1.0, 1.0e99], give_cs_calc_flag=False):
        """Return converged time-step solution indices within the requested time range.

        :param time_filter: Inclusive minimum and maximum physical times.
        :type time_filter: Sequence[float]
        :param give_cs_calc_flag: Whether to enable give cs calc.
        :type give_cs_calc_flag: bool
        :return: The converged time-step solution indices within the requested time range.
        :rtype: Union[list, tuple[list, list]]
        """
        cs = []
        cs_calc = []
        count = 0
        for i in range(len(self.data)):
            row = self.data[i]

            if row[self.DATA_CONV_STATUS] == self.CONVERGED_STATUS and row[self.DATA_PLOT_FLAG] == 1:  # converged step
                if abs(row[self.DATA_SF]) < 1.0e-12:
                    if row[self.DATA_PUSH_LABEL] == "":
                        if row[self.DATA_ARC_LENGTH_FLAG] == 0:
                            if row[self.DATA_TIME] >= time_filter[0] and row[self.DATA_TIME] <= time_filter[1]:
                                cs_calc.append(i)
                                cs.append(count)

            if row[self.DATA_PLOT_FLAG] == 1:  # converged step
                count = count + 1

        if give_cs_calc_flag:
            return cs, cs_calc
        else:
            return cs

    def give_sols_for_time(self, time):
        """Return the two converged solution indices that bracket the requested time.

        :param time: Physical time at which the value or state is requested.
        :type time: float
        :return: The two converged solution indices that bracket the requested time.
        :rtype: list
        """
        good_sols, good_sols_calc = self.give_converged_time_solutions(give_cs_calc_flag=True)
        for i in range(len(good_sols)-1):
            sol_t1_index = good_sols[i]
            sol_t2_index = good_sols[i+1]
            row1 = self.data[good_sols_calc[i]]
            row2 = self.data[good_sols_calc[i+1]]
            t1 = row1[self.DATA_TIME]
            t2 = row2[self.DATA_TIME]
            if time >= t1 and time <= t2:
                aux = [sol_t1_index, sol_t2_index]
                return aux
        return []

    def give_closest_sol_index_for_time(self, time):
        """Return the converged solution index closest to the requested time.

        :param time: Physical time at which the value or state is requested.
        :type time: float
        :return: The converged solution index closest to the requested time.
        :rtype: int
        """
        good_sols, good_sols_calc = self.give_converged_time_solutions(give_cs_calc_flag=True)
        closest = -1
        dt_min = 1.0e38

        for i in range(len(good_sols)):
            row = self.data[good_sols_calc[i]]
            t = row[self.DATA_TIME]
            if abs(t-time) < dt_min:
                dt_min = abs(t-time)
                closest = good_sols[i]
        return closest

    def give_converged_sol_index_for_time(self, time, good_sols_in):
        """Return a converged solution index matching the requested time, or -1.

        :param time: Physical time at which the value or state is requested.
        :type time: float
        :param good_sols_in: Good sols in used by this operation.
        :type good_sols_in: Sequence[int]
        :return: A converged solution index matching the requested time, or -1.
        :rtype: int
        """
        #        if len(good_sols_in)>0:
        #            good_sols = good_sols_in
        #        else:
        #            good_sols = self.give_converged_time_solutions()
        good_sols, good_sols_calc = self.give_converged_time_solutions(give_cs_calc_flag=True)

        closest = -1
        dt_min = 1.0e38

        for i in range(len(good_sols)):
            row = self.data[good_sols_calc[i]]
            t = row[self.DATA_TIME]
            if abs(t-time) < dt_min:
                dt_min = abs(t-time)
                closest = good_sols[i]

        if dt_min < 0.001:
            return closest
        return -1


def main():
    """Run this module's example or command-line workflow."""
    his = HistoryOfExecution('D:/v14_ztuser_vs2012/INP/bugs/menard-16-08-2015/nasyp Gdansk-at1-dense')
    ret = his.give_converged_time_solutions()
    print(ret)


if __name__ == '__main__':
    main()
