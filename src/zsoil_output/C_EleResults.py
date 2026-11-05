"""Provide shared binary-result reading operations for finite elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`EleResults` base locates integration-point records in ZSoil result files,
buffers open files, and returns raw values for selected elements and solution steps.
"""

import os
import numpy
import numpy as np
from .C_Element import *


class EleResults ():
    BUFFER_F_HANDLE = 0
    BUFFER_SUM_ALL_GP = 1
    BUFFER_GP_SIZE_INT = 2
    BUFFER_MAX = 3

    def __init__(self, mesh, his, rcf):
        """Initialize the element results object for the supplied project data.

        :param mesh: Mesh containing the model nodes, elements, and materials.
        :type mesh: Mesh
        :param his: Execution history that maps solution indices to times and states.
        :type his: HistoryOfExecution
        :param rcf: Result configuration describing stored nodal and element fields.
        :type rcf: RCF_info
        """
        self.mesh = mesh
        self.his = his
        self.rcf = rcf

        self.group_buffers = []
        for i in range(Element.GROUP_MAX):
            buffer = []
            for j in range(self.BUFFER_MAX):
                buffer.append(None)
            self.group_buffers.append(buffer)

    def get_element_results_ex(self, sel_elements, solution_indices):
        """Read element results for the requested groups, steps, and filters.

        .. note::
           All selected elements must belong to the same element group.

        :param sel_elements: Elements to process.
        :type sel_elements: Sequence[Element]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The element results for the requested groups, steps, and filters.
        :rtype: list
        """
        mesh = self.mesh
        his = self.his
        rcf = self.rcf

        out = []

        if len(sel_elements) <= 0:
            return out

        ext = sel_elements[0].rsl_ext
        group = sel_elements[0].group
        group_index = sel_elements[0].group_dict_inv[group]

        if self.group_buffers[group_index][self.BUFFER_F_HANDLE] == None:
            f = open(mesh.my_project+ext, "rb")
            sum_of_all_gp_in_rsl_file = mesh.get_sum_of_gauss_points(group)
            gp_size_int = rcf.give_one_gp_size(group)
            self.group_buffers[group_index][self.BUFFER_F_HANDLE] = f
            self.group_buffers[group_index][self.BUFFER_SUM_ALL_GP] = sum_of_all_gp_in_rsl_file
            self.group_buffers[group_index][self.BUFFER_GP_SIZE_INT] = gp_size_int
        else:
            f = self.group_buffers[group_index][self.BUFFER_F_HANDLE]
            sum_of_all_gp_in_rsl_file = self.group_buffers[group_index][self.BUFFER_SUM_ALL_GP]
            gp_size_int = self.group_buffers[group_index][self.BUFFER_GP_SIZE_INT]

        gp_size = gp_size_int * 4  # in bytes

        for e in sel_elements:
            ngaus = len(e.xsiGP)
            for igaus in range(ngaus):
                gauss_rsl = [e.index, igaus+1]
                for sol_index in solution_indices:
                    offs = sol_index * sum_of_all_gp_in_rsl_file * gp_size
                    seek_pos = offs + e.seek_in_rsl + igaus * gp_size  # + col_rsl * 4
                    f.seek(seek_pos, os.SEEK_SET)  # seek
                    arr = numpy.fromfile(f, dtype=numpy.float32, count=gp_size_int, sep="")  # read the data into numpy
                    for i in range(len(arr)):
                        gauss_rsl.append(arr[i])
                out.append(gauss_rsl)
        # do not close the file as it takes lot of time
        #f.close ()
        return out

    def get_element_results_ex(self, sel_elements, solution_indices):
        """Read element results for the requested groups, steps, and filters.

        .. note::
           All selected elements must belong to the same element group.

        :param sel_elements: Elements to process.
        :type sel_elements: Sequence[Element]
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :return: The element results for the requested groups, steps, and filters.
        :rtype: list
        """
        mesh = self.mesh
        his = self.his
        rcf = self.rcf

        out = []

        if len(sel_elements) <= 0:
            return out

        ext = sel_elements[0].rsl_ext
        group = sel_elements[0].group
        group_index = sel_elements[0].group_dict_inv[group]

        if self.group_buffers[group_index][self.BUFFER_F_HANDLE] == None:
            f = open(mesh.my_project+ext, "rb")
            sum_of_all_gp_in_rsl_file = mesh.get_sum_of_gauss_points(group)
            gp_size_int = rcf.give_one_gp_size(group)
            self.group_buffers[group_index][self.BUFFER_F_HANDLE] = f
            self.group_buffers[group_index][self.BUFFER_SUM_ALL_GP] = sum_of_all_gp_in_rsl_file
            self.group_buffers[group_index][self.BUFFER_GP_SIZE_INT] = gp_size_int
        else:
            f = self.group_buffers[group_index][self.BUFFER_F_HANDLE]
            sum_of_all_gp_in_rsl_file = self.group_buffers[group_index][self.BUFFER_SUM_ALL_GP]
            gp_size_int = self.group_buffers[group_index][self.BUFFER_GP_SIZE_INT]

        gp_size = gp_size_int * 4  # in bytes

        for e in sel_elements:
            ngaus = len(e.xsiGP)
            for igaus in range(ngaus):
                gauss_rsl = [e.index, igaus+1]
                for sol_index in solution_indices:
                    offs = sol_index * sum_of_all_gp_in_rsl_file * gp_size
                    seek_pos = offs + e.seek_in_rsl + igaus * gp_size  # + col_rsl * 4
                    f.seek(seek_pos, os.SEEK_SET)  # seek
                    arr = numpy.fromfile(f, dtype=numpy.float32, count=gp_size_int, sep="")  # read the data into numpy
                    for i in range(len(arr)):
                        gauss_rsl.append(arr[i])
                out.append(gauss_rsl)
        # do not close the file as it takes lot of time
        #f.close ()
        return out

    def get_element_results(self, group, solution_indices, mat_filter=[], zoom_filter=[0, [[], [], []]]):
        """Return element results for the requested groups and solution steps.

        :param group: Group used by this operation.
        :type group: str or int
        :param solution_indices: Solution indices defined by :class:`HistoryOfExecution`.
        :type solution_indices: Sequence[int]
        :param mat_filter: Material indices to include; an empty sequence includes all materials.
        :type mat_filter: Sequence[int]
        :param zoom_filter: Spatial filtering settings; the default disables filtering.
        :type zoom_filter: Sequence
        :return: The element results for the requested groups and solution steps.
        :rtype: list[list[float]]
        """
        mesh = self.mesh
        his = self.his
        rcf = self.rcf

        out = []

        sel_elements = mesh.get_list_of_elements(group, 0.0, False, mat_filter=mat_filter,
                                                 ele_class_filter=[], zoom_filter=zoom_filter)
        print(len(sel_elements))
        if len(sel_elements) <= 0:
            return out
        out = self.get_element_results_ex(sel_elements, solution_indices)
        return out
