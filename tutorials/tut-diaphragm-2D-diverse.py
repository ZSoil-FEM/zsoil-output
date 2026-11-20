#-------------------------------------------------------------------------------
# Name:        tut-diaphragm-2D-diverse
# Purpose:     extract and plot MNT+their envelopes, settlements behind the wall, wall deflections
#
# Author:      Andrzej Truty
#
# Created:     2020
# Copyright:   (c) ZACE Services Ltd 2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import matplotlib.pyplot  as pyplot
import matplotlib.pylab   as pylab
from math import *
import numpy as np
import xlwt
import sys,os

#miktex is needed to use the following package
from matplotlib import rc
# those who have installed miktex environment can deactivate the next statement to have nice
# Latex style axes descriptions, equations etc...
#rc('text',usetex=True)


#path to SDK
# Only needed if you want to use a local clone of zsoil_output instead of the pip-installed package:
# import sys
# sys.path.insert(0, '..\\src')   # the folder that CONTAINS zsoil_output\

from zsoil_output import (Mesh, HistoryOfExecution, RCF_info, Beam_EleResults,
                           NodalResults, Continuum_EleResults, Diagram, Contact_2D_Results, Element)

def CSV_reader (csv_name,separator,skip_nr_of_rows = 0):
    f = open (csv_name)
    data = []
    lines = f.readlines()
    echo = 0
    #if csv_name == "../csv/HB-circ-tunnel-axs-stress.csv":
        #echo = 1
    for i in range (0,len(lines)):
        if i+1 > skip_nr_of_rows:
            tmp  = lines [i]
            tmp = tmp.rstrip()
            tmp = tmp.replace(";;",";0.0;")
            if tmp[len(tmp)-1] == ";":
                tmp  = tmp + "0.0"
            #if echo:
                #print tmp
            row_data = []
            aux = tmp.split (separator);
            for j in range (0,len(aux)):
                a = aux [j]
                if echo == 1:
                    print (j,a)
                b = float (a)
                row_data.append (b)
            data.append(row_data)
    return data


# =================================================
def Wall_defl (data_file, time_instances, out_png):
# =================================================
    # goal:
    # this function will plot wall deflections at selected time instances
    #
    # data_file      - name of the project (without file extension !)
    # time_instances - list of real time instances (float values) for which results are to be plot
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    #create instance of diagram class which will contain plot of wall deflections
    dg = Diagram ()
    dg.set_figure_size ((6,12))
    markers = dg.get_markers (len(time_instances))
    colors  = dg.get_colors  (len(time_instances))

    # ------------------------------------------------------------------
    #next 3 calls must always be done to be ablr to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution (data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh (data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info (data_file)
    # ------------------------------------------------------------------

    # here we will plot wall deflections

    # in this project beam elements are used only for beams therefore we can select them by element group ID
    # on the output beam element objects are returned
    beams = mesh.get_group_of_elements_ex (Element.GROUP_BEAM)

    # let us collect all nodes that belong to given elements list
    # on the output node objects are returned
    beam_nodes = mesh.get_list_of_nodes (ele_filter = beams)

    # we do not know in which order nodes are returned therefore we need to sort them out
    # here we define the point
    point_top = [30.0,0.0 ,0.0]
    # here we define the vector to be used for sorting purposes
    dir_vec   = [0.0,-1.0,0.0]

    # now we will get sorted beam nodes in such a way that the first node will be on the top of the wall
    beam_nodes_sorted = mesh.sort_sel_nodes_by_dist_along_dir\
                             (beam_nodes, point_top, dir_vec)

    # now we will create a list of selected nodes indices (indexing starts from 1 !)
    beam_nodes_indices = []
    for k in range (len(beam_nodes_sorted)):
        node = beam_nodes_sorted [k]
        beam_nodes_indices.append (node.index)

    #now we will collect nodes vertical coordinates
    y_node_tab = []
    for node in beam_nodes_sorted:
        # get node coordinates
        xyz = node.get_xyz()
        # add y-coordinate to the list
        y_node_tab.append (xyz [1])

    # ---------------------------
    # now action : get ux profile from results
    # ---------------------------

    #create nodal results class instance
    nodal_rsl = NodalResults (mesh, his, rcf)

    # loop over time instances
    for i, time in enumerate(time_instances):

        #get closest solution index foir a given real time value
        solution_index = his.give_closest_sol_index_for_time (time)

        #get displacements profile for selected nodes
        #on the output we get a numpy array
        rsl = nodal_rsl.get_nodes_displacements (beam_nodes_indices, solution_index, 'X')

        #convert result from [m] to [mm]
        rsl = rsl * 1000.0

        # add a diagram to the plot
        label = str(time)
        dg.add_plot (rsl, y_node_tab, marker_=markers [i], markevery_=5, markersize_=6, \
                     markerfacecolor_=colors[i], markeredgecolor_='black', \
                     linestyle_='-',label_= label, color_=colors[i])

    dg.set_x_axes_tick_params (labelsize=16)
    dg.set_x_min   (-40.0)
    dg.set_x_max   (  0.0)
    dg.set_y_max   (0.0)
    dg.set_y_min   (-32.0)
    #dg.set_x_log_scale ()
    dg.set_x_label (r'$u_x$ [mm]')
    dg.set_y_label (r'$y$ [m]')
    dg.add_grid ()
    dg.add_legend()
    dg.save_fig(out_png)



# ==========================================================
def Wall_settl (zsoil_data, time_instances, out_png):
# ==========================================================
    # goal:
    # this function will plot settlements behind the wall at selected time instances
    #
    # data_file      - name of the project (without file extension !)
    # time_instances - list of real time instances (float values) for which results are to be plot
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    #create instance of diagram class which will contain plot of settlement behinf the wall
    dg = Diagram ()
    dg.set_figure_size ((12,6))
    markers = dg.get_markers (len(time_instances))
    colors  = dg.get_colors  (len(time_instances))

    # ------------------------------------------------------------------
    #next 3 calls must always be done to be ablr to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution (data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh (data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info (data_file)
    # ------------------------------------------------------------------

    # here we will plot settlements (we need to select nodes at the top starting from x=30m)
    # but at x=30 we have 3 nodes (left side of the wall, wall, right side of the wall)

    # here we select all continuum elements
    # on output we get list of continuum elements (objects)
    cont_elements = mesh.get_list_of_elements (Element.GROUP_CONTINUUM, time=0.0, \
                                               only_active = False)

    # here we select nodes that belong to the continuum elements but we put also a zoom in which nodes should be located
    # on output we get list of selected nodes (objects)
    x_wall = 30.0
    dx     = 0.01
    dy     = 0.001
    x_max  = 120.0
    # first parameter in zoom structure is 0 because we select nodes here not elements
    # first parameter in zoom filter can be
    # mesh.ZOOM_FILTER_OFF                    (= 0)
    # mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN  (= 1)
    # mesh.ZOOM_FILTER_ALL_NODES_MUST_BE_IN   (= 2)
    nodes_sel     = mesh.get_list_of_nodes (ele_filter = cont_elements,\
                                            zoom_filter = [0,[[x_wall-dx,x_max],[-dy,dy],[]]])



    # as we do not know the order we need to sort these nodes
    point_top = [30.0,0.0 ,0.0]
    dir_vec   = [1.0,0.0,0.0]

    #get list of sorted nodes
    nodes_sorted = mesh.sort_sel_nodes_by_dist_along_dir (nodes_sel, point_top, dir_vec)

    #let us get selected adn sorted nodes indices
    node_indices = []
    for k in range (len(nodes_sorted)):
        node_indices.append (nodes_sorted [k].index)

    #let us make list of x coordinates of selected nodes
    x_node_tab = []
    for node in nodes_sorted:
        xyz = node.get_xyz()
        x_node_tab.append (xyz [0])

    #get uy profile from results
    #create nodal rsl instance
    nodal_rsl = NodalResults (mesh, his, rcf)

    for i, time in enumerate(time_instances):
        #get index of solution for t=time_ref
        solution_index = his.give_closest_sol_index_for_time (time)
        #get displ. profile for selected nodes
        #on the output we get a numpy array
        rsl = nodal_rsl.get_nodes_displacements (node_indices, solution_index, 'Y')
        #let us switch from [m] to [mm]
        rsl = rsl * 1000.0
        label = str(time)
        dg.add_plot (x_node_tab, rsl ,marker_=markers [i], markevery_=5, markersize_=6, \
                 markerfacecolor_=colors[i], markeredgecolor_='black', \
                 linestyle_='-',label_= label, color_=colors[i])

    dg.set_x_axes_tick_params (labelsize=16)
    dg.set_x_min   ( 30.0)
    dg.set_x_max   (120.0)
    dg.set_y_max   (  0.0)
    dg.set_y_min   (-20.0)
    #dg.set_x_log_scale ()
    dg.set_x_label (r'x [m]')
    dg.set_y_label (r'$s$ [mm]')
    dg.add_grid         ()
    dg.add_legend       ()
    dg.save_fig(out_png)


# =======================================================
def Wall_Mz (zsoil_data, time_instances, out_png):
# =======================================================
    # goal:
    # this function will plot bending moments in the wall at giuven time instances
    #
    # data_file      - name of the project (without file extension !)
    # time_instances - list of real time instances (float values) for which results are to be plot
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    # create instance of diagram class which will contain plot of settlement behinf the wall
    dg = Diagram()
    dg.set_figure_size((6,12))
    markers = dg.get_markers(len(time_instances))
    colors = dg.get_colors(len(time_instances))

    # ------------------------------------------------------------------
    # next 3 calls must always be done to be ablr to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution(data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh(data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info(data_file)
    # ------------------------------------------------------------------

    #let us select beam elements in the diaphragm wall
    beams = mesh.get_group_of_elements_ex (Element.GROUP_BEAM)

    #one could also select these elements by user material index
    # beams = mesh.get_list_of_elements (Element.GROUP_BEAM, time=0.0, only_active=False, \
    #                                    mat_filter=[3], ele_class_filter=[],
    #                                    zoom_filter=[0, [[], [], []]])

    # sort beam elements by Y coordinate
    point_top = [30.0,0.0 ,0.0]
    dir_vec   = [0.0,-1.0,0.0]
    beams_sorted = mesh.sort_sel_elements_by_dist_along_dir (beams,point_top,\
                                                                   dir_vec)

    # get beam element indices
    beam_ele_indices = []
    for k in range (len(beams_sorted)):
        beam_ele_indices.append (beams_sorted [k].index)

    # set array with sorted centers of beam elements
    y_tab = []
    for beam in beams_sorted:
        xyz = beam.get_center ()
        y_tab.append (xyz [1])

    #get Mz profile from results
    #create beam element results instance
    beams_rsl = Beam_EleResults (mesh, his, rcf)

    for i, time in enumerate(time_instances):
        #get index of solution for t=time_ref
        solution_index = his.give_closest_sol_index_for_time (time)
        #get Mz profile for selected beams
        rsl = beams_rsl.get_moment_for_sel_elements \
                    (beam_ele_indices, solution_index,'Z')
        rsl = rsl * (-1.0)
        label = str(time)
        dg.add_plot (rsl,y_tab,marker_=markers [i], markevery_=5, markersize_=8, \
                 markerfacecolor_=colors[i], markeredgecolor_='black', \
                 linestyle_='-',label_= label,color_=colors[i])

    dg.set_x_axes_tick_params (labelsize=16)
    dg.set_x_min   (-800.0)
    dg.set_x_max   ( 400.0)
    dg.set_y_max   (0.0)
    dg.set_y_min   (-32.0)
    #dg.set_x_log_scale ()
    dg.set_x_label (r'$M_z$ [kNm/m]')
    dg.set_y_label (r'$y$ [m]')
    dg.add_grid ()
    dg.add_legend()
    dg.save_fig(out_png)


# =======================================================
def Wall_Mz_envelope (zsoil_data, time_instances, outM_png, outN_png):
# =======================================================
    # goal:
    # this function will plot settlements behind the wall at selected time instances
    #
    # data_file      - name of the project (without file extension !)
    # time_instances - list of real time instances (float values) for which results are to be plot
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    # create instance of diagram class which will contain plot of settlement behinf the wall
    dg = Diagram()
    dg.set_figure_size((6,12))

    dg1 = Diagram()
    dg1.set_figure_size((6, 12))

    markers = dg.get_markers(len(time_instances))
    colors = dg.get_colors(len(time_instances))

    # ------------------------------------------------------------------
    # next 3 calls must always be done to be ablr to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution(data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh(data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info(data_file)
    # ------------------------------------------------------------------

    #let us select beam elements in the diaphragm wall
    beams = mesh.get_group_of_elements_ex (Element.GROUP_BEAM)

    #one could also select these elements by user material index
    # beams = mesh.get_list_of_elements (Element.GROUP_BEAM, time=0.0, only_active=False, \
    #                                    mat_filter=[3], ele_class_filter=[],
    #                                    zoom_filter=[0, [[], [], []]])

    # sort beam elements by Y coordinate
    point_top = [30.0,0.0 ,0.0]
    dir_vec   = [0.0,-1.0,0.0]
    beams_sorted = mesh.sort_sel_elements_by_dist_along_dir (beams,point_top,\
                                                                   dir_vec)

    # get beam element indices
    beam_ele_indices = []
    for k in range (len(beams_sorted)):
        beam_ele_indices.append (beams_sorted [k].index)

    # set array with sorted beam elements nodes top-down
    y_tab = []
    for beam in beams_sorted:
        nodes = beam.get_nodes ()
        xyz_1 = (nodes [0]).get_xyz()
        xyz_2 = (nodes [1]).get_xyz()
        y_tab.append (xyz_1 [1])
        y_tab.append (xyz_2 [1])

    #get Mz profile from results
    #create beam element results instance
    beams_rsl = Beam_EleResults (mesh, his, rcf)

    solution_indices = []
    for time in time_instances:
        solution_index = his.give_closest_sol_index_for_time (time)
        solution_indices.append (solution_index)

    #get envelope at nodes of beam elements
    # each envelope is an array [N_ele*2,6] # because we have a leading force + corresponding remaining ones
    Mz_minus, Mz_plus = \
        beams_rsl.get_moment_envelope_at_nodes_for_sel_elements (beam_ele_indices,\
                                                                solution_indices, 'Z')
    #index of leading force
    MZ = 5
    NX = 0

    #modify Mz points by y-sort order
    Mzp = []
    Nxp = []
    Mzm = []
    Nxm = []
    y   = []

    count = 0
    for k in range (len(y_tab)//2):
        if y_tab [count+0] < y_tab [count+1]:
            k1 = 1
            k2 = 0
        else:
            k1 = 0
            k2 = 1
        y.append (y_tab [count+k1])
        y.append (y_tab [count+k2])

        Mzm.append (Mz_minus[count+k1,MZ] * (-1.0))
        Mzm.append (Mz_minus[count+k2,MZ] * (-1.0))

        Nxm.append (Mz_minus[count+k1,NX])
        Nxm.append (Mz_minus[count+k2,NX])

        Mzp.append (Mz_plus [count+k1,MZ] * (-1.0))
        Mzp.append (Mz_plus [count+k2,MZ] * (-1.0))

        Nxp.append (Mz_plus [count+k1,NX])
        Nxp.append (Mz_plus [count+k2,NX])

        count = count + 2

    dg.add_plot (Mzm, y, marker_= None, markevery_=5, markersize_=8, \
                 markerfacecolor_= 'red', markeredgecolor_='black', \
                 linestyle_='-',label_= '',color_='red')


    dg1.add_plot(Nxm, y, marker_=None, markevery_=5, markersize_=8, \
                markerfacecolor_='red', markeredgecolor_='black', \
                linestyle_='-', label_='', color_='red')

    dg.add_plot (Mzp, y, marker_=None, markevery_=5, markersize_=8, \
                markerfacecolor_='blue', markeredgecolor_='black', \
                linestyle_='-', label_='', color_='blue')


    dg1.add_plot(Nxp, y, marker_=None, markevery_=5, markersize_=8, \
                markerfacecolor_='blue', markeredgecolor_='black', \
                linestyle_='-', label_='', color_='blue')

    dg.set_x_axes_tick_params (labelsize=16)
    dg.set_x_min   (-800.0)
    dg.set_x_max   ( 400.0)
    dg.set_y_max   (0.0)
    dg.set_y_min   (-32.0)
    #dg.set_x_log_scale ()
    dg.set_x_label (r'$M_z$ [kNm/m]')
    dg.set_y_label (r'$y$ [m]')
    dg.add_grid ()
    dg.add_legend()
    dg.save_fig(outM_png)


    dg1.set_x_axes_tick_params(labelsize=16)
    # dg1.set_x_min(-800.0)
    # dg1.set_x_max(0.0)
    dg1.set_y_max(0.0)
    dg1.set_y_min(-32.0)
    # dg.set_x_log_scale ()
    dg1.set_x_label(r'$N_x$ [kN/m]')
    dg1.set_y_label(r'$y$ [m]')
    dg1.add_grid()
    dg1.add_legend()
    dg1.save_fig(outN_png)


# ==========================================================
def Wall_cnt_stresses (zsoil_data, real_time, out_png):
# ==========================================================
    # goal:
    # this function will plot contact normal streses on both sides of the wall at selected time instances
    #
    # zsoil_data     - name of the project (without file extension !)
    # real_time      - real time for which results are to be plot
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    #create instance of diagram class which will contain plot contact stresses
    dg = Diagram ()
    dg.set_figure_size ((12,12))
    markers = dg.get_markers (len(time_instances))
    colors  = dg.get_colors  (len(time_instances))

    # ------------------------------------------------------------------
    #next 3 calls must always be done to be able to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution (data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh (data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info (data_file)
    # ------------------------------------------------------------------

    # here in this example we have only one material for the left and right contact hence we have to identify which is left and which is right

    # set graph in which each node keeps adjacent elements of given classes (here continuum and contact)
    mesh.set_adj_elements_to_nodes([Element.GROUP_BEAM, Element.GROUP_CONTACT])

    x = 30.0
    dx= 1.0e-3

    # select beam elements
    zoom_filter = [mesh.ZOOM_FILTER_ALL_NODES_MUST_BE_IN,[[x-dx,x+dx],[-32.001,0.001],[]]]
    beams = mesh.get_list_of_elements (Element.GROUP_BEAM, time=real_time, \
                                               only_active = True, zoom_filter=zoom_filter)

    #find contact elements adjacent to beams
    mat_cnt_filter = [7]
    adj_contact_elements = mesh.get_adj_elements_to (beams, Element.GROUP_CONTACT, at_time=time,target_ele_mat_filter=mat_cnt_filter)

    left_adj_cnt_elements  = []
    right_adj_cnt_elements = []
    #split one list to two lists using the fact that contact is oriented in opposite direction on the left and right side
    #the easiest way to handle this kind of the problem is to declare two different contact materials on the lft and right side of the wall
    for ele in adj_contact_elements:
        xyz = ele.get_ele_coord()
        dy  = xyz [1][1] - xyz [0][1]
        if dy < 0.0:
            left_adj_cnt_elements.append (ele)
        else:
            right_adj_cnt_elements.append (ele)

    # # as we do not know the order we need to sort these nodes
    point = [0.0,0.0 ,0.0]
    dir   = [0.0,1.0,0.0]

    sorted_sel_cnt_ele_L = mesh.sort_sel_elements_by_dist_along_dir(left_adj_cnt_elements, point, dir)
    indices_of_cnt_ele_L = mesh.get_ele_indices_for_sel_elements(sorted_sel_cnt_ele_L)

    sorted_sel_cnt_ele_R = mesh.sort_sel_elements_by_dist_along_dir(right_adj_cnt_elements, point, dir)
    indices_of_cnt_ele_R = mesh.get_ele_indices_for_sel_elements(sorted_sel_cnt_ele_R)

    cnt_rsl = Contact_2D_Results (mesh,his,rcf)
    solution_index = his.give_closest_sol_index_for_time (real_time)
    sig_nn_L = cnt_rsl.get_rsl_vec_for_sel_elements_ex(indices_of_cnt_ele_L, solution_index, \
                                                       Contact_2D_Results.stress_key, 'YY')
    sig_nn_L = sig_nn_L.reshape (sig_nn_L.size)

    sig_nn_R = cnt_rsl.get_rsl_vec_for_sel_elements_ex(indices_of_cnt_ele_R, solution_index, \
                                                Contact_2D_Results.stress_key, 'YY')
    sig_nn_R = sig_nn_R.reshape(sig_nn_R.size)

    # these results are stored at contact element nodes as nodal integration is used in the interfaces
    # hence jumps between elements are expected due to different penalty factors aet at each interface element
    # to avoid stress oscillations here we will compute averaged normal stresses at element centers
    # a singular situation appears at the bottom of the wall as two interface nodes are the same and therefore
    # artifficial zero stress is generated there

    y_L = []
    sig_nn_L_center = np.zeros (len(sorted_sel_cnt_ele_L))
    offs = 0
    for i, ele in enumerate(sorted_sel_cnt_ele_L):
        xyz = ele.get_center ()
        y_L.append (xyz [1])
        sig_nn = 0.5*(sig_nn_L [offs]+sig_nn_L [offs+1])
        sig_nn_L_center [i] = sig_nn
        offs = offs + ele.get_ngaus()

    y_R = []
    sig_nn_R_center = np.zeros (len(sorted_sel_cnt_ele_R))
    offs = 0
    for i, ele in enumerate(sorted_sel_cnt_ele_R):
        xyz = ele.get_center ()
        y_R.append (xyz [1])
        sig_nn = 0.5*(sig_nn_R [offs]+sig_nn_R [offs+1])
        sig_nn_R_center [i] = sig_nn
        offs = offs + ele.get_ngaus()

    label = 'left side'
    i = 0
    dg.add_plot (sig_nn_L_center, y_L,marker_=markers [i], markevery_=1, markersize_=3, \
                 markerfacecolor_=colors[i], markeredgecolor_='black', \
                 linestyle_='-',label_= label, color_=colors[i])

    i = 1
    label = 'right side'
    dg.add_plot(sig_nn_R_center*(-1.0), y_R, marker_=markers[i], markevery_=1, markersize_=3, \
            markerfacecolor_=colors[i], markeredgecolor_='black', \
            linestyle_='-', label_=label, color_=colors[i])

    # dg.set_x_axes_tick_params (labelsize=16)
    # dg.set_x_min   ( 30.0)
    # dg.set_x_max   (120.0)
    # dg.set_y_max   (  0.0)
    # dg.set_y_min   (-20.0)
    # #dg.set_x_log_scale ()
    # dg.set_x_label (r'x [m]')
    dg.set_y_label (r'$y$ [m]')
    dg.add_grid         ()
    dg.add_legend       ()
    dg.save_fig(out_png)


# ==========================================================
def Wall_eff_stresses_in_soil (zsoil_data, real_time, out_png):
# ==========================================================
    # goal:
    # this function will plot horizontal effective stresses on both sides of the wall at selected time instance
    #
    # zsoil_data     - name of the project (without file extension !)
    # real_time      - time at which the
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    #create instance of diagram class which will contain plot contact stresses
    dg = Diagram ()
    dg.set_figure_size ((12,12))
    markers = dg.get_markers (len(time_instances))
    colors  = dg.get_colors  (len(time_instances))

    # ------------------------------------------------------------------
    #next 3 calls must always be done to be able to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution (data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh (data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info (data_file)
    # ------------------------------------------------------------------

    x = 30.0
    dx= 1.0
    ddx= 0.001

    # select beam elements
    zoom_filter = [mesh.ZOOM_FILTER_ALL_NODES_MUST_BE_IN,[[x-dx-ddx,x+ddx],[-32.001,0.001],[]]]
    cont_L = mesh.get_list_of_elements (Element.GROUP_CONTINUUM, time=real_time, \
                                               only_active = True, zoom_filter=zoom_filter)


    zoom_filter = [mesh.ZOOM_FILTER_ALL_NODES_MUST_BE_IN, [[x -ddx, x + dx + ddx], [-32.001, 0.001], []]]
    cont_R = mesh.get_list_of_elements(Element.GROUP_CONTINUUM, time=real_time, \
                                   only_active=True, zoom_filter=zoom_filter)

    # # as we do not know the order we need to sort these nodes
    point = [0.0,0.0 ,0.0]
    dir   = [0.0,1.0,0.0]

    sorted_sel_ele_L = mesh.sort_sel_elements_by_dist_along_dir(cont_L, point, dir)
    indices_of_ele_L = mesh.get_ele_indices_for_sel_elements(sorted_sel_ele_L)

    sorted_sel_ele_R = mesh.sort_sel_elements_by_dist_along_dir(cont_R, point, dir)
    indices_of_ele_R = mesh.get_ele_indices_for_sel_elements(sorted_sel_ele_R)

    cont_rsl = Continuum_EleResults (mesh,his,rcf)
    solution_index = his.give_closest_sol_index_for_time (real_time)
    sig_xx_L = cont_rsl.get_rsl_for_sel_elements_ex(indices_of_ele_L, solution_index, \
                                                     Continuum_EleResults.STRESS_KEY, 'XX')
    sig_xx_L = sig_xx_L.reshape (sig_xx_L.size)


    sig_xx_R = cont_rsl.get_rsl_for_sel_elements_ex(indices_of_ele_R, solution_index, \
                                                     Continuum_EleResults.STRESS_KEY, 'XX')
    sig_xx_R = sig_xx_R.reshape(sig_xx_R.size)

    y_L = []
    for ele in sorted_sel_ele_L:
        xyz = ele.get_center ()
        y_L.append (xyz [1])

    y_R = []
    for ele in sorted_sel_ele_R:
        xyz = ele.get_center ()
        y_R.append (xyz [1])

    label = 'left side'
    i = 0
    dg.add_plot (sig_xx_L, y_L,marker_=markers [i], markevery_=1, markersize_=3, \
                 markerfacecolor_=colors[i], markeredgecolor_='black', \
                 linestyle_='-',label_= label, color_=colors[i])

    i = 1
    label = 'right side'
    dg.add_plot(sig_xx_R*(-1.0), y_R, marker_=markers[i], markevery_=1, markersize_=3, \
            markerfacecolor_=colors[i], markeredgecolor_='black', \
            linestyle_='-', label_=label, color_=colors[i])

    # dg.set_x_axes_tick_params (labelsize=16)
    # dg.set_x_min   ( 30.0)
    # dg.set_x_max   (120.0)
    # dg.set_y_max   (  0.0)
    # dg.set_y_min   (-20.0)
    # #dg.set_x_log_scale ()
    dg.set_x_label (r'$\sigma_{xx}$ [kPa]')
    dg.set_y_label (r'$y$ [m]')
    dg.add_grid         ()
    dg.add_legend       ()
    dg.save_fig(out_png)


if __name__ == '__main__':
    data_file = "..\\tests\\diaphragm-wall\\HS-Brick-Exc-Berlin-Sand-2phase"
    plots_dir = "..\\tutorials\\plots-tut-diaphragm-2D-diverse\\"
    his = HistoryOfExecution (data_file)
    list_of_solution_indices = his.give_converged_time_solutions()
    time_instances = his.get_data_value_for_sol_instances \
                     (list_of_solution_indices, his.DATA_TIME)
    time_instances = [3.0, 4.0, 6.0, 8.0, 9.0]
    Wall_defl  (data_file, time_instances, plots_dir+"Wall-Deflections.png")
    Wall_settl (data_file, time_instances, plots_dir+"Wall-settlements.png")
    Wall_Mz    (data_file, time_instances, plots_dir+"Wall-Mz.png")
    Wall_Mz_envelope (data_file, time_instances, plots_dir+"Wall-Mz-envelope.png", plots_dir+"Wall-Nx-corresponding.png")

    #plot contact stresses on both sides
    time = 9.0 #this is a real time analysis not a solution index
    Wall_cnt_stresses (data_file, time, plots_dir+"Wall-cnt-stresses.png")

    #plot horizntal effective stresses profiles on left/right sides of the wall
    time = 9.0 #this is a real time analysis not a solution index
    Wall_eff_stresses_in_soil (data_file, time, plots_dir+"Wall_eff_stresses_in_soil.png")

    pyplot.show ()
