from .C_Element import *
from .C_Node import *
from .C_Ltf import *
from .C_Exf import *
from .C_Material import *

from .C_Element_beam_L2 import *
from .C_Element_s_beam_L2 import *
from .C_Element_cont_Q4 import *
from .C_Element_cont_T3 import *
from .C_Element_cont_B8 import *
from .C_Element_cont_TH4 import *
from .C_Element_cont_W6 import *
from .C_Element_shell_Q4 import *
from .C_Element_shell_SHQ4 import *
from .C_Element_memb_L2 import *
from .C_Element_memb_T3 import *
from .C_Element_memb_Q4 import *
from .C_Element_truss_L2 import *
from .C_Element_cnt_L2 import *
from .C_Element_cnt_Q4 import *
from .C_Element_cnt_T3 import *
from .C_Element_cnt_bond_L2 import *
from .C_Element_cnt_tip import *
from .C_Element_heat_exch_L2 import *
from .C_Rcf_info import *
from .C_Rcf_layers_info import *
from math import *


# =====================================================
class Mesh():
    # =====================================================

    # ====================================================
    def __init__(self, project, analyze_rsl_files=True):
        # =====================================================

        self.my_project = project
        self.dict_ndm = {'PLANESTRAIN': 2, '3D': 3, 'AXISYMETRY': 2}
        self.ZOOM_FILTER_OFF = 0
        self.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN = 1
        self.ZOOM_FILTER_ALL_NODES_MUST_BE_IN = 2
        self.debug = False
        self.superconv_patch_recovery = {}

        ##        self.dict_elements_for_group = {'VOLUMICS' :[['Q4','Q4ES','T3'],['B8','B8ES','W6','TH4']],
        ##                                        'SHELLS'   :[[],['SXQ4','SHQ4','SXQ8','SXQ4']],
        ##                                        'CONTACT'  :[['C_L2','CB_L2'],['C_Q4','C_T3']],
        ##                                        'MEMBRANE' :[['M_L2'],['M_T3','M_Q4']],
        ##                                        'BEAMS'    :[['BEL2','SBEL2'],['BEL2','SBEL2']],
        ##                                        'TRUSSES'  :[['TRS2'],['TRS2']]
        ##                                        }
        self.groups_with_rsl = ['VOLUMICS', 'SHELLS', 'TRUSSES', 'BEAMS', 'CONTACT', 'MEMBRANE', 'HEAT_EXCH']
        f = open(project + ".dat", "rt")

        self.init_ndm(f, "")

        if self.debug:
            print("reading elements...")
        self.elements = self.init_elements(f, "")
        if self.debug:
            print("nr of read elements from dat file = ", len(self.elements))

        if self.debug:
            print("")
            print("reading materials...")
        self.materials = self.init_materials(f, "")
        if self.debug:
            print("nr of read materials from dat file = ", len(self.materials))

        if self.debug:
            print("")
            print("reading nodes...")
        self.nodes = self.init_nodes(f, "")
        if self.debug:
            print("nr of read nodes from dat file = ", len(self.nodes))
            print("")

        if self.debug:
            print("reading load time function...")
        self.ltfs = self.init_ltfs(f, "")
        if self.debug:
            print("nr of read load time functions from dat file = ", len(self.ltfs))
            print("")

        if self.debug:
            print("reading existence functions...")
        self.exfs = self.init_exfs(f, "")
        if self.debug:
            print("nr of read existence functions from dat file = ", len(self.exfs))
            print("")

        f.close()

        if analyze_rsl_files:
            self.set_seek_pos_in_rsl_files()
            self.set_seek_pos_in_lay_rsl_files()

    def create_superconv_patch_recovery(self, ele_group_index, at_time,
                                        **options):
        """Create/own an SPR manager for continuum or shells at physical time.

        The returned manager is also stored by group in
        self.superconv_patch_recovery. Call manager.update(time) to reuse it.
        """
        from .C_SuperConvPatchRecovery import SuperConvPatchRecovery
        recovery = SuperConvPatchRecovery(self, ele_group_index, at_time,
                                          **options)
        self.superconv_patch_recovery[ele_group_index] = recovery
        return recovery

    # =========================
    def init_ndm(self, f_dat, dat_filename):
        # =========================
        f = f_dat
        if f == None:
            f = open(dat_filename, "rt")

        line = ""
        done = False
        while not done:
            line = f.readline()
            if "JOB_TYPE" in line:
                done = True
        texts = line.split()
        job = (texts[1]).strip()

        self.ndm = self.dict_ndm[job]

    # =========================
    def init_elements(self, f_dat, dat_filename):
        # =========================
        f = f_dat
        if f == None:
            f = open(dat_filename, "rt")

        done = False
        elements = []

        line = ""
        while not done:
            pos = f.tell()
            line = f.readline()
            if "ELEM" in line:
                done = True
        texts = line.split()

        n = int(texts[1])

        for i in range(n):
            pos = f.tell()
            line = f.readline()
            texts = line.split()
            label = (texts[1]).strip()
            f.seek(pos)
            if label == "Q4" or label == "Q4ES":
                e = Element_cont_Q4(self)
            elif label == "T3":
                e = Element_cont_T3(self)
            elif label == "B8" or label == "B8ES":
                e = Element_cont_B8(self)
            elif label == "W6":
                e = Element_cont_W6(self)
            elif label == "TH4":
                e = Element_cont_TH4(self)
            elif label == "SXQ4":
                e = Element_shell_Q4(self)
            elif label == "SHQ4":
                e = Element_shell_SHQ4(self)
            elif label == "TRS2":
                e = Element_truss_L2(self)
            elif label == "BEL2":
                e = Element_beam_L2(self)
            elif label == "SBEL2":
                e = Element_s_beam_L2(self)
            elif label == "M_L2":
                e = Element_memb_L2(self)
            elif label == "M_Q4":
                e = Element_memb_Q4(self)
            elif label == "M_T3":
                e = Element_memb_Q4(self)
            elif label == "C_L2":
                e = Element_cnt_L2(self)
            elif label == "C_Q4":
                e = Element_cnt_Q4(self)
            elif label == "C_T3":
                e = Element_cnt_T3(self)
            elif label == "CB_L2":
                e = Element_cnt_bond_L2(self)
            elif label == "CNNPL":
                e = Element_cnt_tip(self)
            elif label == "HEXL2":
                e = Element_heat_exch_L2(self)
            elif label == "BHPPE":
                e = Element_heat_exch_L2(self)
            elif label == "BHE1U":
                e = Element_heat_exch_L2(self)
            elif label == "BHE2U":
                e = Element_heat_exch_L2(self)
            elif label == "BHCXA":
                e = Element_heat_exch_L2(self)
            elif label == "BHCXC":
                e = Element_heat_exch_L2(self)
            else:
                e = Element(self)

            e.instanciate(f)
            elements.append(e)

        if f_dat == None:
            f.close()

        return elements

    # =========================
    def init_nodes(self, f_dat, dat_filename):
        # =========================
        f = f_dat
        if f == None:
            f = open(dat_filename, "rt")

        done = False
        nodes = []

        line = ""
        while not done:
            line = f.readline()
            if "NODE" in line:
                done = True
        texts = line.split()

        n = int(texts[1])

        for i in range(n):
            node = Node(self)
            node.instanciate(f)
            nodes.append(node)

        if f_dat == None:
            f.close()

        return nodes

    # =========================
    def init_ltfs(self, f_dat, dat_filename):
        # =========================
        f = f_dat
        if f == None:
            f = open(dat_filename, "rt")

        done = False
        ltfs = []

        line = ""
        while not done:
            line = f.readline()
            if "LOADTIME" in line:
                done = True
        texts = line.split()

        n = int(texts[1])

        for i in range(n):
            ltf = Ltf(self)
            ltf.instanciate(f)
            ltfs.append(ltf)

        if f_dat == None:
            f.close()

        return ltfs

    # =========================
    def init_exfs(self, f_dat, dat_filename):
        # =========================
        f = f_dat
        if f == None:
            f = open(dat_filename, "rt")

        done = False
        exfs = []

        line = ""
        while not done:
            line = f.readline()
            if "EXISTFUN" in line:
                done = True
        texts = line.split()

        n = int(texts[1])

        for i in range(n):
            exf = Exf(self)
            exf.instanciate(f)
            exfs.append(exf)

        if f_dat == None:
            f.close()

        return exfs

    # =========================
    def init_materials(self, f_dat, dat_filename):
        # =========================
        f = f_dat
        if f == None:
            f = open(dat_filename, "rt")

        done = False
        materials = []

        line = ""
        while not done:
            line = f.readline()
            if "PROP" in line:
                done = True
        texts = line.split()

        n = int(texts[1])

        for i in range(n):
            mat = Material(self)
            mat.instanciate(f)
            materials.append(mat)

        if f_dat == None:
            f.close()

        return materials

    # =========================
    def set_seek_pos_in_rsl_files(self):
        # =========================
        rcf = RCF_info(self.my_project)
        for group in self.groups_with_rsl:
            last_pos = 0
            aux = self.get_group_of_elements(group)
            gp_size = rcf.give_one_gp_size(group)
            for e in aux:
                e.seek_in_rsl = last_pos
                last_pos = last_pos + len(e.xsiGP) * gp_size * 4  # in bytes

    # =========================
    def set_seek_pos_in_lay_rsl_files(self):
        # =========================
        rcf_lay = RCF_layers_info(self.my_project)
        for group in self.groups_with_rsl:
            last_pos = 0
            aux = self.get_group_of_elements(group)
            # gp_size= rcf.give_one_gp_size (group)
            lay_size = rcf_lay.give_one_lay_size(group)
            for e in aux:
                e.seek_in_lay_rsl = last_pos
                nlayers_all, nlayers_reinf = e.get_nlayers()
                ngaus = len(e.xsiGP)
                last_pos = last_pos + ngaus * nlayers_all * lay_size * 4  # in bytes

    # =========================
    def give_cell_with_string(self, cells, string):
        # =========================
        for i in range(len(cells)):
            a = cells[i]
            a = a.Trim()
            if a == string:
                return i
        return -1

    # =========================
    def give_Node_at_XYZ(self, X):
        # =========================
        # this function returs node index (starting from 1)
        dist_min = 1.0e38
        node_min = -1

        for i in range(len(self.nodes)):
            x_node = self.nodes[i].xyz
            dist = 0.0
            for k in range(min(len(x_node), len(X))):
                dist = dist + (X[k] - x_node[k]) * (X[k] - x_node[k])
            if dist < dist_min:
                dist_min = dist
                node_min = i + 1
        return node_min

    # =====================================================
    def get_list_of_elements(self, group_filter_in, time=0.0, only_active=False, mat_filter=[], ele_class_filter=[],
                             zoom_filter=[0, [[], [], []]]):
        # =====================================================
        if type(group_filter_in) == type(""):
            group_filter = group_filter_in
        elif type(group_filter_in) == type(1):
            group_filter = Element.group_dict[group_filter_in]
        else:
            print("undefined group filter; use Element.GROUP_CONTINUUM or other group index (class C_Element.py)")
        selection = []

        for e in self.elements:
            # print e.group, group

            material = self.materials[e.material_index - 1]
            mat_filter_pass = True
            if len(mat_filter) > 0:
                if not material.index_in_inp in mat_filter:
                    mat_filter_pass = False

            group_ele_filter_pass = e.group == group_filter

            exf_filter_pass = True

            if only_active:
                exf = e.give_exf()
                if exf != None:
                    if not exf.is_ON(time):
                        exf_filter_pass = False

            ele_class_filter_pass = True
            if len(ele_class_filter) > 0:
                if not e.label in ele_class_filter:
                    ele_class_filter_pass = False

            pass_test = group_ele_filter_pass and mat_filter_pass and exf_filter_pass and ele_class_filter_pass
            if pass_test:
                selection.append(e)

        if zoom_filter[0] == self.ZOOM_FILTER_OFF:
            return selection

        aux_nodes = []
        aux_elements = []

        for i in range(len(self.nodes)):
            aux_nodes.append(0)

        for e in selection:

            for i in range(e.nen):
                node_index = e.nodes[i] - 1
                aux_nodes[node_index] = 1

        for i in range(len(aux_nodes)):
            if aux_nodes[i] == 1:
                node = self.nodes[i]
                zooms = zoom_filter[1]
                for j in range(len(zooms)):
                    zoom = zooms[j]
                    if len(zoom) == 2:
                        xj_min = min(zoom[0], zoom[1])
                        xj_max = max(zoom[0], zoom[1])
                        if node.xyz[j] < xj_min or node.xyz[j] > xj_max:
                            aux_nodes[i] = 0

        for e in selection:
            count = 0
            for i in range(e.nen):
                node_index = e.nodes[i] - 1
                count = count + aux_nodes[node_index]
            if zoom_filter[0] == self.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN:
                if count > 0:
                    aux_elements.append(e)
            elif zoom_filter[0] == self.ZOOM_FILTER_ALL_NODES_MUST_BE_IN:
                if count == e.nen:
                    aux_elements.append(e)

        return aux_elements

    # =====================================================
    def get_group_of_elements(self, group_string):
        # =====================================================
        aux = self.get_list_of_elements(group_string)
        return aux

    # =====================================================
    def get_group_of_elements_ex(self, group_index):
        # =====================================================
        aux = self.get_list_of_elements(Element.group_dict[group_index])
        return aux

    # =====================================================
    def get_sum_of_gauss_points(self, group):
        # =====================================================
        count = 0
        for e in self.elements:
            if e.group == group:
                count = count + len(e.xsiGP)
        return count

    # =====================================================
    def get_sum_of_layers(self, group):
        # =====================================================
        count = 0
        for e in self.elements:
            if e.group == group:
                nlayers_all, nlayers_reinf = e.get_nlayers()
                count = count + len(e.xsiGP) * nlayers_all
        return count

    # =====================================================
    def get_list_of_nodes(self, ele_filter=[], zoom_filter=[0, [[], [], []]]):
        # =====================================================
        selection = []
        node_flags = []
        for i in range(len(self.nodes)):
            node_flags.append(1)

        if len(ele_filter) > 0:
            for e in ele_filter:
                for i in range(len(e.nodes)):
                    node_flags[e.nodes[i] - 1] = -1
            for i in range(len(self.nodes)):
                if node_flags[i] != -1:
                    node_flags[i] = 0
                else:
                    node_flags[i] = 1

        out = []

        for i in range(len(node_flags)):
            if node_flags[i] == 1:
                node = self.nodes[i]
                zooms = zoom_filter[1]
                for j in range(len(zooms)):
                    zoom = zooms[j]
                    if len(zoom) == 2:
                        xj_min = min(zoom[0], zoom[1])
                        xj_max = max(zoom[0], zoom[1])
                        if node.xyz[j] < xj_min or node.xyz[j] > xj_max:
                            node_flags[i] = 0
                if node_flags[i] == 1:
                    out.append(self.nodes[i])
        return out

    # =====================================================
    def get_connected_elements(self):
        # =====================================================
        out = []
        for i in range(len(self.nodes)):
            out.append(0)

        for i in range(len(self.elements)):
            e = self.elements[i]
            for inode in range(e.nen):
                node = e.nodes[inode]
                out[node - 1] = out[node - 1] + 1

        return out

    # =====================================================
    def get_node(self, node_index_1):
        # =====================================================
        return self.nodes[node_index_1 - 1]

    # =====================================================
    def get_element(self, ele_index_1):
        # =====================================================
        return self.elements[ele_index_1 - 1]

    # =====================================================
    def sort_sel_elements_by_dist_along_dir(self, sel_elements, point, dir_vec, \
                                            ret_proj_measures=False):
        # =====================================================
        pt = [0.0, 0.0, 0.0]
        dr = [0.0, 0.0, 0.0]
        pt = point[:]
        dr = dir_vec[:]
        # normalize direction vector
        dr_length = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
        dr_length = sqrt(dr_length)
        for i in range(len(dr)):
            dr[i] = dr[i] / dr_length

        proj_measures = []
        for ele in sel_elements:
            xyz = ele.get_center()
            dx = [xyz[0] - pt[0], xyz[1] - pt[1], xyz[2] - pt[2]]
            aux = dx[0] * dir_vec[0] + dx[1] * dir_vec[1] + dx[2] * dir_vec[2]
            proj_measures.append(aux)
        # now sort projection measures
        sort_indices = sorted(list(range(len(proj_measures))), key=lambda x: proj_measures[x])

        sel_ele_sorted = []
        proj_measures_sorted = []
        for i in sort_indices:
            sel_ele_sorted.append(sel_elements[i])
            proj_measures_sorted.append(proj_measures[i])

        if ret_proj_measures:
            return sel_ele_sorted, proj_measures_sorted
        else:
            return sel_ele_sorted

    # =====================================================
    def sort_sel_nodes_by_dist_along_dir(self, sel_nodes, point, dir_vec, \
                                         ret_proj_measures=False):
        # =====================================================
        pt = [0.0, 0.0, 0.0]
        dr = [0.0, 0.0, 0.0]
        pt = point[:]
        dr = dir_vec[:]
        # normalize direction vector
        dr_length = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
        dr_length = sqrt(dr_length)
        for i in range(len(dr)):
            dr[i] = dr[i] / dr_length

        proj_measures = []
        for node in sel_nodes:
            xyz = node.get_xyz()
            dx = [xyz[0] - pt[0], xyz[1] - pt[1], xyz[2] - pt[2]]
            aux = dx[0] * dir_vec[0] + dx[1] * dir_vec[1] + dx[2] * dir_vec[2]
            proj_measures.append(aux)
        # now sort projection measures
        sort_indices = sorted(list(range(len(proj_measures))), key=lambda x: proj_measures[x])

        sel_nodes_sorted = []
        proj_measures_sorted = []        
        for i in sort_indices:
            sel_nodes_sorted.append(sel_nodes[i])
            proj_measures_sorted.append(proj_measures[i])            

        if ret_proj_measures:
            return sel_nodes_sorted, proj_measures_sorted
        else:
            return sel_nodes_sorted

    # =====================================================
    def get_ele_indices_for_sel_elements(self, sel_elements):
        # =====================================================
        ele_indices = []
        for ele in sel_elements:
            ele_indices.append(ele.index)
        return ele_indices

    # =====================================================
    def get_node_indices_for_sel_nodes(self, sel_nodes):
        # =====================================================
        node_indices = []
        for node in sel_nodes:
            node_indices.append(node.index)
        return node_indices

    # =====================================================
    def find_material_by_label(self, label):
        # =====================================================
        for material in self.materials:
            if material.label == label:
                return material
        return None

    # =====================================================
    def find_materials_by_labels(self, labels):
        # =====================================================
        out = []
        for label in labels:
            material = self.find_material_by_label(label)
            if material != None:
                out.append(material)
        return out

    # =====================================================
    def find_first_element_with_mat_index(self, mat_index, use_index_in_inp_flag=True):
        # =====================================================
        for element in self.elements:
            material = self.materials[element.get_material_index()]
            if use_index_in_inp_flag:
                if material.index_in_inp == mat_index:
                    return element
            else:
                if material.index_in_dat == mat_index:
                    return element
        return None

    # =====================================================
    def set_adj_elements_to_nodes(self, ele_groups=[], mat_filter = []):
        # =====================================================
        for node in self.nodes:
            node.clear_adj_elements()

        for element in self.elements:
            attach = True
            if len(ele_groups) > 0:
                if element.get_group_index() not in ele_groups:
                    attach = False
            if attach and len(mat_filter) > 0:
                material = self.materials [element.material_index-1]
                if not material.index_in_inp in mat_filter:
                    attach = False
            if attach:
                for node_index in element.nodes:
                    if node_index != 0:
                        node = self.get_node(node_index)
                        node.attach_element_to_adj_list(element)

    # =====================================================
    def get_adj_elements_to_sel_nodes(self, sel_nodes, ele_groups=[]):
        # =====================================================
        # make copy of list of nodes
        nodes_aux = []
        for node in self.nodes:
            nodes_aux.append(node.deepcopy())

        for node in nodes_aux:
            nodes_aux.clear_adj_elements()

        flags = []
        for i in range(len(nodes_aux)):
            flags.append(0)

        for node in sel_nodes:
            flags[node.index - 1] = 1

        for element in self.elements:
            attach = True
            if len(ele_groups) > 0:
                if element.get_group_index() not in ele_groups:
                    attach = False
            if attach:
                for node_index in element.nodes:
                    if node_index != 0:
                        if flags[node_index - 1] == 1:
                            node = nodes_aux[node_index - 1]
                            node.attach_element_to_adj_list(element)

    # # =====================================================
    # def get_adj_elements_to(self, src_list, target_ele_group_index, at_time=None):
    #     # =====================================================
    #     out = []
    #     for element in src_list:
    #         out.append(None)
    #         aux = []
    #         for node_index in element.nodes:
    #             node = self.get_node(node_index)
    #             adj_ele_to_node = node.get_adj_elements(target_ele_group_index, at_time)
    #             tmp = []
    #             for adj_ele in adj_ele_to_node:
    #                 tmp.append(adj_ele.index)
    #             if tmp != []:
    #                 aux.append(tmp)
    #         common = set(aux[0])
    #         for i in range(len(aux) - 1):
    #             next = aux[i + 1]
    #             tmp = list(set(common) & set(next))
    #             common = tmp[:]
    #         if len(common) == 1:
    #             out[-1] = self.get_element(common[0])
    #     return out

    # =====================================================
    def get_adj_elements_to(self, src_list, target_ele_group_index, at_time=None, target_ele_mat_filter=[]):
        # =====================================================
        out = []
        for element in src_list:
            aux = []
            for node_index in element.nodes:
                node = self.get_node(node_index)
                adj_ele_to_node = node.get_adj_elements(target_ele_group_index, at_time)
                tmp = []
                for adj_ele in adj_ele_to_node:
                    if target_ele_mat_filter == []:
                        tmp.append(adj_ele.index)
                    else:
                        material = self.materials [adj_ele.get_material_index()]
                        if material.get_user_index() in target_ele_mat_filter:
                            tmp.append(adj_ele.index)
                if tmp != []:
                    aux.append(tmp)
            if aux != []:
                common = set(aux[0])
                for i in range(len(aux) - 1):
                    next = aux[i + 1]
                    tmp = list(set(common) & set(next))
                    common = tmp[:]
                for i in range (len(common)):
                #if len(common) == 1:
                    #out[-1] = self.get_element(common[0])
                    out.append (self.get_element(common[i]))
        return out

    # =====================================================
    def get_ltf(self, ltf_index_1):
        # =====================================================
        if ltf_index_1 <= 0:
            return None
        if ltf_index_1 > len(self.ltfs):
            return None
        return self.ltfs[ltf_index_1 - 1]

    # =====================================================
    def get_continuum_ext_faces(self, only_active=False, time=0.0):
        # =====================================================
        # return list of [ [ele_index_1, face_index_1],.....] being on the external contour
        cont_ele = self.get_list_of_elements(Element.GROUP_CONTINUUM, only_active=only_active, time=time)
        ext_faces = []
        for ele in cont_ele:
            skip = False
            if only_active:
                if not ele.is_ON(time):
                    skip = True
            if not skip:
                # print ele.index
                ref_ele = ele.get_ref_ele()
                nfaces = ref_ele.get_nr_of_faces()
                ele_nodes = ele.get_nodes()
                for i in range(nfaces):
                    loc_face_nodes = ref_ele.get_face_node_indices(i + 1)
                    set_aux = []
                    for j in range(len(loc_face_nodes)):
                        node = ele_nodes[loc_face_nodes[j] - 1]
                        adj_ele = node.get_adj_elements(Element.GROUP_CONTINUUM, time)
                        adj_ele_indices = []
                        for adj_ele_aux in adj_ele:
                            adj_ele_indices.append(adj_ele_aux.index)
                        if j == 0:
                            set_aux[:] = adj_ele_indices[:]
                        else:
                            # common set
                            c = list(set(set_aux) & set(adj_ele_indices))
                            set_aux = []
                            set_aux[:] = c[:]
                    if len(set_aux) == 1:
                        ext_faces.append([ele.index, i + 1])

        # ext_faces_dict = {}
        # for ele in cont_ele:
        #     print ele.index
        #     ref_ele = ele.get_ref_ele()
        #     nfaces = ref_ele.get_nr_of_faces()
        #     ele_nodes = ele.get_nodes()
        #     for i in range (nfaces):
        #         loc_face_nodes = ref_ele.get_face_node_indices(i+1)
        #         glo_face_nodes = []
        #         for j in range (len(loc_face_nodes)):
        #             node = ele_nodes [loc_face_nodes [j]-1]
        #             glo_face_nodes.append (node.index)
        #         glo_face_nodes.sort()
        #         #print glo_face_nodes
        #         key = ''
        #         for item in glo_face_nodes:
        #             key = key + str(item)
        #         if key != '':
        #             if key in ext_faces_dict.keys():
        #                 data  = ext_faces_dict [key]
        #                 data [2] = data [2] + 1
        #                 ext_faces_dict[key] = data
        #             else:
        #                 ext_faces_dict [key] = [ele.index, i+1, 1]
        #
        # ext_faces = []
        # for key in ext_faces_dict.iterkeys ():
        #     data = ext_faces_dict [key]
        #     if data [2] == 1:
        #         ext_faces.append ( [data [0], data [1]])
        return ext_faces


# =====================================================
def main():
    # =====================================================
    mesh = Mesh('..\ExamplesData\Plate-1')
    mesh.set_adj_elements_to_nodes([Element.GROUP_SHELL])
    node = mesh.nodes[0]
    x = node.get_xyz()
    # aux = mesh.get_list_of_elements ("VOLUMICS",5.0,True,[1,2],["Q4","Q4ES"])
    # take beams only with material 5, active at time = 4.0 of any element class
    # aux = mesh.get_list_of_elements ("BEAMS",4.0,True,[5],[])
    ##aux = mesh.get_list_of_elements ("BEAMS",4.0,True,[5],[],[mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN,[[14.29,14.31],[-7.499,-0.001],[]]])
    ##print len(aux)


##    for e in mesh.elements:
##        print e.index,e.label,e.nodes
##    for node in mesh.nodes:
##        print node.index,node.xyz


if __name__ == '__main__':
    main()
