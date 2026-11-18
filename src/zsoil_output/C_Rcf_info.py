"""Read result-configuration metadata for nodal and element result files.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`RCF_info` reader maps result names and components to binary-file columns,
reports record sizes, and compares or maps configurations from different ZSoil projects.
"""


class RCF_info ():
    def __init__(self, project):
        """Initialize the result configuration info object for the supplied project data.

        :param project: Path or filename prefix of the ZSoil project.
        :type project: str
        """
        self.my_project = project
        f = open(project+".rcf", "rt")
        lines = f.readlines()

        n_ele_groups = int(lines[0])
        count = 1

        self.ele_rcf = []

        for i in range(n_ele_groups):

            a = lines[count]
            a = a.split()

            group_id = a[0]
            nitems = int(a[1])

            aux1 = []

            for j in range(nitems):
                count = count + 1
                a = lines[count]
                a = a.split()

                rsl_item = a[0]

                # correcting misprints in older versions
                if rsl_item == "STRESESS":
                    rsl_item = "STRESSES"

                n_comp = int(a[1])
                rsl_items = []

                #print 'item = ',j,rsl_item,n_comp

                if n_comp > 1:
                    for i in range(n_comp):
                        if i+2 >= len(a):
                            rsl_items.append(" ")
                        else:
                            rsl_items.append(a[i+2])
                else:
                    if len(a) > 2:
                        rsl_items.append(a[2])
                    else:
                        rsl_items.append('')

                aux2 = [rsl_item, n_comp, rsl_items]
                aux1.append(aux2)
            count = count + 1

            aux3 = [group_id, aux1]
            self.ele_rcf.append(aux3)

        # nodes
        self.nod_rcf = []
        a = lines[count]
        a = a.split()

        nitems = int(a[0])

        for j in range(nitems):
            count = count + 1
            a = lines[count]
            a = a.split()

            rsl_item = a[0]
            n_comp = int(a[1])
            rsl_items = []

            #print 'item = ',j,rsl_item,n_comp

            # AT 15-07-2020
            # if n_comp > 1:
            if (len(a) > 2):
                for i in range(n_comp):
                    rsl_items.append(a[i+2])
            else:
                rsl_items.append('')

            aux2 = [rsl_item, n_comp, rsl_items]
            self.nod_rcf.append(aux2)

        # nodal velocities
        count = count + 1
        self.nod_rcf_v = []
        a = lines[count]
        a = a.split()

        nitems = int(a[0])

        for j in range(nitems):
            count = count + 1
            a = lines[count]
            a = a.split()

            rsl_item = a[0]
            n_comp = int(a[1])
            rsl_items = []

            #print 'item = ',j,rsl_item,n_comp

            if n_comp > 1:
                for i in range(n_comp):
                    rsl_items.append(a[i+2])
            else:
                rsl_items.append('')

            aux2 = [rsl_item, n_comp, rsl_items]
            self.nod_rcf_v.append(aux2)

        count = count + 1
        # nodal accelerations
        self.nod_rcf_a = []
        a = lines[count]
        a = a.split()

        nitems = int(a[0])

        for j in range(nitems):
            count = count + 1
            a = lines[count]
            a = a.split()

            rsl_item = a[0]
            n_comp = int(a[1])
            rsl_items = []

            #print 'item = ',j,rsl_item,n_comp

            if n_comp > 1:
                for i in range(n_comp):
                    rsl_items.append(a[i+2])
            else:
                rsl_items.append('')

            aux2 = [rsl_item, n_comp, rsl_items]
            self.nod_rcf_a.append(aux2)

        count = count + 1
        # nodal modes
        self.nod_rcf_modes = []
        a = lines[count]
        a = a.split()

        nitems = int(a[0])

        for j in range(nitems):
            count = count + 1
            a = lines[count]
            a = a.split()

            rsl_item = a[0]
            n_comp = int(a[1])
            rsl_items = []

            #print 'item = ',j,rsl_item,n_comp

            if n_comp > 1:
                for i in range(n_comp):
                    rsl_items.append(a[i+2])
            else:
                rsl_items.append('')

            aux2 = [rsl_item, n_comp, rsl_items]
            self.nod_rcf_modes.append(aux2)

        f.close()

    def give_ele_rcf_items_for_group(self, group_id_):
        """Return element result configuration items for group.

        :param group_id_: Group id used by this operation.
        :type group_id_: int
        :return: The element result configuration items for group.
        :rtype: list
        """
        for i in range(len(self.ele_rcf)):
            group_rcf = self.ele_rcf[i]
            group_id = (group_rcf[0]).strip()
            if group_id == group_id_.strip():
                list_items = group_rcf[1]
                return list_items
        return []

    def give_col_for_ele_rsl(self, group_id_, rsl_item_, comp_id_):
        """Return column for element results.

        :param group_id_: Group id used by this operation.
        :type group_id_: int
        :param rsl_item_: Result item used by this operation.
        :type rsl_item_: str
        :param comp_id_: Component id used by this operation.
        :type comp_id_: str
        :return: The column for element results.
        :rtype: tuple[int, int]
        """
        for i in range(len(self.ele_rcf)):
            group_rcf = self.ele_rcf[i]
            group_id = group_rcf[0]
            if group_id.find(group_id_) != -1:
                list_items = group_rcf[1]
                col = 0
                for j in range(len(list_items)):
                    item = list_items[j]
                    if item[0].find(rsl_item_) != -1:
                        n_comps = item[1]
                        comps = item[2]
                        for k in range(n_comps):
                            if comps[k].find(comp_id_) != -1:
                                return col, n_comps
                            else:
                                col = col + 1
                    else:
                        n_comps = item[1]
                        col = col + n_comps
        return -1, 0

    def give_col_for_nod_rsl(self, rsl_item_, comp_id_):
        """Return the storage column and component count for a nodal result item.

        :param rsl_item_: Result item used by this operation.
        :type rsl_item_: str
        :param comp_id_: Component id used by this operation.
        :type comp_id_: str
        :return: The storage column and component count for a nodal result item.
        :rtype: tuple[int, int]
        """
        col = 0
        for j in range(len(self.nod_rcf)):
            item = self.nod_rcf[j]
##                    print item [0]
##                    print item [1]
##                    print item [2]
            if item[0].find(rsl_item_) != -1:
                n_comps = item[1]
                comps = item[2]
                for k in range(n_comps):
                    if comps[k].find(comp_id_) != -1:
                        return col, n_comps
                    else:
                        col = col + 1
            else:
                n_comps = item[1]
                col = col + n_comps
        return -1, 0

    def give_one_gp_size(self, group_id_):
        """Return one integration point size.

        :param group_id_: Group id used by this operation.
        :type group_id_: int
        :return: One integration point size.
        :rtype: int
        """
        for i in range(len(self.ele_rcf)):
            group_rcf = self.ele_rcf[i]
            group_id = (group_rcf[0]).strip()
            if group_id == group_id_.strip():
                list_items = group_rcf[1]
                count = 0
                for j in range(len(list_items)):
                    item = list_items[j]
                    n_comps = item[1]
                    count = count + n_comps
                return count
        return -1

    def give_one_node_size(self):
        """Return the number of stored result values per node.

        :return: The number of stored result values per node.
        :rtype: int
        """
        count = 0
        for j in range(len(self.nod_rcf)):
            item = self.nod_rcf[j]
            n_comps = item[1]
            count = count + n_comps
        return count

    def give_one_node_size_v(self):
        """Return the number of stored velocity values per node.

        :return: The number of stored velocity values per node.
        :rtype: int
        """
        count = 0
        for j in range(len(self.nod_rcf_v)):
            item = self.nod_rcf_v[j]
            n_comps = item[1]
            count = count + n_comps
        return count

    def give_one_node_size_a(self):
        """Return the number of stored acceleration values per node.

        :return: The number of stored acceleration values per node.
        :rtype: int
        """
        count = 0
        for j in range(len(self.nod_rcf_a)):
            item = self.nod_rcf_a[j]
            n_comps = item[1]
            count = count + n_comps
        return count

    def give_one_node_size_modes(self):
        """Return the number of stored modal-result values per node.

        :return: The number of stored modal-result values per node.
        :rtype: int
        """
        count = 0
        for j in range(len(self.nod_rcf_modes)):
            item = self.nod_rcf_modes[j]
            n_comps = item[1]
            count = count + n_comps
        return count

    def is_equivalent_to(self, rcf_):
        """Return whether this object stores the same configuration or records as another instance.

        :param rcf_: Result configuration to compare or map against.
        :type rcf_: RCF_info
        :return: Whether this object stores the same configuration or records as another instance.
        :rtype: bool
        """
        # check nodal results
        if self.nod_rcf != rcf_.nod_rcf:
            return False
        if self.nod_rcf_v != rcf_.nod_rcf_v:
            return False
        if self.nod_rcf_a != rcf_.nod_rcf_a:
            return False
        if self.nod_rcf_modes != rcf_.nod_rcf_modes:
            return False
        if self.ele_rcf != rcf_.ele_rcf:
            return False

        return True

    def get_items_map(self, rcf_group_this, rcf_group_, excluded_items_):
        """Map each result component in this configuration to its column in another configuration.

        :param rcf_group_this: Rcf group this used by this operation.
        :type rcf_group_this: Sequence
        :param rcf_group_: Rcf group used by this operation.
        :type rcf_group_: Sequence
        :param excluded_items_: Excluded items used by this operation.
        :type excluded_items_: Sequence[str]
        :return: Column index in the other configuration for each result component, or -1.
        :rtype: list
        """
        items_list_this = []
        items_list_all_this = []

        for j in range(len(rcf_group_this)):
            item = rcf_group_this[j]
            n_comps = item[1]
            comps = item[2]
            if not item[0].strip() in excluded_items_:
                for comp in comps:
                    items_list_this.append(item[0].strip() + "->"+comp.strip())
            for comp in comps:
                items_list_all_this.append(item[0].strip() + "->"+comp.strip())

        items_list_other = []
        items_list_all_other = []

        for j in range(len(rcf_group_)):
            item = rcf_group_[j]
            n_comps = item[1]
            comps = item[2]
            if not item[0].strip() in excluded_items_:
                for comp in comps:
                    items_list_other.append(item[0].strip() + "->" + comp.strip())
            for comp in comps:
                items_list_all_other.append(item[0].strip() + "->" + comp.strip())

        # common part of two sets of rcf items and their components
        common_set = (set(items_list_this)).intersection(set(items_list_other))
        common_list = list(common_set)

        cols = []
        for i in range(len(items_list_all_this)):
            cols.append(-1)
        for i in range(len(items_list_all_this)):
            if items_list_all_this[i] in common_list:
                cols[i] = items_list_all_other.index(items_list_all_this[i])
        return cols

    def get_nod_items_map(self, rcf_, excluded_items_=[]):
        """Return the column mapping for standard nodal results in another configuration.

        :param rcf_: Result configuration to compare or map against.
        :type rcf_: RCF_info
        :param excluded_items_: Excluded items used by this operation.
        :type excluded_items_: Sequence[str]
        :return: The column mapping for standard nodal results in another configuration.
        :rtype: list
        """
        return self.get_items_map(self.nod_rcf, rcf_.nod_rcf, excluded_items_)

    def get_nod_v_items_map(self, rcf_, excluded_items_=[]):
        """Return the column mapping for nodal velocity results in another configuration.

        :param rcf_: Result configuration to compare or map against.
        :type rcf_: RCF_info
        :param excluded_items_: Excluded items used by this operation.
        :type excluded_items_: Sequence[str]
        :return: The column mapping for nodal velocity results in another configuration.
        :rtype: list
        """
        return self.get_items_map(self.nod_rcf_v, rcf_.nod_rcf_v, excluded_items_)

    def get_nod_a_items_map(self, rcf_, excluded_items_=[]):
        """Return the column mapping for nodal acceleration results in another configuration.

        :param rcf_: Result configuration to compare or map against.
        :type rcf_: RCF_info
        :param excluded_items_: Excluded items used by this operation.
        :type excluded_items_: Sequence[str]
        :return: The column mapping for nodal acceleration results in another configuration.
        :rtype: list
        """
        return self.get_items_map(self.nod_rcf_a, rcf_.nod_rcf_a, excluded_items_)

    def get_nod_modes_items_map(self, rcf_, excluded_items_=[]):
        """Return the column mapping for nodal modal results in another configuration.

        :param rcf_: Result configuration to compare or map against.
        :type rcf_: RCF_info
        :param excluded_items_: Excluded items used by this operation.
        :type excluded_items_: Sequence[str]
        :return: The column mapping for nodal modal results in another configuration.
        :rtype: list
        """
        return self.get_items_map(self.nod_rcf_modes, rcf_.nod_rcf_modes, excluded_items_)

    def get_ele_items_map(self, rcf_, ele_group_id_, excluded_items_=[]):
        """Return the column mapping for one element group in another result configuration.

        :param rcf_: Result configuration to compare or map against.
        :type rcf_: RCF_info
        :param ele_group_id_: Element group id used by this operation.
        :type ele_group_id_: int
        :param excluded_items_: Excluded items used by this operation.
        :type excluded_items_: Sequence[str]
        :return: The column mapping for one element group in another result configuration.
        :rtype: list
        """
        rcf_items_this = self.give_ele_rcf_items_for_group(ele_group_id_)
        rcf_items_other = rcf_.give_ele_rcf_items_for_group(ele_group_id_)
        return self.get_items_map(rcf_items_this, rcf_items_other, excluded_items_)


def main():
    """Run this module's example or command-line workflow."""
    rcf = RCF_info('D:\\v20\\inp\\testy\\MC-EXC-LONDON-CLAY-2PHASE-E6000Z')
    rcf1 = RCF_info('D:\\v20\\inp\\testy\\MC-EXC-LONDON-CLAY-2PHASE-E3600Z')
    print(rcf.is_equivalent_to(rcf1))
    cols = rcf.get_nod_items_map(rcf1, excluded_items_=[])
    print(cols)
    cols = rcf1.get_nod_items_map(rcf, excluded_items_=[])
    print(cols)

    lista = ["a", "b", "d", "c", "f", "e"]
    seta = set(lista)
    listb = ["a", "c", "e", "f"]
    setb = set(listb)
    print(seta.intersection(setb))
    print(setb.intersection(seta))

    print(rcf.give_one_gp_size('CONTACT'))
    print(rcf.give_one_gp_size('NS-CONTACT'))


#     for i in range (len(rcf.ele_rcf)):
#         print rcf.ele_rcf [i]
#         print '======================='
#     print len(rcf.nod_rcf)
#     print len(rcf.nod_rcf_v)
#     print len(rcf.nod_rcf_a)
#     print len(rcf.nod_rcf_modes)
# ##    print rcf.nod_rcf
# ##    print rcf.give_col_for_ele_rsl ('VOLUMICS','NINT','')
# ##    print rcf.give_col_for_ele_rsl ('VOLUMICS','STRESESS','XX')
# ##    print rcf.give_col_for_ele_rsl ('VOLUMICS','STRAINS','XY')
# ##    print rcf.give_col_for_ele_rsl ('VOLUMICS','HARD_GAMMA','')
# ##    print rcf.give_col_for_ele_rsl ('VOLUMICS','EPSVP_EQ','')
# ##    print rcf.give_col_for_ele_rsl ('VOLUMICS','PLA_CODE','')
# ##
# ##    print rcf.give_col_for_nod_rsl ('PPRESS','')
# ##    print rcf.give_col_for_nod_rsl ('DISP_TRA','Y')
#
#     #print rcf.give_one_gp_size ('VOLUMICS')
#     #print rcf.give_one_node_size ()

if __name__ == '__main__':
    main()
