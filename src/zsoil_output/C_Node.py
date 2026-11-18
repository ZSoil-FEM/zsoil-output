"""Represent a mesh node and its coordinates and connectivity.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Node` model stores a 1-based ZSoil node index, global coordinates, and
references to adjacent elements used by mesh-selection and recovery operations.
"""

class Node ():
    def __init__(self, my_mesh):
        """Initialize the node object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        self.mesh_ref = my_mesh
        self.index = -1
        self.load_records = []
        self.xyz = []
        self.adj_elements = []

    def instanciate(self, f):
        """Load this object's data from the associated ZSoil project input.

        :param f: Open input stream positioned at the record to read.
        :type f: IO
        """
        line = f.readline()
        texts = line.split()
        self.index = int(texts[0])
        offs = 1
        n_load_records = 0
        if "L_REC" in line:
            n_load_records = int(texts[len(texts)-1])
        n_xyz = len(texts)-1-(n_load_records//max(1, n_load_records)) * 2

        for i in range(n_xyz):
            xi = float(texts[offs+i])
            self.xyz.append(xi)

        if len(self.xyz) < 3:
            self.xyz.append(0.0)

        if n_load_records > 0:
            line = f.readline()
            texts = line.split()
            for text in texts:
                self.load_records.append(int(text))

    def get_xyz(self):
        """Return this node's global coordinates.

        :return: This node's global coordinates.
        :rtype: Sequence[float]
        """
        return self.xyz

    def attach_element_to_adj_list(self, element):
        """Attach an element to this node's adjacency list.

        :param element: Element used by this operation.
        :type element: Element
        """
        self.adj_elements.append(element)

    def get_adj_elements(self, ele_group_index=None, at_time=None):
        """Return elements adjacent to this node, optionally filtered by group and time.

        :param ele_group_index: Numeric element-group identifier.
        :type ele_group_index: Optional[int]
        :param at_time: Physical time at which the model state is evaluated.
        :type at_time: Optional[float]
        :return: The elements adjacent to this node, optionally filtered by group and time.
        :rtype: list
        """
        aux = []
        for element in self.adj_elements:
            take_it = True
            if ele_group_index != None:
                if element.get_group_index() != ele_group_index:
                    take_it = False
            if at_time != None:
                if not element.is_ON(at_time):
                    take_it = False
            if take_it:
                aux.append(element)
        return aux

    def clear_adj_elements(self):
        """Remove all elements from this node's adjacency list."""
        self.adj_elements = []
