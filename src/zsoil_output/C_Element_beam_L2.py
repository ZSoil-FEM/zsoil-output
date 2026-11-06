"""Represent standard two-node beam elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_beam_L2` model stores beam connectivity, geometry, material information,
result-file locations, and layer metadata read from a ZSoil project.
"""

from .C_Material import *
from .C_Exf import *
from .C_Element import *
from .C_Node import *
import numpy


class Element_beam_L2 (Element):
    def __init__(self, my_mesh):
        """Initialize the element beam L2 object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)

        # self.centroids = self.nodes [2:3]
        # self.directors = self.nodes [4:5]
        # del self.nodes [2:]
        self.group = "BEAMS"
        self.rsl_ext = ".s04"
        self.rsl_lay_ext = ".l04"
        self.xsiGP = [[0.0]]
        self.nen = 2

    def instanciate(self, my_mesh):
        """Load this object's data from the associated ZSoil project input.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.instanciate(self, my_mesh)

        self.centroids = self.nodes[2:4]
        self.directors = self.nodes[4:6]
        del self.nodes[2:]

    def get_xyz_gp(self, igaus):
        """Return global coordinates for a requested integration point.

        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The global coordinates for a requested integration point.
        :rtype: list
        """
        xyz = [0.0, 0.0, 0.0]
        for i in range(self.nen):
            node = self.mesh_ref.nodes[self.nodes[i]-1]
            for j in range(len(node.xyz)):
                xyz[j] = xyz[j] + node.xyz[j] / self.nen
        return xyz

    def get_ele_dim(self):
        """Return the element dimensions derived from its node coordinates.

        :return: The element dimensions derived from its node coordinates.
        :rtype: float
        """
        node1 = self.mesh_ref.nodes[self.nodes[0]-1]
        node2 = self.mesh_ref.nodes[self.nodes[1]-1]
        vec = np.zeros(3)
        for i in range(3):
            vec[i] = node2.xyz[i] - node1.xyz[i]
        norm = np.linalg.norm(vec)
        return norm

# def get_nlayers (self):
##        mat_index = self.material_index
##        material  = self.mesh_ref.materials [mat_index-1]
##        main_data = material.data [material.dict['MAIN']]
##        nlayers_all    = int(main_data [0])
##        nlayers_reinf  = int(main_data [3])
# return nlayers_all,nlayers_reinf

    # fixed by MP (18.09.2026): when both linear and nonlinear
    # beams exist, SDK erroneously read layers also for elastic beams
    def get_nlayers(self):
        """Return the element layer count and related layer metadata.

        :return: The element layer count and related layer metadata.
        :rtype: tuple[int, int]
        """
        mat_index = self.material_index
        material = self.mesh_ref.materials[mat_index - 1]
        if material.id != 'LAYERED_B':
            return 0, 0
        main_data = material.data[material.dict['MAIN']]
        nlayers_all = int(main_data[0])
        nlayers_reinf = int(main_data[3])
        return nlayers_all, nlayers_reinf
