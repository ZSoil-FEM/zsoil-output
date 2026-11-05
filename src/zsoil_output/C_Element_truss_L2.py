"""Represent two-node truss elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_truss_L2` model stores truss connectivity, geometry, material and
existence data, integration points, and result-file metadata.
"""

from .C_Material import *
from .C_Exf import *
from .C_Element import *


class Element_truss_L2 (Element):
    def __init__(self, my_mesh):
        """Initialize the element truss L2 object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)
        self.group = "TRUSSES"
        self.rsl_ext = ".s03"
        self.xsiGP = [[0.0]]
        self.nen = 2

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
