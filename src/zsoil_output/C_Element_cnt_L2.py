"""Represent two-node, two-dimensional contact elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_cnt_L2` model provides contact-line geometry, integration-point data,
connectivity, and result-file metadata read from a ZSoil project.
"""

from .C_Material import *
from .C_Exf import *
from .C_Element import *


class Element_cnt_L2 (Element):
    def __init__(self, my_mesh):
        """Initialize a two-node contact element for the supplied mesh.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)
        self.group = "CONTACT"
        self.rsl_ext = ".s07"
        self.xsiGP = [[-1.0], [1.0]]
        self.nen = 2

    def get_xyz_gp(self, igaus):
        """Return global coordinates for a requested integration point.

        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The global coordinates for a requested integration point.
        :rtype: list
        """
        node = self.mesh_ref.nodes[self.nodes[igaus]]

        xyz = [0.0, 0.0, 0.0]
        for i in range(len(self.xsiGP)):
            xyz[i] = node.xyz[i]

        return xyz
