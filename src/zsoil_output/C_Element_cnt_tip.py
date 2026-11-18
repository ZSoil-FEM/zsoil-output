"""Represent contact-tip elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_cnt_tip` model stores the connectivity, material data, integration-point
information, and result-file metadata associated with a ZSoil contact tip.
"""

from .C_Material import *
from .C_Exf import *
from .C_Element import *


class Element_cnt_tip (Element):
    def __init__(self, my_mesh):
        """Initialize a contact-tip element for the supplied mesh.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)
        self.group = "CONTACT"
        self.rsl_ext = ".s07"
        self.xsiGP = [[0.0]]
        self.nen = 1
        self.T_VEC = None

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

    def get_T_VEC(self):
        """Return the vector-component transformation matrix.

        :return: The vector-component transformation matrix.
        :rtype: numpy.ndarray
        """
        return self.T_VEC
