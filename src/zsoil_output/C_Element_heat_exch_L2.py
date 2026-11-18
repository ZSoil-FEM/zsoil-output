"""Represent two-node heat-exchange elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_heat_exch_L2` model stores line connectivity, material and existence
data, integration-point information, and heat-result file locations.
"""

from .C_Material import *
from .C_Exf import *
from .C_Element import *


class Element_heat_exch_L2 (Element):
    def __init__(self, my_mesh):
        """Initialize the element heat exch L2 object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)
        self.group = "HEAT_EXCH"
        self.rsl_ext = ".s22"
        self.xsiGP = [0.0]
        self.nen = 4

    def get_xyz_gp(self, igaus):
        """Return global coordinates for a requested integration point.

        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The global coordinates for a requested integration point.
        :rtype: Sequence[float]
        """
        return self.get_center()
