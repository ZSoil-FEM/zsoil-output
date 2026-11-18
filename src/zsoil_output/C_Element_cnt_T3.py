"""Represent three-node, three-dimensional contact elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_cnt_T3` model provides interface geometry, integration-point data, local
coordinate transformations, connectivity, and result-file metadata.
"""

from .C_Material import *
from .C_Exf import *
from .C_Element import *


class Element_cnt_T3 (Element):
    def __init__(self, my_mesh):
        """Initialize a three-node triangular contact element for the supplied mesh.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)
        self.group = "CONTACT"
        self.rsl_ext = ".s07"
        self.xsiGP = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
        self.nen = 3
        self.ref_ele = Element_T3_generic()
        self.user_vec = None
        self.T_TNS = None
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

    def set_transf_matrices(self, user_vec):
        """Build the element transformation matrices using the supplied orientation vector.

        :param user_vec: User-supplied orientation vector in global coordinates.
        :type user_vec: Sequence[float]
        """
        self.user_vec = user_vec
        self.T_VEC = None
        self.T_TNS = None

        ele_coord = self.get_ele_coord()

        self.T_LG0 = self.ref_ele.get_T_LG0(ele_coord)
        self.T_GL0 = self.ref_ele.get_T_GL0(ele_coord)

        self.T_VEC = self.ref_ele.get_T_user_VEC(self.T_GL0, user_vec, axis='T')
        if isinstance(self.T_VEC, np.ndarray):
            self.T_TNS = self.ref_ele.get_T_user_TNS(self.T_VEC)

    def get_T_TNS(self):
        """Return the tensor-component transformation matrix.

        :return: The tensor-component transformation matrix.
        :rtype: numpy.ndarray
        """
        return self.T_TNS

    def get_T_VEC(self):
        """Return the vector-component transformation matrix.

        :return: The vector-component transformation matrix.
        :rtype: numpy.ndarray
        """
        return self.T_VEC

    def get_T_GL(self):
        """Return the global-to-local coordinate transformation matrix.

        :return: The global-to-local coordinate transformation matrix.
        :rtype: numpy.ndarray
        """
        return self.T_GL0

    def get_T_LG(self):
        """Return the local-to-global coordinate transformation matrix.

        :return: The local-to-global coordinate transformation matrix.
        :rtype: numpy.ndarray
        """
        return self.T_LG0

    def get_surface(self, quadr_enum):
        """Return the element surface area using the requested quadrature rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :return: The element surface area using the requested quadrature rule.
        :rtype: float
        """
        ele_coord = self.get_ele_coord()
        return self.ref_ele.get_surface(ele_coord, quadr_enum)
