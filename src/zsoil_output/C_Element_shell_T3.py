"""Represent three-node triangular shell elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_shell_T3` model combines triangular interpolation with shell
connectivity, local coordinate transformations, material data, integration points, and
result-file metadata.
"""


from .C_Element import *
from .C_Element_T3_generic import *


class Element_shell_T3 (Element):
    def __init__(self, my_mesh):
        """Initialize the element shell T3 object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        Element.__init__(self, my_mesh)
        self.group = "SHELLS"
        self.rsl_ext = ".s02"
        self.rsl_lay_ext = ".l02"
        self.xsiGP = [[1.0/3.0, 1.0/3.0]]
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
        xyz = [0.0, 0.0, 0.0]
        for i in range(self.nen):
            node = self.mesh_ref.nodes[self.nodes[i]-1]
            for j in range(len(node.xyz)):
                xyz[j] = xyz[j] + node.xyz[j] / self.nen

        return xyz

    def get_nlayers(self):
        """Return the element layer count and related layer metadata.

        :return: The element layer count and related layer metadata.
        :rtype: tuple[int, int]
        """
        mat_index = self.material_index
        material = self.mesh_ref.materials[mat_index-1]
        main_data = material.data[material.dict['MAIN']]
        nlayers_all = int(main_data[0])
        nlayers_reinf = int(main_data[3])
        return nlayers_all, nlayers_reinf

    def set_transf_matrices(self, user_vec):
        """Build the element transformation matrices using the supplied orientation vector.

        :param user_vec: User-supplied orientation vector in global coordinates.
        :type user_vec: Sequence[float]
        """
        self.user_vec = user_vec
        self.T_VEC = None
        self.T_TNS = None

        ele_coord = self.get_ele_coord()
        T_GL0 = self.ref_ele.get_T_GL0(ele_coord)
        T_VEC = self.ref_ele.get_T_user_VEC(T_GL0, user_vec, axis='T')
        if T_VEC != None:
            T_TNS = self.ref_ele.get_T_user_TNS(T_VEC)

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

    def get_surface(self, quadr_enum):
        """Return the element surface area using the requested quadrature rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :return: The element surface area using the requested quadrature rule.
        :rtype: float
        """
        ele_coord = self.get_ele_coord()
        return self.ref_ele.get_surface(ele_coord, quadr_enum)
