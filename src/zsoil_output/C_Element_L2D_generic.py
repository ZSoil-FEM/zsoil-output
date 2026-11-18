"""Provide geometry and interpolation operations for generic two-node line elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_L2D_generic` implementation supplies line shape functions, derivatives,
quadrature data, natural coordinates, and length-related operations shared by line
formulations.
"""

from .C_Element_Surf_generic import *
import numpy as np


class Element_L2D_generic(Element_Surf_generic):
    def __init__(self):
        """Initialize the element L2 D generic object for the supplied project data."""
        # Element_Surf_generic.__init__(self)
        self.nen = 2
        self.xsi_gp_1 = [0.0]
        self.W_gp_1 = [2.0]
        sqrt_3_by_3 = 1.7320508075688772 / 3.0
        self.xsi_gp_2 = [-sqrt_3_by_3, sqrt_3_by_3]
        self.W_gp_2 = [1.0, 1.0]

    def get_quadr(self, quadr_enum):
        """Return quadrature data for the requested integration rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :return: The quadrature data for the requested integration rule.
        :rtype: tuple
        """
        if quadr_enum == self.QUADR_STD:
            return 2, self.xsi_gp_2, self.W_gp_2
        elif quadr_enum == self.QUADR_CENTRAL:
            return 1, self.xsi_gp_1, self.W_gp_1

    def get_xsi(self, quadr_enum, igaus):
        """Return natural coordinates for a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The natural coordinates for a requested quadrature point.
        :rtype: Optional[Sequence[float]]
        """
        if quadr_enum == self.QUADR_CENTRAL:
            return self.xsi_gp_1
        elif quadr_enum == self.QUADR_STD:
            return self.xsi_gp_2[igaus]
        else:
            return None

    def get_N(self, xsi):
        """Return shape-function values at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function values at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        s = xsi
        N = np.zeros(2)

        N[0] = (1.0 - s) / 2.0
        N[1] = (1.0 + s) / 2.0

        return N

    def get_N0(self):
        """Return shape-function values at the element center.

        :return: The shape-function values at the element center.
        :rtype: numpy.ndarray
        """
        N = np.zeros(2)
        N[0] = 0.5
        N[1] = 0.5
        return N

    def get_DN(self, xsi):
        """Return shape-function derivatives at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function derivatives at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        s = xsi
        DN = np.zeros(2).reshape(2, 1)

        DN[0][0] = -0.5
        DN[1][0] = 0.5

        return DN

    def get_DN0(self):
        """Return shape-function derivatives at the element center.

        :return: The shape-function derivatives at the element center.
        :rtype: numpy.ndarray
        """
        DN = np.zeros(2).reshape(2, 1)

        DN[0][0] = -0.5
        DN[1][0] = 0.5

        return DN

    def get_DN_ex(self, quadr_enum, igaus):
        """Return shape-function derivatives at a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The shape-function derivatives at a requested quadrature point.
        :rtype: numpy.ndarray
        """
        if quadr_enum == self.QUADR_CENTRAL:
            return self.get_DN0()
        else:
            xsi = self.get_xsi(quadr_enum, igaus)
            return self.get_DN(xsi)

    def get_N_ex(self, quadr_enum, igaus):
        """Return shape-function values at a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The shape-function values at a requested quadrature point.
        :rtype: numpy.ndarray
        """
        if quadr_enum == self.QUADR_CENTRAL:
            return self.get_N0()
        else:
            xsi = self.get_xsi(quadr_enum, igaus)
            return self.get_N(xsi)

    def get_nr_of_faces(self):
        """Return the number of faces defined for the reference element.

        :return: The number of faces defined for the reference element.
        :rtype: int
        """
        return 1

    def get_face_node_indices(self, face_index_1):
        """Return the local-node indices of a requested element face.

        :param face_index_1: 1-based face index.
        :type face_index_1: int
        :return: The local-node indices of a requested element face.
        :rtype: list
        """
        return [1, 2]

    def get_ndmG(self):
        """Return the number of global spatial dimensions.

        :return: The number of global spatial dimensions.
        :rtype: int
        """
        return 3

    def get_ndmL(self):
        """Return the number of local element dimensions.

        :return: The number of local element dimensions.
        :rtype: int
        """
        return 1


def main():
    """Run this module's example or command-line workflow."""
    l2 = Element_L2D_generic()
    ele_coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    S = l2.get_surface(ele_coord, l2.QUADR_STD)
    thick_at_nodes = np.zeros(2)
    thick_at_nodes.fill(1.0)
    V = l2.get_volume(ele_coord, thick_at_nodes, l2.QUADR_STD)
    print(S, V)

    # #print ele_coord
    # T_LG = q4.get_T_LG0 (ele_coord)
    # print T_LG
    # T_GL = q4.get_T_GL0 (ele_coord)
    # print T_GL
    #
    # user_vec = np.array ([1.0,1.0,0.0])
    # T_VEC = q4.get_T_user_VEC (T_GL,user_vec,'T')
    # print T_VEC
    # print T_VEC.dot (user_vec)
    # print q4.get_T_user_TNS (T_VEC)


if __name__ == '__main__':
    main()
