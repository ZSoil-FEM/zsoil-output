"""Provide geometry and interpolation operations for generic four-node quadrilaterals.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_Q4_generic` implementation supplies shape functions, derivatives,
quadrature data, natural coordinates, local bases, and surface-related operations.
"""

from .C_Element_Surf_generic import *
import numpy as np


class Element_Q4_generic(Element_Surf_generic):
    def __init__(self):
        """Initialize the element Q4 generic object for the supplied project data."""
        # Element_Surf_generic.__init__(self)
        self.nen = 4
        self.xsi_gp_1 = [0.0, 0.0]
        self.W_gp_1 = [4.0]
        sqrt_3_by_3 = 1.7320508075688772 / 3.0
        self.xsi_gp_4 = [[-sqrt_3_by_3, -sqrt_3_by_3], [sqrt_3_by_3, -sqrt_3_by_3], [sqrt_3_by_3, sqrt_3_by_3],
                         [-sqrt_3_by_3, sqrt_3_by_3]]
        self.W_gp_4 = [1.0, 1.0, 1.0, 1.0]

    def get_quadr(self, quadr_enum):
        """Return quadrature data for the requested integration rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :return: The quadrature data for the requested integration rule.
        :rtype: tuple
        """
        if quadr_enum == self.QUADR_STD:
            return 4, self.xsi_gp_4, self.W_gp_4
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
            return self.xsi_gp_4[igaus]
        else:
            return None

    def get_N(self, xsi):
        """Return shape-function values at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function values at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        s = xsi[0]
        t = xsi[1]

        N = np.zeros(4)

        N[0] = (1.0 - s) * (1.0 - t) * 0.25
        N[1] = (1.0 + s) * (1.0 - t) * 0.25
        N[2] = (1.0 + s) * (1.0 + t) * 0.25
        N[3] = (1.0 - s) * (1.0 + t) * 0.25

        return N

    def get_N0(self):
        """Return shape-function values at the element center.

        :return: The shape-function values at the element center.
        :rtype: numpy.ndarray
        """
        N = np.zeros(4)
        N[0] = 0.25
        N[1] = 0.25
        N[2] = 0.25
        N[3] = 0.25
        return N

    def get_DN(self, xsi):
        """Return shape-function derivatives at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function derivatives at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        s = xsi[0]
        t = xsi[1]

        DN = np.zeros(8).reshape(4, 2)

        DN[0, 0] = -0.25 * (1.0 - t)
        DN[0, 1] = -0.25 * (1.0 - s)
        DN[1, 0] = 0.25 * (1.0 - t)
        DN[1, 1] = -0.25 * (1.0 + s)
        DN[2, 0] = 0.25 * (1.0 + t)
        DN[2, 1] = 0.25 * (1.0 + s)
        DN[3, 0] = -0.25 * (1.0 + t)
        DN[3, 1] = 0.25 * (1.0 - s)

        return DN

    def get_DN0(self):
        """Return shape-function derivatives at the element center.

        :return: The shape-function derivatives at the element center.
        :rtype: numpy.ndarray
        """
        DN = np.zeros(8).reshape(4, 2)
        DN[0, 0] = -0.25
        DN[0, 1] = -0.25
        DN[1, 0] = 0.25
        DN[1, 1] = -0.25
        DN[2, 0] = 0.25
        DN[2, 1] = 0.25
        DN[3, 0] = -0.25
        DN[3, 1] = 0.25

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
        return 4

    def get_face_node_indices(self, face_index_1):
        """Return the local-node indices of a requested element face.

        :param face_index_1: 1-based face index.
        :type face_index_1: int
        :return: The local-node indices of a requested element face.
        :rtype: list
        """
        n1 = face_index_1
        n2 = n1 + 1
        if n2 > 4:
            n2 = 1
        return [n1, n2]


def main():
    """Run this module's example or command-line workflow."""
    q4 = Element_Q4_generic()
    ele_coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    thick_at_nodes = np.array([1.0, 1.0, 2.0, 2.0])

    S = q4.get_surface(ele_coord, q4.QUADR_STD)
    V = q4.get_volume(ele_coord, thick_at_nodes, q4.QUADR_STD)
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
