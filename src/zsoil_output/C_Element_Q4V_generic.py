"""Provide volumetric operations for generic four-node quadrilaterals.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_Q4V_generic` implementation combines surface and volumetric helpers for
quadrilateral formulations used by ZSoil continuum elements.
"""

from .C_Element_Surf_generic import *
from .C_Element_Volumic_generic import *
from .C_Element_Q4_generic import *
import numpy as np


class Element_Q4V_generic (Element_Volumic_generic):
    def __init__(self):
        """Initialize the element Q4 V generic object for the supplied project data."""
        self.Q4S = Element_Q4_generic()

    def get_quadr(self, quadr_enum):
        """Return quadrature data for the requested integration rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :return: The quadrature data for the requested integration rule.
        :rtype: tuple[int, Sequence, Sequence]
        """
        return self.Q4S.get_quadr(quadr_enum)

    def get_xsi(self, quadr_enum, igaus):
        """Return natural coordinates for a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The natural coordinates for a requested quadrature point.
        :rtype: Optional[Sequence[float]]
        """
        return self.Q4S.get_xsi(quadr_enum, igaus)

    def get_N(self, xsi):
        """Return shape-function values at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function values at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        return self.Q4S.get_N(xsi)

    def get_N0(self):
        """Return shape-function values at the element center.

        :return: The shape-function values at the element center.
        :rtype: numpy.ndarray
        """
        return self.Q4S.get_N0()

    def get_DN(self, xsi):
        """Return shape-function derivatives at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function derivatives at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        return self.Q4S.get_DN(xsi)

    def get_DN0(self):
        """Return shape-function derivatives at the element center.

        :return: The shape-function derivatives at the element center.
        :rtype: numpy.ndarray
        """
        return self.Q4S.get_DN0()

    def get_DN_ex(self, quadr_enum, igaus):
        """Return shape-function derivatives at a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The shape-function derivatives at a requested quadrature point.
        :rtype: numpy.ndarray
        """
        return self.Q4S.get_DN_ex(quadr_enum, igaus)

    def get_N_ex(self, quadr_enum, igaus):
        """Return shape-function values at a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        :return: The shape-function values at a requested quadrature point.
        :rtype: numpy.ndarray
        """
        return self.Q4S.get_N_ex(quadr_enum, igaus)

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

    def get_face_info(self, face_index_1):
        """Return the type and local-node indices of a requested element face.

        :param face_index_1: 1-based face index.
        :type face_index_1: int
        :return: The type and local-node indices of a requested element face.
        :rtype: tuple[str, list[int]]
        """
        return 'L2', self.get_face_node_indices(face_index_1)


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
