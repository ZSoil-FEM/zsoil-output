"""Provide geometry and interpolation operations for generic eight-node brick elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_B8_generic` implementation supplies shape functions, derivatives,
quadrature data, natural coordinates, faces, and volume-related operations shared by brick
formulations.
"""

from .C_Element_Volumic_generic import *
import numpy as np


class Element_B8_generic (Element_Volumic_generic):
    def __init__(self):
        """Initialize the element B8 generic object for the supplied project data."""
        self.nen = 8
        self.xsi_gp_1 = [0.0, 0.0, 0.0]
        self.W_gp_1 = [4.0]
        sqrt_3_by_3 = 1.7320508075688772 / 3.0
        self.xsi_gp_8 = [[-sqrt_3_by_3, -sqrt_3_by_3, -sqrt_3_by_3],
                         [sqrt_3_by_3, -sqrt_3_by_3, -sqrt_3_by_3],
                         [sqrt_3_by_3, sqrt_3_by_3, -sqrt_3_by_3],
                         [-sqrt_3_by_3, sqrt_3_by_3, -sqrt_3_by_3],
                         [-sqrt_3_by_3, -sqrt_3_by_3, sqrt_3_by_3],
                         [sqrt_3_by_3, -sqrt_3_by_3, sqrt_3_by_3],
                         [sqrt_3_by_3, sqrt_3_by_3, sqrt_3_by_3],
                         [-sqrt_3_by_3, sqrt_3_by_3, sqrt_3_by_3]]
        self.W_gp_8 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        self.faces = [[1, 4, 3, 2],
                      [1, 2, 6, 5],
                      [1, 5, 8, 4],
                      [2, 3, 7, 6],
                      [3, 4, 8, 7],
                      [5, 6, 7, 8]]

    def get_quadr(self, quadr_enum):
        """Return quadrature data for the requested integration rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :return: The quadrature data for the requested integration rule.
        :rtype: tuple
        """
        if quadr_enum == self.QUADR_STD:
            return 8, self.xsi_gp_8, self.W_gp_8
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
            return self.xsi_gp_8[igaus]
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
        u = xsi[2]
        N = np.zeros(8)

        N[0] = 0.125 * (1.0 - s) * (1.0 - t) * (1.0 - u)
        N[1] = 0.125 * (1.0 + s) * (1.0 - t) * (1.0 - u)
        N[2] = 0.125 * (1.0 + s) * (1.0 + t) * (1.0 - u)
        N[3] = 0.125 * (1.0 - s) * (1.0 + t) * (1.0 - u)
        N[4] = 0.125 * (1.0 - s) * (1.0 - t) * (1.0 + u)
        N[5] = 0.125 * (1.0 + s) * (1.0 - t) * (1.0 + u)
        N[6] = 0.125 * (1.0 + s) * (1.0 + t) * (1.0 + u)
        N[7] = 0.125 * (1.0 - s) * (1.0 + t) * (1.0 + u)

        return N

    def get_N0(self):
        """Return shape-function values at the element center.

        :return: The shape-function values at the element center.
        :rtype: numpy.ndarray
        """
        return self.get_N(self.xsi_gp_1)

    def get_DN(self, xsi):
        """Return shape-function derivatives at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        :return: The shape-function derivatives at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        s = xsi[0]
        t = xsi[1]
        u = xsi[2]

        DN = np.zeros(24).reshape(8, 3)

        DN[0, 0] = -0.125 * (1.0 - t) * (1.0 - u)
        DN[0, 1] = -0.125 * (1.0 - s) * (1.0 - u)
        DN[0, 2] = -0.125 * (1.0 - s) * (1.0 - t)

        DN[1, 0] = 0.125 * (1.0 - t) * (1.0 - u)
        DN[1, 1] = -0.125 * (1.0 + s) * (1.0 - u)
        DN[1, 2] = -0.125 * (1.0 + s) * (1.0 - t)

        DN[2, 0] = 0.125 * (1.0 + t) * (1.0 - u)
        DN[2, 1] = 0.125 * (1.0 + s) * (1.0 - u)
        DN[2, 2] = -0.125 * (1.0 + s) * (1.0 + t)

        DN[3, 0] = -0.125 * (1.0 + t) * (1.0 - u)
        DN[3, 1] = 0.125 * (1.0 - s) * (1.0 - u)
        DN[3, 2] = -0.125 * (1.0 - s) * (1.0 + t)

        DN[4, 0] = -0.125 * (1.0 - t) * (1.0 + u)
        DN[4, 1] = -0.125 * (1.0 - s) * (1.0 + u)
        DN[4, 2] = 0.125 * (1.0 - s) * (1.0 - t)

        DN[5, 0] = 0.125 * (1.0 - t) * (1.0 + u)
        DN[5, 1] = -0.125 * (1.0 + s) * (1.0 + u)
        DN[5, 2] = 0.125 * (1.0 + s) * (1.0 - t)

        DN[6, 0] = 0.125 * (1.0 + t) * (1.0 + u)
        DN[6, 1] = 0.125 * (1.0 + s) * (1.0 + u)
        DN[6, 2] = 0.125 * (1.0 + s) * (1.0 + t)

        DN[7, 0] = -0.125 * (1.0 + t) * (1.0 + u)
        DN[7, 1] = 0.125 * (1.0 - s) * (1.0 + u)
        DN[7, 2] = 0.125 * (1.0 - s) * (1.0 + t)

        return DN

    def get_DN0(self):
        """Return shape-function derivatives at the element center.

        :return: The shape-function derivatives at the element center.
        :rtype: numpy.ndarray
        """
        return self.get_DN(self.xsi_gp_1)

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
        return len(self.faces)

    def get_face_info(self, face_index_1):
        """Return the type and local-node indices of a requested element face.

        :param face_index_1: 1-based face index.
        :type face_index_1: int
        :return: The type and local-node indices of a requested element face.
        :rtype: tuple
        """
        return 'Q4', self.faces[face_index_1-1]
