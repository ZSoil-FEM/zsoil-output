"""Provide shared geometry operations for volumetric elements.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_Volumic_generic` base supplies element-centre, extent, face, and volume
operations used by continuum formulations.
"""

import numpy as np


class Element_Volumic_generic ():
    QUADR_CENTRAL = 0
    QUADR_STD = 1
    QUADR_NODAL = 2

    def __init__(self,):
        """Initialize the element volumic generic object for the supplied project data."""
        self.nen = 0

    def get_quadr(self, quadr_enum):
        """Return quadrature data for the requested integration rule.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        """
        pass

    def get_xsi(self, ngaus, igaus):
        """Return natural coordinates for a requested quadrature point.

        :param ngaus: Ngaus used by this operation.
        :type ngaus: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        """
        pass

    def get_N_ex(self, quadr_enum, igaus):
        """Return shape-function values at a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        """
        pass

    def get_N(self, xsi):
        """Return shape-function values at the requested natural coordinates.

        :param xsi: Natural coordinates inside the reference element.
        :type xsi: Sequence[float]
        """
        pass

    def get_DN0(self):
        """Return shape-function derivatives at the element center."""
        pass

    def get_DN_ex(self, quadr_enum, igaus):
        """Return shape-function derivatives at a requested quadrature point.

        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param igaus: Zero-based integration-point index.
        :type igaus: int
        """
        pass

    def get_volume(self, ele_coord, quadr_enum, ret_V_igaus=False):
        """Return the element volume using the requested quadrature rule.

        :param ele_coord: Element coordinates used by this operation.
        :type ele_coord: numpy.ndarray
        :param quadr_enum: Quadrature-rule identifier.
        :type quadr_enum: int
        :param ret_V_igaus: Ret v igaus used by this operation.
        :type ret_V_igaus: bool
        :return: The element volume using the requested quadrature rule.
        :rtype: float
        """
        ngaus, xsi, Wi = self.get_quadr(quadr_enum)
        Vi = np.zeros(len(xsi))

        V = 0.0

        ndmG = 3
        ndmL = 3

        dxG_dxL = np.zeros(ndmG * ndmL).reshape(ndmG, ndmL)

        for igaus in range(ngaus):
            DN = self.get_DN_ex(quadr_enum, igaus)
            # dx/dxsi
            dxG_dxL.fill(0.0)
            for inode in range(self.nen):
                for idmL in range(ndmL):
                    for idmG in range(ndmG):
                        dxG_dxL[idmG, idmL] = dxG_dxL[idmG, idmL] + \
                            DN[inode, idmL] * ele_coord[inode, idmG]
            aux = norm = np.linalg.det(dxG_dxL)
            V = V + aux * Wi[igaus]
            if ret_V_igaus:
                Vi[igaus] = aux * Wi[igaus]

        if ret_V_igaus:
            return Vi
        else:
            return V
