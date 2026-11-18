"""NumPy port of Ref_Element_EvaluateLocCoordForGloCoord_ex.

Coordinates are node-major, unlike the Fortran arrays. The Newton and
Gauss-Newton equations (including the optional global-coordinate constraint)
follow the supplied Fortran. Failed solves have an explicit status; simplex
membership also checks the sum of barycentric coordinates.
"""
from dataclasses import dataclass
import numpy as np

OUTSIDE_ROUGH_CHECK = 0
OUTSIDE_LOCAL_COORDINATES = 1
OUTSIDE_DISTANCE = 2
INSIDE = 3
PROCEDURE_FAILURE = -1


@dataclass
class LocalCoordinateResult:
    status: int
    xsi: np.ndarray
    point: np.ndarray
    distance: float
    iterations: int


class ReferenceGeometry:
    """Reference coordinates and topology for SDK linear continuum elements.

    Q4/B8 shape functions are reused from the SDK. T3, TH4 and W6 follow
    QS.FOR (the SDK T3 methods are currently shadowed by duplicate definitions).
    L2 is provided for 2D face projection.
    """
    def __init__(self, label):
        self.label = {'Q4ES': 'Q4', 'B8ES': 'B8'}.get(label, label)
        centers = {'L2': [0.], 'T3': [1/3, 1/3], 'Q4': [0., 0.],
                   'TH4': [1/3, 1/3, 1/3], 'W6': [1/3, 1/3, 0.],
                   'B8': [0., 0., 0.]}
        if self.label not in centers:
            raise ValueError('Unsupported element geometry: %s' % label)
        self.center = np.array(centers[self.label])
        self.dimension = len(self.center)
        self.nen = {'L2': 2, 'T3': 3, 'Q4': 4, 'TH4': 4, 'W6': 6, 'B8': 8}[self.label]
        self.reference = None
        if self.label == 'Q4':
            from C_Element_Q4_generic import Element_Q4_generic
            self.reference = Element_Q4_generic()
        elif self.label == 'B8':
            from C_Element_B8_generic import Element_B8_generic
            self.reference = Element_B8_generic()

    def shape(self, xsi):
        if self.reference is not None:
            return self.reference.get_N(xsi), self.reference.get_DN(xsi)
        if self.label == 'L2':
            return np.array([(1-xsi[0])/2, (1+xsi[0])/2]), np.array([[-.5], [.5]])
        if self.label in ('T3', 'TH4'):
            return np.r_[1-np.sum(xsi), xsi], np.vstack((-np.ones(self.dimension), np.eye(self.dimension)))
        r, s, t = xsi
        triangle = np.array([1-r-s, r, s])
        derivative = np.array([[-1., -1.], [1., 0.], [0., 1.]])
        n = np.r_[triangle*(1-t)/2, triangle*(1+t)/2]
        dn = np.vstack((np.column_stack((derivative*(1-t)/2, -triangle/2)),
                        np.column_stack((derivative*(1+t)/2, triangle/2))))
        return n, dn

    def inside(self, xsi, tolerance):
        if not np.all(np.isfinite(xsi)):
            return False
        if self.label in ('T3', 'TH4'):
            return bool(np.min(xsi) >= -tolerance and np.sum(xsi) <= 1+tolerance)
        if self.label == 'W6':
            return bool(np.min(xsi[:2]) >= -tolerance and np.sum(xsi[:2]) <= 1+tolerance
                        and abs(xsi[2]) <= 1+tolerance)
        return bool(np.max(np.abs(xsi)) <= 1+tolerance)

    def faces(self):
        """Zero-based local node numbers, matching Ref_element.for."""
        return {
            'T3': [(0, 1), (1, 2), (2, 0)],
            'Q4': [(0, 1), (1, 2), (2, 3), (3, 0)],
            'TH4': [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)],
            'W6': [(0, 2, 1), (3, 4, 5), (0, 3, 5, 2), (4, 1, 2, 5), (0, 1, 4, 3)],
            'B8': [(0, 3, 2, 1), (0, 1, 5, 4), (0, 4, 7, 3),
                   (1, 2, 6, 5), (2, 3, 7, 6), (4, 5, 6, 7)],
        }[self.label]


def evaluate_local_coordinates(point, geometry, coordinates, tol_xsi=1.e-6,
                               rough_check=True, tol_distance=1.e-10,
                               constrained_global_coord=0, constraint_value=0.,
                               max_iterations=10):
    """Return status, xsi, projected point, distance and iteration count.

    ``geometry`` is a ReferenceGeometry or a supported label. For continuum
    elements use coordinates[:, :mesh.ndm]; use three columns for surfaces in
    3D. tol_distance is relative to the largest bounding-box extent, as in
    Fortran. constrained_global_coord is 0 (disabled) or a 1-based global axis;
    constraints apply only when global and local dimensions differ.
    """
    g = ReferenceGeometry(geometry) if isinstance(geometry, str) else geometry
    x = np.asarray(coordinates, dtype=float)
    p = np.asarray(point, dtype=float)
    if (x.ndim != 2 or x.shape[0] != g.nen or x.shape[1] not in (2, 3)
            or x.shape[1] < g.dimension or p.shape != (x.shape[1],)
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(p))):
        raise ValueError('Expected finite node-major coordinates and a matching point')
    if (not np.isfinite(tol_xsi) or tol_xsi < 0 or not np.isfinite(tol_distance)
            or tol_distance < 0 or not isinstance(max_iterations, int) or max_iterations < 1
            or constrained_global_coord not in range(x.shape[1]+1)
            or not np.isfinite(constraint_value)):
        raise ValueError('Invalid solver options')
    lo, hi = x.min(axis=0), x.max(axis=0)
    scale = float(np.max(hi-lo))
    xsi = g.center.copy()
    # Translate before multiplying: avoid cancellation for large world origins.
    origin = x[0].copy()
    x = x-origin
    target = p-origin
    def result(status, iterations):
        mapped = g.shape(xsi)[0] @ x
        return LocalCoordinateResult(status, xsi.copy(), mapped+origin,
                                     float(np.linalg.norm(mapped-target)), iterations)
    if scale == 0:
        return result(PROCEDURE_FAILURE, 0)
    padding = np.maximum((hi-lo)*.01, scale*.001)
    if rough_check and (np.any(p < lo-padding) or np.any(p > hi+padding)):
        return result(OUTSIDE_ROUGH_CHECK, 0)
    projection = x.shape[1] != g.dimension
    constraint = projection and constrained_global_coord != 0
    multiplier = 0.
    converged = False
    for iteration in range(1, max_iterations+1):
        n, dn = g.shape(xsi)
        residual = n @ x-target
        jacobian = x.T @ dn
        if projection:
            lhs = jacobian.T @ jacobian
            rhs = -jacobian.T @ residual
            if constraint:
                axis = constrained_global_coord-1
                row = jacobian[axis]
                rhs = np.r_[rhs-multiplier*row,
                            constraint_value-origin[axis]-(n @ x)[axis]]
                lhs = np.block([[lhs, row[:, None]], [row[None, :], np.zeros((1, 1))]])
        else:
            lhs, rhs = jacobian, -residual
        try:
            step = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError:
            return result(PROCEDURE_FAILURE, iteration)
        if not np.all(np.isfinite(step)):
            return result(PROCEDURE_FAILURE, iteration)
        xsi += step[:g.dimension]
        if constraint:
            multiplier += step[-1]
        if np.max(np.abs(step[:g.dimension])) <= 1.e-14 and (not constraint or iteration > 1):
            converged = True
            break
    if not converged:
        return result(PROCEDURE_FAILURE, iteration)
    if not g.inside(xsi, tol_xsi):
        return result(OUTSIDE_LOCAL_COORDINATES, iteration)
    answer = result(INSIDE, iteration)
    # Recompute the residual at the final xsi (Fortran uses the previous iterate).
    residual_limit = max(tol_distance, 1.e-15)*scale
    if answer.distance > residual_limit:
        answer.status = OUTSIDE_DISTANCE
    if constraint:
        if abs(answer.point[axis]-constraint_value) > residual_limit:
            answer.status = OUTSIDE_DISTANCE
    return answer
