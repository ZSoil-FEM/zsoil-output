"""Point location using the SDK local-coordinate solver.

Shell normals follow cross(dx/dxi, dx/deta) at the projected point.
SHQ4 uses all eight nodes as a volumic B8, rather than a surface.
"""
from dataclasses import dataclass
import numpy as np
from C_LocalCoordinates import (ReferenceGeometry, evaluate_local_coordinates,
                                INSIDE, OUTSIDE_ROUGH_CHECK, OUTSIDE_DISTANCE,
                                PROCEDURE_FAILURE)


@dataclass
class PointLocation:
    inside: bool
    element: object
    xsi: object
    distance: object
    projected_point: object
    normal: object
    status: int

    @property
    def element_index(self):
        return None if self.element is None else self.element.index

    @property
    def xi(self):
        return None if self.xsi is None else float(self.xsi[0])

    @property
    def eta(self):
        return None if self.xsi is None or len(self.xsi) < 2 else float(self.xsi[1])

    @property
    def zeta(self):
        return None if self.xsi is None or len(self.xsi) < 3 else float(self.xsi[2])


def _point(mesh, point):
    if mesh.ndm not in (2, 3):
        raise ValueError('Expected mesh.ndm equal to 2 or 3')
    p = np.asarray(point, dtype=float)
    if mesh.ndm == 2 and p.shape == (3,):
        if p[2] != 0:
            raise ValueError('A 2D point must have zero z coordinate')
        p = p[:2]
    if p.shape != (mesh.ndm,) or not np.all(np.isfinite(p)):
        raise ValueError('Expected a finite point matching mesh.ndm')
    return p


def _geometry(mesh, element):
    if element.group == 'SHELLS' and element.label == 'SHQ4':
        if mesh.ndm != 3 or len(element.nodes) != 8:
            raise ValueError('SHQ4 requires eight nodes and a 3D mesh')
        x = np.array([mesh.nodes[i-1].xyz for i in element.nodes], dtype=float)
        return ReferenceGeometry('B8'), x, False
    if element.group == 'VOLUMICS':
        g = ReferenceGeometry(element.label)
        if g.dimension != mesh.ndm:
            raise ValueError('Element dimension differs from mesh.ndm')
        return g, np.asarray(element.get_ele_coord(), dtype=float)[:, :mesh.ndm], False
    if (mesh.ndm == 3 and element.group in ('SHELLS', 'MEMBRANES', 'CONTACT')
            and element.nen in (3, 4)):
        g = ReferenceGeometry('T3' if element.nen == 3 else 'Q4')
        return g, np.asarray(element.get_ele_coord(), dtype=float), True
    raise ValueError('Unsupported point-location element: %s' % element.label)


def evaluate_point(mesh, point, element, signed_distance=False,
                   tol_xsi=1.e-6, tol_distance=1.e-10, max_iterations=30):
    """Evaluate a point in one SDK element (object or 1-based SDK index).

    Surface distance is unsigned by default; signed_distance=True uses the
    local normal at the projection. Membership requires the projected local
    coordinates to lie in the element and abs(distance) <= tol_distance*size,
    where size is its largest bounding-box extent. No shell thickness test is
    performed. Continuum/SHQ4 distance is 0; the third local coordinate is zeta.
    """
    p = _point(mesh, point)
    if isinstance(element, (int, np.integer)):
        matches = [e for e in mesh.elements if e.index == element]
        if len(matches) != 1:
            raise ValueError('Unknown or ambiguous SDK element index: %s' % element)
        element = matches[0]
    if getattr(element, 'mesh_ref', None) is not mesh:
        raise ValueError('Element does not belong to this mesh')
    g, x, surface = _geometry(mesh, element)
    # Disable the Fortran rough check so off-surface points still produce a
    # projection/distance, even well outside a flat surface's bounding box.
    mapping = evaluate_local_coordinates(p, g, x, tol_xsi=tol_xsi,
                                         rough_check=False, tol_distance=tol_distance,
                                         max_iterations=max_iterations)
    if mapping.status == PROCEDURE_FAILURE:
        return PointLocation(False, element, None, None, None, None, mapping.status)
    normal = None
    distance = 0.
    if surface:
        dn = g.shape(mapping.xsi)[1]
        tangents = (x-x[0]).T @ dn
        normal = np.cross(tangents[:, 0], tangents[:, 1])
        norm = np.linalg.norm(normal)
        if norm == 0 or not np.isfinite(norm):
            return PointLocation(False, element, None, None, None, None, PROCEDURE_FAILURE)
        normal = normal/norm
        # Use the translated solver distance to avoid cancellation at large
        # world origins; the orientation only supplies its sign.
        distance = mapping.distance
        if signed_distance and np.dot(p-mapping.point, normal) < 0:
            distance = -distance
    return PointLocation(mapping.status == INSIDE, element, mapping.xsi,
                         float(distance), mapping.point, normal, mapping.status)


def locate_point(mesh, point, at_time, elements=None, signed_distance=False,
                 tol_xsi=1.e-6, tol_distance=1.e-10, max_iterations=30):
    """Search the shared grid, then test active elements at physical at_time.

    A containing volume takes priority. Otherwise return the nearest active
    surface with an in-domain projection, even if its distance makes it outside.
    No match returns element/xsi/distance=None. Boundaries tie by SDK index.
    """
    from C_MeshSpatialGrid import get_mesh_grid, element_is_active, validate_time
    p = _point(mesh, point)
    validate_time(at_time)
    if (not np.isfinite(tol_xsi) or tol_xsi < 0 or not np.isfinite(tol_distance)
            or tol_distance < 0 or not isinstance(max_iterations, int) or max_iterations < 1):
        raise ValueError('Invalid point-location options')
    grid = get_mesh_grid(mesh)
    allowed = None
    if elements is not None:
        allowed = set()
        for e in elements:
            if getattr(e, 'mesh_ref', None) is not mesh:
                raise ValueError('Element does not belong to this mesh')
            _geometry(mesh, e)  # Preserve explicit-selection validation.
            allowed.add(e.index)
    padding = max((4*tol_xsi+tol_distance)*grid.max_element_size,
                  64*np.finfo(float).eps*grid.max_element_size)
    def evaluate(record):
        e = record.element
        if not element_is_active(e, at_time):
            return None
        return evaluate_point(mesh, p, e, signed_distance=signed_distance,
                              tol_xsi=tol_xsi, tol_distance=tol_distance,
                              max_iterations=max_iterations)
    candidates, failures = [], []
    for record in grid.query_point(p, padding):
        index = record.element.index
        if index not in grid.volume_ids or (allowed is not None and index not in allowed):
            continue
        result = evaluate(record)
        if result is not None:
            if result.inside:
                candidates.append(result)
            elif result.status == PROCEDURE_FAILURE:
                failures.append(result)
    if candidates:
        return min(candidates, key=lambda r: r.element_index)
    best = None
    for lower_distance, record in grid.near_surfaces(p, allowed):
        if best is not None and lower_distance > abs(best.distance)+padding:
            break
        if record is None:
            continue
        result = evaluate(record)
        if result is None:
            continue
        if result.status == PROCEDURE_FAILURE:
            failures.append(result)
        elif result.inside or result.status == OUTSIDE_DISTANCE:
            if best is None or (not result.inside, abs(result.distance), result.element_index) < (
                    not best.inside, abs(best.distance), best.element_index):
                best = result
    if best is not None:
        return best
    if failures:
        return min(failures, key=lambda r: r.element_index)
    return PointLocation(False, None, None, None, None, None, OUTSIDE_ROUGH_CHECK)
