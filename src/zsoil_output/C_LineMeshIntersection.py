"""Finite line / continuum mesh intersections for result extraction.

Based on Ref_element_setGrid, Ref_Element_IntersectWithLine_3D and
Ref_Element_ProjectPointOnFaceInDir. The shared mesh spatial grid supplies
candidate elements for every query. No result files or mesh adjacency lists are changed.
"""
from dataclasses import dataclass
import numpy as np
from C_LocalCoordinates import ReferenceGeometry, evaluate_local_coordinates, INSIDE


@dataclass
class IntersectionSegment:
    element: object
    point1: np.ndarray
    point2: np.ndarray
    xsi1: np.ndarray
    xsi2: np.ndarray
    t1: float
    t2: float

    @property
    def element_index(self):
        return self.element.index

    @property
    def length(self):
        return float(np.linalg.norm(self.point2-self.point1))


@dataclass
class LineMeshIntersection:
    segments: list
    endpoint_elements: tuple

    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return len(self.segments)


def _edge_hits(a, direction, p, q, distance_tol, t_tol):
    """Also handles collinear edges, for a line lying on a face boundary."""
    edge = q-p
    matrix = np.column_stack((direction, -edge))
    sol, _, rank, _ = np.linalg.lstsq(matrix, p-a, rcond=None)
    if rank == 2:
        if (np.linalg.norm(matrix @ sol-(p-a)) <= distance_tol
                and -t_tol <= sol[0] <= 1+t_tol and -1.e-10 <= sol[1] <= 1+1.e-10):
            return [float(np.clip(sol[0], 0., 1.))]
        return []
    norm2 = float(direction @ direction)
    tp = float((p-a) @ direction/norm2)
    tq = float((q-a) @ direction/norm2)
    if max(np.linalg.norm(p-a-tp*direction), np.linalg.norm(q-a-tq*direction)) > distance_tol:
        return []
    lo, hi = max(0., min(tp, tq)), min(1., max(tp, tq))
    return [lo, hi] if hi >= lo else []


def _face_hits(a, direction, x, tol_xsi, distance_tol, t_tol):
    """Fortran directional face projection: [dx/dxsi, -AB] delta = rhs.

    The finite parameter range and final residual are checked explicitly.
    Multiple seeds cover either root on a warped bilinear quadrilateral.
    """
    if len(x) == 2:
        return _edge_hits(a, direction, x[0], x[1], distance_tol, t_tol)
    g = ReferenceGeometry('T3' if len(x) == 3 else 'Q4')
    # Scale the line column by element size, rather than by query length.
    length = np.linalg.norm(direction)
    scale = np.max(np.ptp(x, axis=0))
    if scale == 0:
        return []
    vector = direction/length*scale
    hits = []
    seeds = [g.center]
    if len(x) == 4:
        seeds += [np.array([s, t]) for s in (-.8, .8) for t in (-.8, .8)]
    for seed in seeds:
        xsi, lam = seed.copy(), 0.
        for _ in range(21):
            n, dn = g.shape(xsi)
            lhs = np.column_stack((x.T @ dn, -vector))
            rhs = a+lam*vector-n @ x
            try:
                step = np.linalg.solve(lhs, rhs)
            except np.linalg.LinAlgError:
                break
            if not np.all(np.isfinite(step)):
                break
            xsi += step[:-1]
            lam += step[-1]
            if np.max(np.abs(xsi)) > 10:
                break
            if np.linalg.norm(step[:-1]) <= 1.e-10:
                t = float(lam*scale/length)
                if (g.inside(xsi, tol_xsi) and -t_tol <= t <= 1+t_tol
                        and np.linalg.norm(g.shape(xsi)[0] @ x-a-t*direction) <= distance_tol):
                    hits.append(float(np.clip(t, 0., 1.)))
                break
    # Singular face projection occurs for lines contained in the face. Its
    # boundary edges still locate entry/exit, including lines along an edge.
    for p, q in zip(x, np.roll(x, -1, axis=0)):
        hits.extend(_edge_hits(a, direction, p, q, distance_tol, t_tol))
    return hits


class LineMeshIntersector:
    """Reusable intersector for T3/Q4/Q4ES (2D), B8/B8ES/W6/TH4 (3D).

    Continuum elements and volumic SHQ4 are selected by default; an explicit elements iterable can
    restrict the query (e.g. mesh.get_list_of_elements(..., only_active=True)).
    Geometry uses the shared mesh grid: rebuild it after changing nodes/connectivity.
    Surface/shell crossings are points and are outside this segment API.

    tol_xsi is the dimensionless membership tolerance. tol_distance is relative
    to each element's largest extent. min_length is an absolute physical length.
    Valid, non-inverted linear isoparametric elements are required.
    """
    def __init__(self, mesh, elements=None, tol_xsi=1.e-6,
                 tol_distance=1.e-10, min_length=0.):
        for value in (tol_xsi, tol_distance, min_length):
            if not np.isfinite(value) or value < 0:
                raise ValueError('Tolerances and min_length must be finite and nonnegative')
        self.mesh_ref = mesh
        self.tol_xsi = float(tol_xsi)
        self.tol_distance = float(tol_distance)
        self.min_length = float(min_length)
        self._selection = None if elements is None else tuple(elements)
        self.rebuild(rebuild_mesh_grid=False)

    def rebuild(self, rebuild_mesh_grid=True):
        """Refresh selection; optionally rebuild the shared geometry grid."""
        from C_MeshSpatialGrid import get_mesh_grid, MeshSpatialGrid
        self.dimension = self.mesh_ref.ndm
        if rebuild_mesh_grid:
            old = get_mesh_grid(self.mesh_ref)
            self.mesh_ref.spatial_grid = MeshSpatialGrid(self.mesh_ref, old.max_cells_per_axis)
            self.mesh_ref.bound_box = self.mesh_ref.spatial_grid.bound_box
        grid = get_mesh_grid(self.mesh_ref)
        self._allowed = None
        selected = (grid.records[i].element for i in grid.volume_ids) if self._selection is None else self._selection
        indices = set()
        for e in selected:
            if getattr(e, 'mesh_ref', None) is not self.mesh_ref:
                raise ValueError('Element does not belong to this mesh')
            if e.index not in grid.volume_ids:
                raise ValueError('Explicit selection must contain volumic elements only')
            if e.index in indices:
                raise ValueError('Duplicate element index: %s' % e.index)
            indices.add(e.index)
            g = ReferenceGeometry('B8' if e.label == 'SHQ4' else e.label)
            if g.dimension != self.dimension:
                raise ValueError('Element %s has the wrong dimension' % e.index)
        if self._selection is not None:
            self._allowed = indices

    def _candidate_records(self, a, b, at_time):
        from C_MeshSpatialGrid import get_mesh_grid, element_is_active
        grid = get_mesh_grid(self.mesh_ref)
        padding = max(self.tol_distance, 64*np.finfo(float).eps)*grid.max_element_size
        for record in grid.query_segment(a, b, padding):
            e = record.element
            if e.index not in grid.volume_ids or (self._allowed is not None and e.index not in self._allowed):
                continue
            if not element_is_active(e, at_time):
                continue
            g = ReferenceGeometry('B8' if e.label == 'SHQ4' else e.label)
            if record.size == 0:
                raise ValueError('Degenerate element %s' % e.index)
            pad = max(self.tol_distance, 64*np.finfo(float).eps)*record.size
            yield e, g, record.coordinates, record.lower-pad, record.upper+pad, pad

    def _point(self, value):
        p = np.asarray(value, dtype=float)
        if p.shape == (3,) and self.dimension == 2:
            if p[2] != 0:
                raise ValueError('A 2D line must have zero z coordinate')
            p = p[:2]
        if p.shape != (self.dimension,) or not np.all(np.isfinite(p)):
            raise ValueError('Expected a finite %sD point' % self.dimension)
        return p.copy()

    def intersect(self, a, b, at_time, remove_overlaps=True):
        """Return ordered segments and element indices containing A/B (or None).

        at_time is a required finite physical time; inactive elements are
        excluded before assembling segments and removing overlaps.
        t runs from 0 at A to 1 at B. Partial intersections and a segment wholly
        inside an element are included. Gaps stay gaps; point-only contacts are
        omitted. By default an overlap belongs to the first entered element,
        breaking ties by SDK element index; clipped tails are retained. Set
        remove_overlaps=False to keep a separate segment for each owner.
        """
        from C_MeshSpatialGrid import validate_time
        validate_time(at_time)
        a, b = self._point(a), self._point(b)
        direction = b-a
        length = float(np.linalg.norm(direction))
        if length == 0 or not np.isfinite(length):
            raise ValueError('Line segment must have finite, positive length')
        candidates = []
        endpoints = [None, None]
        for e, g, world, _, _, distance_tol in self._candidate_records(a, b, at_time):
            origin = world[0]
            x, p = world-origin, a-origin
            t_tol = max(64*np.finfo(float).eps, distance_tol/length)
            def local(t):
                return evaluate_local_coordinates(p+t*direction, g, x,
                                                  tol_xsi=self.tol_xsi,
                                                  tol_distance=self.tol_distance)
            for j in (0, 1):
                if local(float(j)).status == INSIDE:
                    endpoints[j] = e.index if endpoints[j] is None else min(endpoints[j], e.index)
            cuts = [0., 1.]
            for face in g.faces():
                cuts.extend(_face_hits(p, direction, x[list(face)], self.tol_xsi,
                                       distance_tol, t_tol))
            cuts.sort()
            unique = [cuts[0]]
            for t in cuts[1:]:
                if t-unique[-1] > t_tol:
                    unique.append(t)
                elif t == 1.:
                    unique[-1] = 1.
            intervals = []
            for t1, t2 in zip(unique, unique[1:]):
                if (t2-t1)*length <= max(self.min_length, distance_tol):
                    continue
                if local((t1+t2)/2).status == INSIDE:
                    if intervals and abs(t1-intervals[-1][1]) <= t_tol:
                        intervals[-1] = (intervals[-1][0], t2)
                    else:
                        intervals.append((t1, t2))
            for t1, t2 in intervals:
                candidates.append((t1, t2, e.index, e, g, x, p, distance_tol))
        candidates.sort(key=lambda item: (item[0], item[2], item[1]))
        segments = []
        front = 0.
        for t1, t2, _, e, g, x, p, distance_tol in candidates:
            if remove_overlaps:
                t1 = max(t1, front)
            if (t2-t1)*length <= max(self.min_length, distance_tol):
                continue
            loc = [evaluate_local_coordinates(p+t*direction, g, x,
                                              tol_xsi=self.tol_xsi,
                                              tol_distance=self.tol_distance) for t in (t1, t2)]
            if any(r.status != INSIDE for r in loc):
                raise RuntimeError('Endpoint inverse mapping failed at element %s' % e.index)
            segments.append(IntersectionSegment(e, a+t1*direction, a+t2*direction,
                                                 loc[0].xsi, loc[1].xsi, t1, t2))
            front = max(front, t2)
        return LineMeshIntersection(segments, tuple(endpoints))
