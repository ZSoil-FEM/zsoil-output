"""Reusable sparse uniform grid containing ALL mesh elements, at ALL times."""
from dataclasses import dataclass
from itertools import product
import heapq
import numpy as np


def element_is_active(element, at_time):
    if at_time is None:
        return True
    exf = element.give_exf()
    return exf is None or bool(exf.is_ON(at_time))


def validate_time(at_time):
    if at_time is None or not np.isscalar(at_time) or not np.isfinite(at_time):
        raise ValueError('A finite physical at_time is required')


def box_interval(a, direction, lower, upper):
    start, end = 0., 1.
    for i in range(len(a)):
        if direction[i] == 0:
            if a[i] < lower[i] or a[i] > upper[i]:
                return None
        else:
            t0, t1 = sorted(((lower[i]-a[i])/direction[i], (upper[i]-a[i])/direction[i]))
            start, end = max(start, t0), min(end, t1)
            if start > end:
                return None
    return start, end


@dataclass
class GridElement:
    element: object
    coordinates: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    size: float


class MeshSpatialGrid:
    """Node bounding box and element bounding boxes in a sparse cubic grid.

    The longest axis uses ceil(n_elements**(1/ndm)) cells, capped by
    max_cells_per_axis (50 by default); flat axes use one cell. Elements are
    inserted into every touched cell, including both sides of grid boundaries.
    Large boxes covering >4096 cells use an overflow list to bound storage.
    Existence functions are deliberately never consulted during construction.
    """
    def __init__(self, mesh, max_cells_per_axis=50):
        if not isinstance(max_cells_per_axis, int) or max_cells_per_axis < 1:
            raise ValueError('max_cells_per_axis must be a positive integer')
        self.mesh_ref = mesh
        self.dimension = mesh.ndm
        if self.dimension not in (2, 3):
            raise ValueError('Expected mesh.ndm equal to 2 or 3')
        self.max_cells_per_axis = max_cells_per_axis
        self.records, self.cells, self.overflow = {}, {}, set()
        self.volume_ids, self.surface_ids = set(), set()
        nodes = np.asarray([n.xyz[:self.dimension] for n in mesh.nodes], dtype=float)
        if len(nodes):
            if nodes.shape != (len(mesh.nodes), self.dimension) or not np.all(np.isfinite(nodes)):
                raise ValueError('Invalid mesh node coordinates')
            self.bound_box = np.array([nodes.min(axis=0), nodes.max(axis=0)])
        else:
            self.bound_box = np.zeros((2, self.dimension))
        for e in mesh.elements:
            if e.index in self.records:
                raise ValueError('Duplicate SDK element index: %s' % e.index)
            if not e.nodes:
                continue
            # SHQ4: all eight volumic nodes. Other SDK element types use nen.
            indices = e.nodes if e.label == 'SHQ4' else e.nodes[:e.nen or len(e.nodes)]
            x = np.array([mesh.nodes[i-1].xyz[:self.dimension] for i in indices], dtype=float)
            lo, hi = x.min(axis=0), x.max(axis=0)
            self.records[e.index] = GridElement(e, x, lo, hi, float(np.max(hi-lo)))
            if e.group == 'VOLUMICS' or (self.dimension == 3 and e.label == 'SHQ4'):
                self.volume_ids.add(e.index)
            elif (self.dimension == 3 and e.group in ('SHELLS', 'MEMBRANES', 'CONTACT')
                  and e.nen in (3, 4)):
                self.surface_ids.add(e.index)
        self.max_element_size = max((r.size for r in self.records.values()), default=0.)
        extent = self.bound_box[1]-self.bound_box[0]
        longest = float(np.max(extent))
        count = min(max_cells_per_axis, max(1, int(np.ceil(len(self.records)**(1/self.dimension)))))
        self.cell_size = longest/count if longest > 0 else 1.
        self.n_cells = np.maximum(1, np.ceil(extent/self.cell_size).astype(int))
        self.origin = self.bound_box[0].copy()
        for index, record in self.records.items():
            # nextafter avoids losing cells at exact boundaries through rounding.
            low, high = self._range(np.nextafter(record.lower, -np.inf),
                                    np.nextafter(record.upper, np.inf))
            if np.prod(high-low+1) > 4096:
                self.overflow.add(index)
            else:
                for key in product(*(range(low[i], high[i]+1) for i in range(self.dimension))):
                    self.cells.setdefault(key, set()).add(index)

    def _range(self, lower, upper):
        low = np.floor((lower-self.origin)/self.cell_size).astype(int)
        high = np.floor((upper-self.origin)/self.cell_size).astype(int)
        return np.clip(low, 0, self.n_cells-1), np.clip(high, 0, self.n_cells-1)

    def query_box(self, lower, upper):
        lower, upper = np.asarray(lower), np.asarray(upper)
        if np.any(upper < self.bound_box[0]) or np.any(lower > self.bound_box[1]):
            return []
        low, high = self._range(lower, upper)
        found = self.overflow.copy()
        for key in product(*(range(low[i], high[i]+1) for i in range(self.dimension))):
            found.update(self.cells.get(key, ()))
        return [self.records[i] for i in sorted(found)
                if np.all(self.records[i].upper >= lower) and np.all(self.records[i].lower <= upper)]

    def query_point(self, point, padding=0.):
        return self.query_box(point-padding, point+padding)

    def query_segment(self, a, b, padding=0.):
        direction = b-a
        for r in self.query_box(np.minimum(a, b)-padding, np.maximum(a, b)+padding):
            if box_interval(a, direction, r.lower-padding, r.upper+padding) is not None:
                yield r

    def near_surfaces(self, point, allowed=None):
        """Yield (cell distance lower bound, record) outward from the query.

        The caller may stop once this bound exceeds its best projection distance
        (plus tolerance). Elements are deduplicated across cells. Starting from
        the clamped point cell also handles queries outside the mesh box.
        """
        eligible = self.surface_ids if allowed is None else self.surface_ids.intersection(allowed)
        if not eligible:
            return
        seen = set()
        for index in sorted(self.overflow.intersection(eligible)):
            seen.add(index)
            yield 0., self.records[index]
        first = tuple(self._range(point, point)[0])
        visited = {first}
        def bound(key):
            lo = self.origin+np.array(key)*self.cell_size
            hi = np.minimum(lo+self.cell_size, self.bound_box[1])
            return float(np.linalg.norm(np.maximum(np.maximum(lo-point, point-hi), 0.)))
        heap = [(bound(first), first)]
        while heap:
            distance, key = heapq.heappop(heap)
            # Yield an empty record too: caller can stop before entering cells
            # containing no previously unseen elements.
            yield distance, None
            for index in sorted(self.cells.get(key, set()).intersection(eligible)-seen):
                seen.add(index)
                yield distance, self.records[index]
            if seen == eligible:
                return
            for axis in range(self.dimension):
                for step in (-1, 1):
                    neighbour = list(key)
                    neighbour[axis] += step
                    neighbour = tuple(neighbour)
                    if neighbour not in visited and 0 <= neighbour[axis] < self.n_cells[axis]:
                        visited.add(neighbour)
                        heapq.heappush(heap, (bound(neighbour), neighbour))


def get_mesh_grid(mesh):
    """Lazy compatibility for manually constructed Mesh/SimpleNamespace objects."""
    grid = getattr(mesh, 'spatial_grid', None)
    if grid is None:
        grid = MeshSpatialGrid(mesh)
        mesh.spatial_grid = grid
    return grid
