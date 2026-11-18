"""Linear, element-centre patch recovery for ZSoil continuum and shell results.

This module is part of ``zsoil-output``, a Python SDK developed for use
with ZSoil FEM to read model data and output-result files.

The :class:`SuperConvPatchRecovery` manager implements the linear patch
recovery corresponding to the ZSoil ``SPR_ex`` workflow. Indices in this
API are SDK indices starting at 1. Supplied values contain one scalar or
vector per element centre in a common component frame; Gauss-point
averaging and shell-tensor rotations remain the caller's responsibility.

Only active, same-material elements enter a target element/local-node
patch. Shell-normal compatibility and geometric rank tests control patch
membership and enlargement. Singular fits retain the constant term while
discarding weak coordinate terms.
"""

from collections import defaultdict
from collections.abc import Mapping

import numpy as np

from .C_Element import Element


class SuperConvPatchRecovery:
    """Geometry/patch manager for one element group at one physical time.

    Only same-material active elements enter a target-element/local-node patch.
    Shell normals must have signed dot product >= normal_cos_min with the target;
    face neighbours used in enlargement also obey that criterion pairwise.
    Small or geometrically rank-deficient patches get ONE neighbour ring.
    This includes corner nodes with just one active incident element. The
    extension is simultaneous from the original patch, never recursive.
    Singular fits cancel weak coordinate terms, retaining the constant term.
    Shell fits use two coordinates in a best-fit plane through patch nodes
    (or the target-element plane when shell_plane='element').
    """

    def __init__(self, mesh, ele_group_index, at_time, consider_existence=True,
                 mat_filter=None, inertia_tol=1.0e-2, rank_tol=1.0e-12,
                 shell_plane='patch', normal_cos_min=0.867):
        """Initialize the super conv patch recovery object for the supplied project data.

        :param mesh: Mesh containing the model nodes, elements, and materials.
        :type mesh: Mesh
        :param ele_group_index: Numeric element-group identifier.
        :type ele_group_index: int
        :param at_time: Physical time at which the model state is evaluated.
        :type at_time: float
        :param consider_existence: Consider existence used by this operation.
        :type consider_existence: bool
        :param mat_filter: Material indices to include; an empty sequence includes all materials.
        :type mat_filter: Optional[Sequence[int]]
        :param inertia_tol: Inertia tol used by this operation.
        :type inertia_tol: float
        :param rank_tol: Rank tol used by this operation.
        :type rank_tol: float
        :param shell_plane: Shell plane used by this operation.
        :type shell_plane: str
        :param normal_cos_min: Normal cos min used by this operation.
        :type normal_cos_min: float
        """
        if ele_group_index not in (Element.GROUP_CONTINUUM, Element.GROUP_SHELL):
            raise ValueError('SPR supports continuum and shell groups only')
        if not 0 <= inertia_tol < 1 or not 0 < rank_tol < 1:
            raise ValueError('Invalid inertia_tol or rank_tol')
        if shell_plane not in ('patch', 'element'):
            raise ValueError("shell_plane must be 'patch' or 'element'")
        if not -1.0 <= normal_cos_min <= 1.0:
            raise ValueError('normal_cos_min must be in [-1,1]')
        self.normal_cos_min = float(normal_cos_min)
        self.shell_plane = shell_plane
        self.mesh_ref = mesh
        self.ele_group_index = ele_group_index
        self.consider_existence = bool(consider_existence)
        self.mat_filter = tuple(() if mat_filter is None else mat_filter)
        self.inertia_tol = float(inertia_tol)
        self.rank_tol = float(rank_tol)
        self.update(at_time)

    def update(self, at_time):
        """Rebuild after changing physical time, connectivity or coordinates.

        :param at_time: Physical time at which the model state is evaluated.
        :type at_time: float
        :return: This recovery manager after rebuilding its geometry and patches.
        :rtype: SuperConvPatchRecovery
        """
        if self.consider_existence and at_time is None:
            raise ValueError('Physical at_time is required for existence filtering')
        # Node.adj_elements contains topology, never a time-filtered snapshot.
        # Do not clear existing lists: other postprocessors may own their use.
        for element in self.mesh_ref.elements:
            if element.get_group_index() == self.ele_group_index:
                for node_index in element.nodes:
                    node = self.mesh_ref.get_node(node_index)
                    if not any(e is element for e in node.adj_elements):
                        node.attach_element_to_adj_list(element)
        elements = {}
        for e in self.mesh_ref.elements:
            if e.get_group_index() != self.ele_group_index:
                continue
            if self.mat_filter:
                mat = self.mesh_ref.materials[e.material_index - 1]
                if mat.index_in_inp not in self.mat_filter:
                    continue
            if self.consider_existence:
                # SDK Element.is_ON assumes a non-null existence function.
                exf = e.give_exf()
                if exf is not None and not exf.is_ON(at_time):
                    continue
            elements[e.index] = e
        dimension = 2 if self.ele_group_index == Element.GROUP_SHELL else self.mesh_ref.ndm
        if dimension not in (2, 3):
            raise ValueError('Expected a 2D or 3D mesh')
        faces = defaultdict(list)
        centres = {}
        for index, e in elements.items():
            centres[index] = np.asarray(e.get_center(), dtype=float)
            if not np.all(np.isfinite(centres[index])):
                raise ValueError('Non-finite element coordinates')
            for face in self._faces(e, dimension):
                faces[tuple(sorted(face))].append(index)
        normals = {}
        if self.ele_group_index == Element.GROUP_SHELL:
            for index, e in elements.items():
                normal = np.asarray(e.ref_ele.get_normal(e.get_ele_coord()), dtype=float)
                length = np.linalg.norm(normal)
                if not np.isfinite(length) or length == 0:
                    raise ValueError('Degenerate shell normal at element %s' % index)
                normals[index] = normal / length

        def compatible(target, candidate):
            """Return whether a candidate element can belong to the current recovery patch.

            :param target: Target used by this operation.
            :type target: Element
            :param candidate: Candidate used by this operation.
            :type candidate: Element
            :return: Whether a candidate element can belong to the current recovery patch.
            :rtype: bool
            """
            if target.material_index != candidate.material_index:
                return False
            return (not normals or float(np.dot(normals[target.index], normals[candidate.index]))
                    >= self.normal_cos_min)

        neighbours = {index: set() for index in elements}
        for owners in faces.values():
            for index in owners:
                neighbours[index].update(i for i in owners if i != index
                                         and compatible(elements[index], elements[i]))
        patches = {}
        for e in elements.values():
            patches[e.index] = []
            for local, node_index in enumerate(e.nodes, 1):
                node = self.mesh_ref.get_node(node_index)
                original = {a.index for a in node.get_adj_elements(self.ele_group_index)
                            if a.index in elements and compatible(e, elements[a.index])}
                patch = set(original)
                sample_xyz = np.array([centres[i] for i in sorted(original)])
                if dimension == 2 and not normals:
                    sample_xyz = sample_xyz[:, :2]
                spread = sample_xyz - sample_xyz.mean(axis=0)
                singular = np.linalg.svd(spread, compute_uv=False)
                rank = int(np.count_nonzero(singular > self.rank_tol * singular[0])) if singular[0] > 0 else 0
                if len(original) < dimension + 1 or rank < dimension:
                    for index in original:
                        patch.update(i for i in neighbours[index]
                                     if compatible(e, elements[i]))
                patches[e.index].append(sorted(patch))
        # Commit element-owned lists only after the new geometry is built.
        for e in self.mesh_ref.elements:
            if e.get_group_index() == self.ele_group_index:
                e.spr_patches = patches.get(e.index, [[] for _ in range(e.nen)])
        self.at_time = at_time
        self.dimension = dimension
        self.elements = elements
        self.centres = centres
        self.normals = normals
        self.face_neighbours = {i: tuple(sorted(n)) for i, n in neighbours.items()}
        self.ele_patches = patches
        return self

    @staticmethod
    def _faces(element, dimension):
        """Return face-node groups used to connect neighboring elements.

        :param element: Element used by this operation.
        :type element: Element
        :param dimension: Dimension used by this operation.
        :type dimension: int
        :return: The face-node groups used to connect neighboring elements.
        :rtype: list[tuple[int, ...]]
        """
        nodes = element.nodes
        if dimension == 2:
            if len(nodes) not in (3, 4):
                raise NotImplementedError('SPR supports T3/Q4 surface topology')
            return [(nodes[i], nodes[(i + 1) % len(nodes)]) for i in range(len(nodes))]
        if element.ref_ele is not None:
            ref = element.ref_ele
            return [tuple(nodes[j - 1] for j in ref.get_face_info(i + 1)[1])
                    for i in range(ref.get_nr_of_faces())]
        if len(nodes) == 4:  # TH4
            topology = ((0, 1, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))
        elif len(nodes) == 6:  # W6: bottom/top triangles, three quadrilateral faces
            topology = ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3),
                        (1, 2, 5, 4), (2, 0, 3, 5))
        else:
            raise NotImplementedError('Unsupported 3D element topology')
        return [tuple(nodes[j] for j in face) for face in topology]

    def get_patch(self, ele_index_1, local_node_index_1):
        """Sorted element IDs for an element-local node (both indices 1-based).

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param local_node_index_1: 1-based local node index.
        :type local_node_index_1: int
        :return: The sorted element IDs for an element-local node (both indices 1-based).
        :rtype: tuple[int, ...]
        """
        patches = self.ele_patches[ele_index_1]
        if not 1 <= local_node_index_1 <= len(patches):
            raise IndexError('Local node index must be in 1..nen')
        return tuple(patches[local_node_index_1 - 1])

    def _axes(self, element, selected):
        """Return coordinate axes for fitting values over the recovery patch.

        :param element: Element used by this operation.
        :type element: Element
        :param selected: Selected used by this operation.
        :type selected: Sequence[int]
        :return: The coordinate axes for fitting values over the recovery patch.
        :rtype: numpy.ndarray
        """
        if self.ele_group_index != Element.GROUP_SHELL:
            return np.eye(3)[:, :self.dimension]
        if self.shell_plane == 'patch':
            node_ids = sorted({n for i in selected for n in self.elements[i].nodes})
            xyz = np.array([self.mesh_ref.get_node(n).xyz for n in node_ids], dtype=float)
            _, singular, directions = np.linalg.svd(xyz - xyz.mean(axis=0), full_matrices=False)
            if len(singular) < 2 or singular[1] <= self.rank_tol * singular[0]:
                raise ValueError('Shell patch nodes do not define a plane')
            normal = directions[-1]
            if np.dot(normal, self.normals[element.index]) < 0:
                normal = -normal
            tangent = np.asarray(element.ref_ele.get_T_LG0(element.get_ele_coord()))[:, 0]
            tangent = tangent - np.dot(tangent, normal) * normal
            length = np.linalg.norm(tangent)
            if length <= self.rank_tol:
                tangent = directions[0]
            else:
                tangent = tangent / length
            return np.column_stack((tangent, np.cross(normal, tangent)))
        # Reference surface basis columns are tangent directions in global xyz.
        ref = element.ref_ele
        if ref is None:
            raise NotImplementedError('Shell requires a surface reference element')
        axes = np.asarray(ref.get_T_LG0(element.get_ele_coord()))[:, :2]
        if not np.all(np.isfinite(axes)):
            raise ValueError('Degenerate shell element')
        return axes

    def _fit(self, offsets, values):
        """Return value at origin and diagnostics; solve LS without normal equations.

        :param offsets: Offsets used by this operation.
        :type offsets: numpy.ndarray
        :param values: Values used by this operation.
        :type values: Union[Mapping, array_like]
        :return: The value at origin and diagnostics; solve LS without normal equations.
        :rtype: tuple[numpy.ndarray, dict]
        """
        centred = offsets - offsets.mean(axis=0)
        inertia = np.sum(centred ** 2, axis=0)
        largest = float(inertia.max())
        keep = list(np.flatnonzero(inertia > self.inertia_tol * largest)) if largest > 0 else []
        # Cancel smallest inertia directions until the sample count is sufficient.
        while len(keep) + 1 > len(values):
            keep.remove(min(keep, key=lambda i: (inertia[i], i)))
        while True:
            if keep:
                # One isotropic scale keeps directional inertia ratios unchanged.
                scale = float(np.sqrt(largest / len(values)))
                design = np.column_stack((np.ones(len(values)), offsets[:, keep] / scale))
            else:
                design = np.ones((len(values), 1))
            coefficients, _, rank, singular = np.linalg.lstsq(design, values, rcond=self.rank_tol)
            if rank == design.shape[1]:
                return coefficients[0], {'rank': int(rank), 'kept_directions': tuple(keep),
                                         'inertia': inertia, 'singular_values': singular,
                                         'sample_count': len(values)}
            if not keep:
                raise np.linalg.LinAlgError('Constant patch fit is singular')
            keep.remove(min(keep, key=lambda i: (inertia[i], i)))

    @staticmethod
    def _values(values):
        """Validate and normalize the supplied element-center values.

        :param values: Values used by this operation.
        :type values: Union[Mapping, array_like]
        :return: Element-indexed mapping of normalized NumPy value arrays.
        :rtype: dict[int, numpy.ndarray]
        """
        if not isinstance(values, Mapping):
            raise TypeError('values must map 1-based element IDs to centre scalar/vectors')
        result = {}
        size = None
        for index, value in values.items():
            array = np.asarray(value, dtype=float)
            if array.ndim > 1:
                raise ValueError('Supply one centre scalar/vector per element, not Gauss arrays')
            array = np.atleast_1d(array)
            if not array.size or not np.all(np.isfinite(array)):
                raise ValueError('Element values must be non-empty and finite')
            if size is None:
                size = array.size
            if array.size != size:
                raise ValueError('Inconsistent component counts')
            result[index] = array
        return result

    def recover_element_nodes(self, ele_index_1, values, valid_elements=None, return_diagnostics=False):
        """Return (nen,ncomps) recovered values in element connectivity order.

                valid_elements is an optional iterable of usable element IDs (the
                Fortran ele_ok_in_patch counterpart). Missing values without an explicit
                validity filter are errors, preventing accidental incomplete fits.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param values: Values used by this operation.
        :type values: Union[Mapping, array_like]
        :param valid_elements: Optional element indices allowed to contribute values.
        :type valid_elements: Optional[Collection[int]]
        :param return_diagnostics: Return diagnostics used by this operation.
        :type return_diagnostics: bool
        :return: The (nen,ncomps) recovered values in element connectivity order.
        :rtype: Union[numpy.ndarray, tuple[numpy.ndarray, list[dict]]]
        """
        data = self._values(values)
        return self._recover_nodes(ele_index_1, data, valid_elements, return_diagnostics)

    def recover_node(self, ele_index_1, local_node_index_1, values, valid_elements=None, return_diagnostics=False):
        """Recover just one element-local node; return a vector of length ncomps.

                Both indices start at 1. Only values for this node's patch are needed.
                For a scalar result, use returned_vector[0].

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param local_node_index_1: 1-based local node index.
        :type local_node_index_1: int
        :param values: Values used by this operation.
        :type values: Union[Mapping, array_like]
        :param valid_elements: Optional element indices allowed to contribute values.
        :type valid_elements: Optional[Collection[int]]
        :param return_diagnostics: Return diagnostics used by this operation.
        :type return_diagnostics: bool
        :return: The recover just one element-local node; return a vector of length ncomps.
        :rtype: Union[numpy.ndarray, tuple[numpy.ndarray, dict]]
        """
        data = self._values(values)
        allowed = None if valid_elements is None else set(valid_elements)
        recovered, info = self._recover_node_data(ele_index_1, local_node_index_1,
                                                  data, allowed)
        return (recovered, info) if return_diagnostics else recovered

    def _recover_node_data(self, ele_index_1, local_node_index_1, data, allowed):
        """Recover and return node data.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param local_node_index_1: 1-based local node index.
        :type local_node_index_1: int
        :param data: Data used by this operation.
        :type data: list
        :param allowed: Allowed used by this operation.
        :type allowed: Optional[Collection[int]]
        :return: The node data.
        :rtype: tuple[numpy.ndarray, dict]
        """
        element = self.elements[ele_index_1]
        patch = self.get_patch(ele_index_1, local_node_index_1)
        selected = [i for i in patch if allowed is None or i in allowed]
        if not selected:
            raise ValueError('No usable samples in patch at element %s, node %s' %
                             (ele_index_1, local_node_index_1))
        missing = [i for i in selected if i not in data]
        if missing:
            raise ValueError('Missing centre values for elements %s' % missing)
        axes = self._axes(element, selected)
        node_index = element.nodes[local_node_index_1 - 1]
        origin = np.asarray(self.mesh_ref.get_node(node_index).xyz, dtype=float)
        offsets = (np.array([self.centres[i] for i in selected]) - origin) @ axes
        recovered, info = self._fit(offsets, np.array([data[i] for i in selected]))
        info['element_indices'] = tuple(selected)
        return recovered, info

    def _recover_nodes(self, ele_index_1, data, valid_elements, diagnostics):
        """Recover and return nodes.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param data: Data used by this operation.
        :type data: list
        :param valid_elements: Optional element indices allowed to contribute values.
        :type valid_elements: Optional[Collection[int]]
        :param diagnostics: Diagnostics used by this operation.
        :type diagnostics: bool
        :return: The nodes.
        :rtype: Union[numpy.ndarray, tuple[numpy.ndarray, list[dict]]]
        """
        element = self.elements[ele_index_1]
        allowed = None if valid_elements is None else set(valid_elements)
        output, info = [], []
        for local in range(1, element.nen + 1):
            recovered, fit_info = self._recover_node_data(ele_index_1, local, data, allowed)
            output.append(recovered)
            info.append(fit_info)
        result = np.array(output)
        return (result, info) if diagnostics else result

    def recover(self, values, valid_elements=None):
        """Map active element ID to (nen,ncomps) recovered values.

        :param values: Values used by this operation.
        :type values: Union[Mapping, array_like]
        :param valid_elements: Optional element indices allowed to contribute values.
        :type valid_elements: Optional[Collection[int]]
        :return: Recovered nodal values keyed by active element index.
        :rtype: dict[int, numpy.ndarray]
        """
        data = self._values(values)
        allowed = None if valid_elements is None else set(valid_elements)
        return {i: self._recover_nodes(i, data, allowed, False) for i in self.elements}

    def set_values_at_pt(self, ele_index_1, xi, values, valid_elements=None):
        """Shape-function interpolation of nodal recovered values at natural xi.

        :param ele_index_1: 1-based element index.
        :type ele_index_1: int
        :param xi: Xi used by this operation.
        :type xi: Sequence[float]
        :param values: Values used by this operation.
        :type values: Union[Mapping, array_like]
        :param valid_elements: Optional element indices allowed to contribute values.
        :type valid_elements: Optional[Collection[int]]
        :return: Recovered value interpolated at the requested natural coordinates.
        :rtype: numpy.ndarray
        """
        element = self.elements[ele_index_1]
        nodal = self.recover_element_nodes(ele_index_1, values, valid_elements)
        if element.ref_ele is not None:
            shape = element.ref_ele.get_N(xi)
        elif self.dimension == 3 and element.nen == 4:
            r, s, t = xi
            shape = np.array([1 - r - s - t, r, s, t])
        elif self.dimension == 3 and element.nen == 6:
            r, s, t = xi
            shape = np.array([(1-r-s)*(1-t), r*(1-t), s*(1-t),
                              (1-r-s)*(1+t), r*(1+t), s*(1+t)]) * 0.5
        else:
            raise NotImplementedError('Element has no supported shape functions')
        shape = np.asarray(shape, dtype=float)
        if shape.shape != (element.nen,):
            raise ValueError('Incompatible element shape-function vector')
        return shape @ nodal
