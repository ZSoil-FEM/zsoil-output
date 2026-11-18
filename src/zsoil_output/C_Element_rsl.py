"""Map element groups to their ZSoil result-file extensions.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Element_rsl` mapping associates each supported element family with the binary
file that stores its integration-point results.
"""

class Element_rsl ():
    def __init__(self, my_mesh):
        """Initialize the element results object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        self.mesh_ref = my_mesh
        # file extensions for results in elements
        self.dict = {'VOLUMICS': '.s01', 'SHELLS': '.s02', 'TRUSSES': '.s03', 'BEAMS': '.s04',
                     'CONTACT': '.s07', 'MEMBRANE': '.s15', 'NS-CONTACT': '.s17', 'HEAT_EXCH': '.s22'}
        self.dict_inv = {v: k for k, v in list(self.dict.items())}
