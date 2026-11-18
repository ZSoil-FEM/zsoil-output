from .C_Material import *
from .C_Exf import *
from .C_Element import *


# =====================================================
class Element_cont_TH4(Element):
    # =====================================================

    # =====================================================
    def __init__(self, my_mesh):
        # =====================================================
        Element.__init__(self, my_mesh)
        self.group = "VOLUMICS"
        self.rsl_ext = ".s01"
        self.xsiGP = [[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]]
        self.nen = 4

    # =====================================================
    def get_xyz_gp(self, igaus):
        # =====================================================
        xyz = [0.0, 0.0, 0.0]
        for i in range(self.nen):
            node = self.mesh_ref.nodes[self.nodes[i] - 1]
            for j in range(len(node.xyz)):
                xyz[j] = xyz[j] + node.xyz[j] / self.nen

        return xyz
