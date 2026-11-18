from .C_Material import *
from .C_Exf      import *
from .C_Element  import *



#=====================================================
class Element_cnt_L2 (Element):
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        Element.__init__ (self,my_mesh)
        self.group    = "CONTACT"
        self.rsl_ext  = ".s07"
        self.xsiGP    = [[-1.0],[1.0]]
        self.nen      = 2

    #=====================================================
    def get_xyz_gp (self,igaus):
    #=====================================================
        node = self.mesh_ref.nodes [self.nodes [igaus]]

        xyz = [0.0,0.0,0.0]
        for i in range (len(self.xsiGP)):
            xyz [i] = node.xyz [i]

        return xyz



