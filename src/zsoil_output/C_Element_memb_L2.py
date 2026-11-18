from .C_Material import *
from .C_Exf      import *
from .C_Element  import *



#=====================================================
class Element_memb_L2 (Element):
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        Element.__init__ (self,my_mesh)
        self.group    = "MEMBRANE"
        self.rsl_ext  = ".s15"
        self.xsiGP    = [[0.0]]
        self.nen      = 2


    #=====================================================
    def get_xyz_gp (self,igaus):
    #=====================================================
        xyz = [0.0,0.0,0.0]
        for i in range (self.nen):
            node = self.mesh_ref.nodes [self.nodes [i]-1]
            for j in range (len(node.xyz)):
                xyz [j] = xyz [j] + node.xyz [j] / self.nen

        return xyz



