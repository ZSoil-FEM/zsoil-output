from .C_Material import *
from .C_Exf      import *
from .C_Element  import *



#=====================================================
class Element_heat_exch_L2 (Element):
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        Element.__init__ (self,my_mesh)
        self.group    = "HEAT_EXCH"
        self.rsl_ext  = ".s22"
        self.xsiGP    = [0.0]
        self.nen      = 4

    #=====================================================
    def get_xyz_gp (self,igaus):
    #=====================================================
        return self.get_center ()





