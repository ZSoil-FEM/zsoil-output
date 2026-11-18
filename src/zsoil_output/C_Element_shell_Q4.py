from .C_Material import *
from .C_Exf      import *
from .C_Element            import *
from .C_Element_Q4_generic import *
import numpy as np

#=====================================================
class Element_shell_Q4 (Element):
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        Element.__init__ (self,my_mesh)
        self.group    = "SHELLS"
        self.rsl_ext  = ".s02"
        self.rsl_lay_ext = ".l02"
        self.xsiGP    = [[0.0,0.0]]
        self.nen      = 4
        self.ref_ele  = Element_Q4_generic ()
        self.user_vec = None
        self.T_TNS    = None
        self.T_VEC    = None

    #=====================================================
    def get_xyz_gp (self,igaus):
    #=====================================================
        xyz = [0.0,0.0,0.0]
        for i in range (self.nen):
            node = self.mesh_ref.nodes [self.nodes [i]-1]
            for j in range (len(node.xyz)):
                xyz [j] = xyz [j] + node.xyz [j] / self.nen

        return xyz

    #=====================================================
    def get_nlayers (self):
    #=====================================================
        mat_index = self.material_index
        material  = self.mesh_ref.materials [mat_index-1]
        main_data = material.data [material.dict['MAIN']]
        nlayers_all    = int(main_data [0])
        nlayers_reinf  = int(main_data [3])
        return nlayers_all,nlayers_reinf

    #=====================================================
    def set_transf_matrices (self,user_vec):
    #=====================================================
        self.user_vec = user_vec
        self.T_VEC    = None
        self.T_TNS    = None

        ele_coord = self.get_ele_coord ()

        self.T_LG0 = self.ref_ele.get_T_LG0 (ele_coord)
        self.T_GL0 = self.ref_ele.get_T_GL0 (ele_coord)

        self.T_VEC = self.ref_ele.get_T_user_VEC(self.T_GL0, user_vec, axis='T')
        if isinstance (self.T_VEC,np.ndarray):
            self.T_TNS = self.ref_ele.get_T_user_TNS(self.T_VEC)

    # =====================================================
    def get_T_TNS (self):
    # =====================================================
        return self.T_TNS

    # =====================================================
    def get_T_VEC (self):
    # =====================================================
        return self.T_VEC

    # =====================================================
    def get_T_GL (self):
    # =====================================================
        return self.T_GL0

    # =====================================================
    def get_T_LG (self):
    # =====================================================
        return self.T_LG0

    # =====================================================
    def get_surface (self,quadr_enum):
    # =====================================================
        ele_coord = self.get_ele_coord()
        return self.ref_ele.get_surface (ele_coord, quadr_enum)



