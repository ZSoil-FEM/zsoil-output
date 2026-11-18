from .C_Material import *
from .C_Exf      import *
from .C_Element  import *



#=====================================================
class Element_cnt_Q4 (Element):
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        Element.__init__ (self,my_mesh)
        self.group    = "CONTACT"
        self.rsl_ext  = ".s07"
        self.xsiGP    = [[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0]]
        self.nen      = 4
        self.ref_ele = Element_Q4_generic()
        self.user_vec = None
        self.T_TNS = None
        self.T_VEC = None

    #=====================================================
    def get_xyz_gp (self,igaus):
    #=====================================================
        node = self.mesh_ref.nodes [self.nodes [igaus]]

        xyz = [0.0,0.0,0.0]
        for i in range (len(self.xsiGP)):
            xyz [i] = node.xyz [i]

        return xyz

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



