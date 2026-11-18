from .C_Element_Surf_generic import *
from .C_Element_Volumic_generic import *
from .C_Element_Q4_generic import *
import numpy as np

#=====================================================
class Element_Q4V_generic (Element_Volumic_generic):
#=====================================================

    #=====================================================
    def __init__ (self):
    #=====================================================
        self.Q4S = Element_Q4_generic ()

    #=====================================================
    def get_quadr (self,quadr_enum):
    #=====================================================
        return self.Q4S.get_quadr (quadr_enum)

    #=====================================================
    def get_xsi (self,quadr_enum,igaus):
    #=====================================================
        return self.Q4S.get_xsi (quadr_enum,igaus)

    #=====================================================
    def get_N (self,xsi):
    #=====================================================
        return self.Q4S.get_N (xsi)

    #=====================================================
    def get_N0 (self):
    #=====================================================
        return self.Q4S.get_N0()

    #=====================================================
    def get_DN (self,xsi):
    #=====================================================
        return self.Q4S.get_DN (xsi)

    #=====================================================
    def get_DN0 (self):
    #=====================================================
        return self.Q4S.get_DN0 ()

    #=====================================================
    def get_DN_ex (self,quadr_enum,igaus):
    #=====================================================
        return self.Q4S.get_DN_ex (quadr_enum,igaus)

    #=====================================================
    def get_N_ex (self,quadr_enum,igaus):
    #=====================================================
        return self.Q4S.get_N_ex(quadr_enum, igaus)

    #=====================================================
    def get_nr_of_faces (self):
    #=====================================================
        return 4


    #=====================================================
    def get_face_node_indices(self,face_index_1):
    #=====================================================
        n1 = face_index_1
        n2 = n1 + 1
        if n2 > 4:
            n2 = 1
        return [n1,n2]

    # =============================================
    def get_face_info (self,face_index_1):
    # =============================================
        return 'L2',self.get_face_node_indices (face_index_1)


#=====================================================
def main():
#=====================================================
    q4 = Element_Q4_generic ()
    ele_coord = np.array ([[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]])
    thick_at_nodes = np.array ([1.0,1.0,2.0,2.0])

    S = q4.get_surface(ele_coord,q4.QUADR_STD)
    V = q4.get_volume (ele_coord,thick_at_nodes,q4.QUADR_STD)
    print(S,V)

    # #print ele_coord
    # T_LG = q4.get_T_LG0 (ele_coord)
    # print T_LG
    # T_GL = q4.get_T_GL0 (ele_coord)
    # print T_GL
    #
    # user_vec = np.array ([1.0,1.0,0.0])
    # T_VEC = q4.get_T_user_VEC (T_GL,user_vec,'T')
    # print T_VEC
    # print T_VEC.dot (user_vec)
    # print q4.get_T_user_TNS (T_VEC)


if __name__ == '__main__':
    main()