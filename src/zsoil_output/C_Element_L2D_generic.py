from .C_Element_Surf_generic import *
import numpy as np


# =====================================================
class Element_L2D_generic(Element_Surf_generic):
    # =====================================================

    # =====================================================
    def __init__(self):
        # =====================================================
        # Element_Surf_generic.__init__(self)
        self.nen = 2
        self.xsi_gp_1 = [0.0]
        self.W_gp_1 = [2.0]
        sqrt_3_by_3 = 1.7320508075688772 / 3.0
        self.xsi_gp_2 = [-sqrt_3_by_3, sqrt_3_by_3]
        self.W_gp_2 = [1.0, 1.0]

    # =====================================================
    def get_quadr(self, quadr_enum):
        # =====================================================
        if quadr_enum == self.QUADR_STD:
            return 2, self.xsi_gp_2, self.W_gp_2
        elif quadr_enum == self.QUADR_CENTRAL:
            return 1, self.xsi_gp_1, self.W_gp_1

    # =====================================================
    def get_xsi(self, quadr_enum, igaus):
        # =====================================================
        if quadr_enum == self.QUADR_CENTRAL:
            return self.xsi_gp_1
        elif quadr_enum == self.QUADR_STD:
            return self.xsi_gp_2[igaus]
        else:
            return None

    # =====================================================
    def get_N(self, xsi):
        # =====================================================
        s = xsi
        N = np.zeros(2)

        N[0] = (1.0 - s) / 2.0
        N[1] = (1.0 + s) / 2.0

        return N

    # =====================================================
    def get_N0(self):
        # =====================================================
        N = np.zeros(2)
        N[0] = 0.5
        N[1] = 0.5
        return N

    # =====================================================
    def get_DN(self, xsi):
        # =====================================================
        s = xsi
        DN = np.zeros(2).reshape(2,1)

        DN[0][0] = -0.5
        DN[1][0] =  0.5

        return DN

    # =====================================================
    def get_DN0(self):
        # =====================================================
        DN = np.zeros(2).reshape(2,1)

        DN[0][0] = -0.5
        DN[1][0] =  0.5

        return DN

    # =====================================================
    def get_DN_ex(self, quadr_enum, igaus):
        # =====================================================
        if quadr_enum == self.QUADR_CENTRAL:
            return self.get_DN0()
        else:
            xsi = self.get_xsi(quadr_enum, igaus)
            return self.get_DN(xsi)

    # =====================================================
    def get_N_ex(self, quadr_enum, igaus):
        # =====================================================
        if quadr_enum == self.QUADR_CENTRAL:
            return self.get_N0()
        else:
            xsi = self.get_xsi(quadr_enum, igaus)
            return self.get_N(xsi)

    # =====================================================
    def get_nr_of_faces(self):
        # =====================================================
        return 1

    # =====================================================
    def get_face_node_indices(self, face_index_1):
        # =====================================================
        return [1,2]


    #=====================================================
    def get_ndmG (self):
    #=====================================================
        return 3

    #=====================================================
    def get_ndmL (self):
    #=====================================================
        return 1


# =====================================================
def main():
    # =====================================================
    l2 = Element_L2D_generic()
    ele_coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    S = l2.get_surface(ele_coord, l2.QUADR_STD)
    thick_at_nodes = np.zeros(2)
    thick_at_nodes.fill(1.0)
    V = l2.get_volume (ele_coord, thick_at_nodes, l2.QUADR_STD)
    print(S, V)

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
