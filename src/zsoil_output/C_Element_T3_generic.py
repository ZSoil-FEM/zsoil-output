from .C_Element_Surf_generic import *
import numpy as np

#=====================================================
class Element_T3_generic (Element_Surf_generic):
#=====================================================

    #=====================================================
    def __init__ (self,):
    #=====================================================
        self.nen = 3
        self.xsi_gp_1 = [1.0/3.0,1.0/3.0]
        self.W_gp_1   = [0.5]

    #=====================================================
    def get_quadr (self,quadr_enum):
    #=====================================================
        if quadr_enum == self.QUADR_STD:
            return 1, self.xsi_gp_1, self.W_gp_1
        elif quadr_enum == self.QUADR_CENTRAL:
            return 1, self.xsi_gp_1, self.W_gp_1

    #=====================================================
    def get_xsi (self,quadr_enum,igaus):
    #=====================================================
        if quadr_enum == self.QUADR_STD:
            return self.xsi_gp_1
        else:
            return self.xsi_gp_1

    #=====================================================
    def get_N (self,xsi):
    #=====================================================
        s = xsi [0]
        t = xsi [1]

        N = np.zeros (3)

        N [0] = 1.0-s-t
        N [1] = s
        N [2] = t

        return N

    #=====================================================
    def get_N0 (self):
    #=====================================================
        N = np.zeros (3)
        N [0] = 1.0/3.0
        N [1] = 1.0/3.0
        N [2] = 1.0/3.0
        return N

    #=====================================================
    def get_DN (self,xsi):
    #=====================================================
        s  = xsi [0]
        t  = xsi [1]

        DN = np.zeros (6).reshape (3,2)

        DN [0,0] = -1.0
        DN [0,1] = -1.0

        DN [1,0] =  1.0
        DN [1,1] =  0.0

        DN [2,0] =  0.0
        DN [2,1] =  1.0

        return DN

    #=====================================================
    def get_DN0 (self):
    #=====================================================
        return self.get_DN (self.xsi_gp_1)


    # =====================================================
    def get_DN(self, quadr_enum, igaus):
        # =====================================================
        if quadr_enum == self.QUADR_CENTRAL:
            return self.get_DN0()
        else:
            xsi = self.get_xsi(quadr_enum, igaus)
            return self.get_DN(xsi)


    # =====================================================
    def get_N(self, quadr_enum, igaus):
        # =====================================================
        if quadr_enum == self.QUADR_CENTRAL:
            return self.get_N0()
        else:
            xsi = self.get_xsi(quadr_enum, igaus)
            return self.get_N(xsi)



#=====================================================
def main():
#=====================================================
    t3 = Element_T3_generic ()
    ele_coord = array ([[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0]])
    T_GL = t3.get_T_GL0 (ele_coord)
    user_vec = array ([1.0,1.0,0.0])
    T_VEC = t3.get_T_user_VEC (T_GL,user_vec,'T')
    print(T_VEC)
    if T_VEC != None:
        print(T_VEC.dot (user_vec))
    print(t3.get_T_user_TNS (T_VEC))


if __name__ == '__main__':
    main()