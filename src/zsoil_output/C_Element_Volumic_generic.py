import numpy as np

#=====================================================
class Element_Volumic_generic ():
#=====================================================

    QUADR_CENTRAL = 0
    QUADR_STD     = 1
    QUADR_NODAL   = 2

    #=====================================================
    def __init__ (self,):
    #=====================================================
        self.nen = 0

    #=====================================================
    def get_quadr (self,quadr_enum):
    #=====================================================
        pass

    #=====================================================
    def get_xsi (self,ngaus,igaus):
    #=====================================================
        pass

    #=====================================================
    def get_N_ex (self,quadr_enum,igaus):
    #=====================================================
        pass

    #=====================================================
    def get_N (self,xsi):
    # =====================================================
        pass

    #=====================================================
    def get_DN0 (self):
    # =====================================================
        pass

    #=====================================================
    def get_DN_ex (self,quadr_enum,igaus):
    #=====================================================
        pass

    #=====================================================
    def get_volume (self,ele_coord,quadr_enum,ret_V_igaus=False):
    #=====================================================
        ngaus, xsi, Wi = self.get_quadr (quadr_enum)
        Vi  = np.zeros (len(xsi))

        V = 0.0

        ndmG = 3
        ndmL = 3

        dxG_dxL = np.zeros(ndmG * ndmL).reshape(ndmG,ndmL)

        for igaus in range (ngaus):
            DN = self.get_DN_ex (quadr_enum,igaus)
            # dx/dxsi
            dxG_dxL.fill (0.0)
            for inode in range (self.nen):
                for idmL in range (ndmL):
                    for idmG in range (ndmG):
                        dxG_dxL [idmG,idmL] = dxG_dxL [idmG,idmL] + \
                                              DN [inode,idmL] * ele_coord [inode,idmG]
            aux = norm = np.linalg.det (dxG_dxL)
            V   =  V + aux * Wi [igaus]
            if ret_V_igaus:
                Vi [igaus] = aux * Wi [igaus]

        if ret_V_igaus:
            return Vi
        else:
            return V


