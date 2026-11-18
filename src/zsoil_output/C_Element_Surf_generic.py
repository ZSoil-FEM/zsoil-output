import numpy as np

#=====================================================
class Element_Surf_generic ():
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
    def get_T_LG_private (self,ele_coord,DN):
    #=====================================================
        SQRT_2 = 0.707106781186547524401

        ndmL = self.get_ndmL()
        ndmG = self.get_ndmG()

        dx_dxsi = np.zeros (ndmL*ndmG)
        dx_dxsi = dx_dxsi.reshape (ndmG,ndmL)

        n = 0
        for j in range (ndmL):#NdmL
            for i in range (ndmG):#NdmG
                dx_dxsi [i,j] = 0.0
                for k in range (self.nen):#Nen
                    dx_dxsi [i,j] = dx_dxsi [i,j] + DN [k,j] * ele_coord [k,i]

        if ndmL == 2:
            tmp = np.cross (dx_dxsi [:,0],dx_dxsi [:,1])
            exsi = dx_dxsi [:,0]
            eeta = dx_dxsi [:,1]
        else:
            ez = np.zeros(3)
            ez [2] = 1.0
            tmp = np.cross (dx_dxsi [:,0],ez)
            exsi = dx_dxsi [:,0]
            eeta = ez [:]


        norm = np.linalg.norm (tmp)
        tmp = tmp / norm

        # norm = np.linalg.norm (dx_dxsi [:,0])
        # dx_dxsi [:,0] = dx_dxsi [:,0] / norm
        #
        # norm = np.linalg.norm (dx_dxsi [:,1])
        # dx_dxsi [:,1] = dx_dxsi [:,1] / norm

        norm = np.linalg.norm (exsi)
        exsi [:] = exsi [:] / norm

        norm = np.linalg.norm(eeta)
        eeta [:] = eeta [:] / norm

        T_LG = np.zeros (9)
        T_LG = T_LG.reshape (3,3)

        T_LG [0,2] = tmp [0]
        T_LG [1,2] = tmp [1]
        T_LG [2,2] = tmp [2]

        ea = np.zeros (3)
        for j in range (3):
            # ea [j] = 0.5 * (dx_dxsi [j,0]+dx_dxsi [j,1])
            ea [j] = 0.5 * (exsi [j]+eeta [j])
        norm = np.linalg.norm (ea)
        ea = ea / norm

        eb = np.cross (tmp,ea)
        norm = np.linalg.norm (eb)
        eb = eb / norm

        for j in range (3):
            T_LG [j,0] = SQRT_2 * (ea [j] - eb [j])
            T_LG [j,1] = SQRT_2 * (ea [j] + eb [j])

        return T_LG


    #=====================================================
    def get_T_LG (self,ele_coord,xsi):
    #=====================================================
        DN = self.get_DN (xsi)
        return self.get_T_LG_private (ele_coord,DN)

    #=====================================================
    def get_T_GL (self,ele_coord,xsi):
    #=====================================================
        DN = self.get_DN (xsi)
        T_LG = self.get_T_LG_private (ele_coord,DN)
        T_GL = T_LG.T ()
        return T_GL

    #=====================================================
    def get_T_LG0 (self,ele_coord):
    #=====================================================
        DN = self.get_DN0 ()
        return self.get_T_LG_private (ele_coord,DN)

    #=====================================================
    def get_T_GL0 (self,ele_coord):
    #=====================================================
        DN = self.get_DN0 ()
        T_LG = self.get_T_LG_private (ele_coord,DN)
        return T_LG.T

    #=====================================================
    def get_normal (self,ele_coord):
    #=====================================================
        T_GL0 = self.get_T_GL0 (ele_coord)
        n = np.array (T_GL0 [2,:])
        return n

    #=====================================================
    def get_T_user_VEC (self,T_GL,user_vec,axis='T'):
    #=====================================================
        #it returns the transformation matrix in tensorial sense from local g.p. base
        #to another local base defined with the aid of the "user vector"
        #so   Muser = T * M g.p. * T'
        #     Nuser = T * N g.p. * T'
        #     Quser = T * Q g.p.
        # axis can be 'T' or 'N'

        #find user vector in local shell cordinate system
        norm = np.linalg.norm (user_vec)
        user_vecx = user_vec / norm

        ##print user_vecx

        tmp = T_GL.dot (user_vecx)  # matrix times vector multiplication in numpy

        if tmp[:2].dot (tmp [0:2]) < 1.0e-6:
            return None # vector perpendicular to the surface

        e   = np.zeros (3) #numpy vector
        #  define e1'  axis for 'T' or e2' axis for 'N'
        e [0] = tmp [0]
        e [1] = tmp [1]
        e [2] = 0.0
        norm = np.linalg.norm (e)
        e =  e / norm

        T = np.zeros (9)
        T = T.reshape (3,3)

        if ( axis == 'T' ):
            row = 0
        else:
            row = 1

        for i in range (3):
            T [row,i] = e [i]
        T [2,2] = 1.0

        if ( axis == 'T' ):
            row = 1
        else:
            row = 0

        T [row,0] = -e [1]
        T [row,1] =  e [0]

        return T

    #=====================================================
    def get_T_user_TNS (self,T_VEC):
    #=====================================================
        # this transformation matrix can be used as follows  M' = T * M_gp
        T = np.zeros (9).reshape (3,3)
        T [0,0] = T_VEC [0,0] * T_VEC [0,0]
        T [0,1] = T_VEC [0,1] * T_VEC [0,1]
        T [0,2] = T_VEC [0,0] * T_VEC [0,1] * 2.0

        T [1,0] = T_VEC [1,0] * T_VEC [1,0]
        T [1,1] = T_VEC [1,1] * T_VEC [1,1]
        T [1,2] = T_VEC [1,0] * T_VEC [1,1] * 2.0

        T [2,0] = T_VEC [0,0] * T_VEC [1,0]
        T [2,1] = T_VEC [0,1] * T_VEC [1,1]
        T [2,2] = T_VEC [0,0] * T_VEC [1,1] + T_VEC [0,1] * T_VEC [1,0]

        return T


    #=====================================================
    def get_surface (self,ele_coord,quadr_enum,ret_S_igaus=False):
    #=====================================================
        ngaus, xsi, Wi = self.get_quadr (quadr_enum)
        Si  = np.zeros (len(xsi))

        S = 0.0

        ndmG = self.get_ndmG()
        ndmL = self.get_ndmL()

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
            dxG_dxsi = dxG_dxL [:,0]
            if ndmL == 2:
                dxG_deta = dxG_dxL [:,1]
            else:
                dxG_deta = np.zeros (ndmG)
                dxG_deta [2] = 1.0 # ez

            #make cross product
            dx_dxsi_X_dx_deta = np.cross (dxG_dxsi,dxG_deta)
            aux = norm = np.linalg.norm (dx_dxsi_X_dx_deta) #sqrt(np.dot (dx_dxsi_X_dx_deta,dx_dxsi_X_dx_deta))
            S   =  S + aux * Wi [igaus]
            if ret_S_igaus:
                Si [igaus] = aux * Wi [igaus]

        if ret_S_igaus:
            return Si
        else:
            return S


    #=====================================================
    def get_volume  (self,ele_coord,thick_at_nodes,quadr_enum):
    #=====================================================
        Si = self.get_surface (ele_coord,quadr_enum,True)
        ngaus, xsi, Wi = self.get_quadr(quadr_enum)

        V = 0.0
        for igaus in range (ngaus):
            N = self.get_N_ex (quadr_enum,igaus)
            h_gp = np.dot (N,thick_at_nodes)
            V = V + Si [igaus] * h_gp
        return V

    #=====================================================
    def get_ndmG (self):
    #=====================================================
        return 3

    #=====================================================
    def get_ndmL (self):
    #=====================================================
        return 2

