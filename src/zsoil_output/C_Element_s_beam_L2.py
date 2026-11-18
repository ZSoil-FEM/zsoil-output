from .C_Material import *
from .C_Exf      import *
from .C_Element  import *



#=====================================================
class Element_s_beam_L2 (Element):
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        Element.__init__ (self,my_mesh)
        self.group    = "BEAMS"
        self.rsl_ext  = ".s04"
        # self.centroids = self.nodes [2:3]
        # self.directors = self.nodes [4:5]
        # del self.nodes [2:]
        self.xsiGP    = [[-1.0],[-0.654653670707977],[0.0],[0.654653670707977],[1.0]]
        self.nen      = 2

    #=====================================================
    def instanciate (self,my_mesh):
    #=====================================================
        Element.instanciate (self,my_mesh)

        self.centroids = self.nodes [2:3]
        self.directors = self.nodes [4:5]
        del self.nodes [2:]


    #=====================================================
    def get_xyz_gp (self,igaus):
    #=====================================================
        node1 = self.mesh_ref.nodes [self.nodes [0]-1]
        node2 = self.mesh_ref.nodes [self.nodes [1]-1]

        xyz1 = node1.get_xyz()
        xyz2 = node2.get_xyz()

        xyz = [0.0,0.0,0.0]
        xsi = self.xsiGP [igaus][0]
        N1  = 0.5*(1.0-xsi)
        N2  = 0.5*(1.0+xsi)
        for i in range (len(xyz1)):
            xyz [i] = xyz [i] + N1 * xyz1 [i] + N2 * xyz2 [i]
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


