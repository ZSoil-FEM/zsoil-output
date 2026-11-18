#=====================================================
class Element_rsl ():
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        self.mesh_ref = my_mesh
        #file extensions for results in elements
        self.dict     = {'VOLUMICS':'.s01','SHELLS':'.s02','TRUSSES':'.s03','BEAMS':'.s04','CONTACT':'.s07','MEMBRANE':'.s15','NS-CONTACT':'.s17','HEAT_EXCH':'.s22'}
        self.dict_inv = {v: k for k, v in list(self.dict.items())}
