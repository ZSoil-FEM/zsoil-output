#=====================================================
class Node ():
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        self.mesh_ref = my_mesh
        self.index = -1
        self.load_records = []
        self.xyz   = []
        self.adj_elements = []

    #=====================================================
    def instanciate (self,f):
    #=====================================================
        line  = f.readline()
        texts = line.split()
        self.index = int(texts[0])
        offs = 1
        n_load_records = 0
        if "L_REC" in line:
            n_load_records = int(texts[len(texts)-1])
        n_xyz = len(texts)-1-(n_load_records//max(1,n_load_records)) * 2

        for i in range (n_xyz):
            xi = float(texts [offs+i])
            self.xyz.append (xi)

        if len(self.xyz) < 3:
            self.xyz.append(0.0)

        if n_load_records > 0:
            line  = f.readline()
            texts = line.split()
            for text in texts:
                self.load_records.append (int(text))

    #=====================================================
    def get_xyz (self):
    # =====================================================
        return self.xyz

    #=====================================================
    def attach_element_to_adj_list (self,element):
    # =====================================================
        self.adj_elements.append (element)

    #=====================================================
    def get_adj_elements (self,ele_group_index=None,at_time = None):
    # =====================================================
        aux = []
        for element in self.adj_elements:
            take_it = True
            if ele_group_index != None:
                if element.get_group_index() != ele_group_index:
                    take_it = False
            if at_time != None:
                if not element.is_ON(at_time):
                    take_it = False
            if take_it:
                aux.append (element)
        return aux

    #=====================================================
    def clear_adj_elements (self):
    #=====================================================
        self.adj_elements = []

