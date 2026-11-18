#=====================================================
class Exf ():
#=====================================================

    #=====================================================
    def __init__ (self,my_mesh):
    #=====================================================
        self.mesh_ref = my_mesh
        self.index    = -1
        self.data     = []

    #=====================================================
    def instanciate (self,f):
    #=====================================================
        line  = f.readline ()
        texts = line.split ()
        self.index = int(texts[0])
        n  = int(texts[1])

        for i in range (n):
            aux = []
            line  = f.readline ()
            texts = line.split ()
            for k in range (len(texts)):
                aux.append ( float (texts [k]))
            self.data.append (aux)

    #=====================================================
    def is_ON (self,time):
    #=====================================================

        for i in range (len(self.data)):
            data = self.data [i]
            t1   = data [0]
            t2   = data [1]
            if time < 1.0e-12:
                t1 = t1 -2.0e-12
            if time > t1+1.0e-12 and time <= t2:
                return True
        return False

