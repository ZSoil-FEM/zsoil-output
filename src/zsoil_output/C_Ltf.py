import numpy as np

#=====================================================
class Ltf ():
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

        aux = []
        for i in range (n):
            line  = f.readline ()
            texts = line.split ()
            aux1 = []
            for k in range (len(texts)):
                aux1.append ( float (texts [k]))
            aux.append (aux1)
        self.data = aux


    #=====================================================
    def interpolate (self,time):
    #=====================================================
        # first make numpy arrays from lists
        size = len (self.data)
        if size == 0:
            return 0.0
        if time <= self.data [0][0]:
            return self.data [0][1]
        if time >= self.data [-1][0]:
            return self.data [-1][1]

        ltf_ti  = np.zeros (size)
        ltf_vi  = np.zeros(size)
        for i in range (size):
            ltf_ti [i] = self.data [i][0]
            ltf_vi [i] = self.data [i][1]
        return np.interp (time,ltf_ti,ltf_vi)


    #=====================================================
    def get_times_for_value (self, value):
    #=====================================================
        # first make numpy arrays from lists
        size = len (self.data)
        if size == 0:
            return []
        out = []
        for i in range (size-1):
            t1 = self.data [i][0]
            t2 = self.data [i+1][0]
            v1 = self.data [i][1]
            v2 = self.data [i+1][1]
            if t2-t1 > 1.0e-14:
                if value >= v1 and value < v2:
                    t = t1 + (value-v1)/(v2-v1) * (t2-t1)
                    out.append (t)
            else:
                if abs(value-v1) < 1.0e-14:
                    out.append (t1)
        #last value exception
        v2 = self.data [-1][1]
        if abs(value-v2) < 1.0e-14:
            t2 = self.data[-1][0]
            if not t2 in out:
                out.append (t2)
        return out