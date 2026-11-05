"""Read and evaluate element existence functions.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Exf` model loads activation intervals from a ZSoil project and determines
whether an element is active at a requested physical time.
"""

class Exf ():
    def __init__(self, my_mesh):
        """Initialize the existence function object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        self.mesh_ref = my_mesh
        self.index = -1
        self.data = []

    def instanciate(self, f):
        """Load this object's data from the associated ZSoil project input.

        :param f: Open input stream positioned at the record to read.
        :type f: IO
        """
        line = f.readline()
        texts = line.split()
        self.index = int(texts[0])
        n = int(texts[1])

        for i in range(n):
            aux = []
            line = f.readline()
            texts = line.split()
            for k in range(len(texts)):
                aux.append(float(texts[k]))
            self.data.append(aux)

    def is_ON(self, time):
        """Return whether this object is active at the requested time.

        :param time: Physical time at which the value or state is requested.
        :type time: float
        :return: Whether this object is active at the requested time.
        :rtype: bool
        """
        for i in range(len(self.data)):
            data = self.data[i]
            t1 = data[0]
            t2 = data[1]
            if time < 1.0e-12:
                t1 = t1 - 2.0e-12
            if time > t1+1.0e-12 and time <= t2:
                return True
        return False
