"""Represent material definitions read from a ZSoil project.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Material` model stores material identifiers, labels, property blocks, and lookup
information used by mesh and result-processing operations.
"""

class Material ():
    def __init__(self, my_mesh):
        """Initialize the material object for the supplied project data.

        :param my_mesh: Mesh that owns this object.
        :type my_mesh: Mesh
        """
        self.mesh_ref = my_mesh

        self.index_in_dat = 0
        self.id = ""
        self.label = ""
        self.exf_index = 0
        self.unlf_index = 0
        self.index_in_inp = 0
        self.phase = 0
        self.layer = 0

        self.dict = {'MAIN': 0, 'DENS': 1, 'ELAS': 2, 'CREEP': 3, 'NONL': 4, 'HEAT': 5, 'HUMID': 6,
                     'INIS': 7, 'STAB': 8, 'DAMP': 9, 'DISC': 10, 'XXXX': 11, 'GEOM': 12, 'FLOW': 13}
        self.dict_inv = {v: k for k, v in list(self.dict.items())}

        self.data = []
        self.data_ltf = []
        self.data_spe = []
        self.data_evf = []

        for i in range(len(self.dict)):
            self.data.append([])
            self.data_ltf.append([])
            self.data_spe.append([])
            self.data_evf.append([])

    def get_copy(self):
        """Return an independent copy of this material definition.

        :return: An independent copy of this material definition.
        :rtype: Material
        """
        mat = Material(self.mesh_ref)

        mat.index_in_dat = self.index_in_dat
        mat.id = self.id
        mat.label = self.label
        mat.exf_index = self.exf_index
        mat.unlf_index = self.unlf_index
        mat.index_in_inp = self.index_in_inp
        mat.phase = self.phase
        mat.layer = self.layer
        mat.dict = {'MAIN': 0, 'DENS': 1, 'ELAS': 2, 'CREEP': 3, 'NONL': 4, 'HEAT': 5, 'HUMID': 6,
                    'INIS': 7, 'STAB': 8, 'DAMP': 9, 'DISC': 10, 'XXXX': 11, 'GEOM': 12, 'FLOW': 13}
        mat.dict_inv = {v: k for k, v in list(self.dict.items())}

        mat.data = []
        for item in self.data:
            mat.data.append(item[:])

        mat.data_ltf = []
        for item in self.data_ltf:
            mat.data_ltf.append(item[:])

        mat.data_spe = []
        for item in self.data_spe:
            mat.data_spe.append(item[:])

        mat.data_evf = []
        for item in self.data_evf:
            mat.data_evf.append(item[:])

        return mat

    def instanciate(self, f):
        """Load this object's data from the associated ZSoil project input.

        :param f: Open input stream positioned at the record to read.
        :type f: IO
        """
        line = f.readline()
        texts = line.split()

        self.label = line[90:]
        self.label = self.label.replace("\n", "")  # remove end of line if needed
        self.label = self.label.lstrip()  # cancel all spaces at the beginning
        self.index_in_dat = int(texts[0])
        #print "reading material : ",self.index_in_dat,"--->",self.label
        #print len(texts)
        #print texts
        self.id = texts[1]
        self.exf_index = int(texts[2])
        self.unlf_index = int(texts[3])
        self.index_in_inp = 0
        if len(texts) > 4:
            if len(texts[4]) > 0:
                self.index_in_inp = int(texts[4])

        if len(texts) > 5:
            if len(texts[5]) > 0:
                self.phase = int(texts[5])
        if len(texts) > 6:
            if len(texts[6]) > 0:
                self.layer = int(texts[6])

        while True:  # little dangerous loop but for correct data set is ok
            last_file_pos = f.tell()
            line = f.readline()
            #print "read line :",line
            # first check for the group label
            group_id = self.group_label_is_in_string(line)
            #next_line = ""
            #print "group id--->",group_id,"----"

            if group_id != "":
                index = self.dict[group_id]
                texts = line.split()
                bit_code = 0
                if len(texts) > 1:
                    bit_code = int(texts[1])
                #set_of_lines = []
                done = False

                next_group = False
                next_model = False
                end_materials = False

                # read group data
                tmp = []
                while not next_group and not next_model and not end_materials:
                    last_file_pos = f.tell()
                    aux_string = f.readline()
                    #print "read line xxxxxxx:",aux_string
                    group_id = self.group_label_is_in_string(aux_string)
                    if group_id != "":
                        next_group = True
                        f.seek(last_file_pos)
                    elif "LINK" in aux_string or "NODE" in aux_string or "PILES" in aux_string or \
                         "NAILS" in aux_string or "ANCHOR_HEADS" in aux_string or \
                         "CABLES" in aux_string or "AUX_ELE" in aux_string:
                        end_materials = True
                        f.seek(last_file_pos)
                    else:
                        aaa = aux_string[:2]
                        #print "---",aux_string[0:1],"---",len(aux_string[0:1])
                        if aaa == "  ":
                            next_model = True
                            f.seek(last_file_pos)
                    if not next_model and not next_group and not end_materials:
                        tmp.append(aux_string)

                if len(tmp) > 0:
                    #data = self.data [index]
                    self.repack_to_data(bit_code, self.data[index], self.data_ltf[index],
                                        self.data_spe[index], self.data_evf[index], tmp)

                if next_model or end_materials:
                    return
            else:
                f.seek(last_file_pos)
                return

    def repack_to_data(self, bit_code, data, data_ltf, data_spe, data_evf, tmp):
        """Split raw material parameters into their standard and optional data arrays.

        :param bit_code: Bit code used by this operation.
        :type bit_code: int
        :param data: Data used by this operation.
        :type data: list
        :param data_ltf: Data load-time function used by this operation.
        :type data_ltf: list
        :param data_spe: Data spe used by this operation.
        :type data_spe: list
        :param data_evf: Data evf used by this operation.
        :type data_evf: list
        :param tmp: Tmp used by this operation.
        :type tmp: Sequence[str]
        """
        n = 1
        if bit_code & 1:
            n = n + 1
        if bit_code & 2:
            n = n + 1
        if bit_code & 4:
            n = n + 1

        nlines = len(tmp) // max(1, n)  # per group of parameters

        offs = 0
        # standard data
        for i in range(nlines):
            line = tmp[i+offs]
            texts = line.split()
            for k in range(len(texts)):
                a = float(texts[k])
                data.append(a)
        offs = offs + nlines

        if bit_code & 1:  # load time functions to parameters
            for i in range(nlines):
                line = tmp[i+offs]
                texts = line.split()
                for k in range(len(texts)):
                    a = int(float(texts[k]))
                    data_ltf.append(a)
            offs = offs + nlines

        if bit_code & 2:  # superelements
            for i in range(nlines):
                line = tmp[i+offs]
                texts = line.split()
                for k in range(len(texts)):
                    a = int(float(texts[k]))
                    data_spe.append(a)
            offs = offs + nlines

        if bit_code & 4:  # evolution functions
            for i in range(nlines):
                line = tmp[i+offs]
                texts = line.split()
                for k in range(len(texts)):
                    a = int(float(texts[k]))
                    data_evf.append(a)
            offs = offs + nlines

    def group_label_is_in_string(self, string):
        """Return the recognized element-group label at the start of the supplied text.

        :param string: String used by this operation.
        :type string: str
        :return: The recognized element-group label at the start of the supplied text.
        :rtype: str
        """
        if (len(string) == 0):
            return ""

        texts = string.split()
        if (len(texts) == 0):
            return ""

        tmp = (texts[0]).strip()
        for i in range(len(self.dict)):
            group_label = self.dict_inv[i]
            if group_label == tmp:
                return group_label
        return ""

    def get_user_index(self):
        """Return the material index used in the input file.

        :return: The material index used in the input file.
        :rtype: int
        """
        return self.index_in_inp

    def get_dat_index(self):
        """Return the material index used in the DAT data block.

        :return: The material index used in the DAT data block.
        :rtype: int
        """
        return self.index_in_dat
