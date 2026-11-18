from .C_Material import *
from .C_Exf      import *
from .C_Element_T3_generic import *
from .C_Element_Q4_generic import *
from .C_Element_L2D_generic import*
# from numpy      import *
import numpy    as np
from math import *


#=====================================================
class Element ():
#=====================================================

	GROUP_CONTINUUM = 0 #implemented in SDK
	GROUP_SHELL     = 1 #implemented
	GROUP_TRUSS     = 2 #implemented
	GROUP_BEAM      = 3 #implemented
	GROUP_DSC       = 4 #not implemented
	GROUP_CONTACT   = 5 #implemented
	GROUP_SEEPAGE   = 6 #not implemented
	GROUP_CONVECTION= 7 #not implemented
	GROUP_FLVOL     = 8 #not implemented
	GROUP_PVOL      = 9 #not implemented
	GROUP_INFINITE  = 10 #not implemented
	GROUP_MEMBRANE  = 11 #implemented
	GROUP_RING      = 12 #not implemented
	GROUP_NSCNT     = 13 #not implemented
	GROUP_MASS      = 14 #not implemented
	GROUP_VISC_DAMP = 15 #not implemented
	GROUP_LIN_HINGE = 16 #not implemented
	GROUP_HEAT_EXCH = 17 #implemented
	GROUP_MAX       = 18

	group_dict_inv = {"BEAMS"    : GROUP_BEAM, \
			  	  "TRUSSES"  : GROUP_TRUSS, \
				  "VOLUMICS" : GROUP_CONTINUUM, \
				  "SHELLS"   : GROUP_SHELL, \
				  "CONTACT"  : GROUP_CONTACT, \
				  "MEMBRANE" : GROUP_MEMBRANE,\
                  "HEAT_EXCH": GROUP_HEAT_EXCH}


	group_dict = {v: k for k, v in list(group_dict_inv.items())}
	#now we can get key strings to element groups via dictionary by
	#calling    Element.group_dict [Element.GROUP_CONTINUUM] instead of "VOLUMICS"
	# calling   Element.group_dict [Element.GROUP_SHELL] instead of "SHELLS"
	# etc....

	#=====================================================
	def __init__ (self,my_mesh):
	#=====================================================
		self.mesh_ref = my_mesh
		self.index    = -1
		self.label    = ""
		self.group    = "UNKNOWN"
		self.rsl_ext  = ""
		self.rsl_lay_ext = None
		self.xsiGP    = []
		self.material_index = -1
		self.nodes    = []
		self.load_records = []
		self.load_faces   = []
		self.seek_in_rsl = 0
		self.seek_in_lay_rsl = 0
		self.nen      = 0
		self.ref_ele  = None


	#=====================================================
	def instanciate (self,f):
	#=====================================================
		line  = f.readline()
		texts = line.split()
		self.index = int(texts[0])
		self.label =     (texts[1]).strip()
		self.material_index = int (texts[2])
		offs = 3
		n_load_records = 0
		if "L_REC" in line:
			n_load_records = int(texts[4])
			offs = 5
		n_nodes = len(texts)-3-(n_load_records//max(1,n_load_records)) * 2
		for i in range (n_nodes):
			node = int(texts [offs+i])
			self.nodes.append (node)
		if n_load_records > 0:
			n_lines_with_load_records = n_load_records // 8 + 1
			if n_load_records % 8 ==0:
				n_lines_with_load_records = n_lines_with_load_records - 1
			#print "n_lines_with_load_records",n_lines_with_load_records,n_load_records
			for j in range (n_lines_with_load_records):
				line  = f.readline()
				texts = line.split()
				count = 0
				for text in texts:
					if count % 2 == 0:
						self.load_records.append (int(text))
					else:
						self.load_faces.append (int(text))
					count = count + 1


	#=====================================================
	def give_exf (self):
	#=====================================================
		mat = self.mesh_ref.materials [self.material_index-1]
		if mat.exf_index == 0:
			return None
		else:
			return self.mesh_ref.exfs [mat.exf_index-1]

	#=====================================================
	def is_ON (self,time):
	#=====================================================
		exf = self.give_exf()
		return exf.is_ON (time)

	#=====================================================
	def get_ngaus (self):
	#=====================================================
		return len(self.xsiGP)

	#=====================================================
	def get_ele_coord (self):
	#=====================================================
		ele_coord = np.zeros (self.nen * 3).reshape (self.nen,3)
		for i in range (self.nen):
			node = self.mesh_ref.nodes [self.nodes [i]-1]
			for j in range (3):
				ele_coord [i,j] = node.xyz [j]
		return ele_coord

	#=====================================================
	def get_ele_dim (self):
	#=====================================================
		return None


	#=====================================================
	def get_nlayers (self):
	#=====================================================
		return 0,0

	#=====================================================
	def get_group (self):
	#=====================================================
		return self.group

	#=====================================================
	def get_group_index (self):
	#=====================================================
		if self.group not in list(self.group_dict_inv.keys()):
			return -1
		else:
			return self.group_dict_inv [self.group]

	#=====================================================
	def get_material_index (self):
	#=====================================================
		return self.material_index


	#=====================================================
	def get_nodes (self):
	#=====================================================
		nodes = []
		for node_index in self.nodes:
			nodes.append (self.mesh_ref.nodes [node_index-1])
		return nodes


	#=====================================================
	def get_center (self):
	#=====================================================
		xyz = [0.0,0.0,0.0]
		for i in range (self.nen):
			node = self.mesh_ref.nodes [self.nodes [i]-1]
			for j in range (len(node.xyz)):
				xyz [j] = xyz [j] + node.xyz [j] / self.nen

		return xyz

	#=====================================================
	def get_size_along_dir (self,vec):
	#=====================================================
		versor = np.array (vec)
		length = sqrt(np.dot (versor,versor))
		versor = versor * (1.0/length)
		min_proj =  1.0e38
		max_proj = -1.0e38
		for i in range (self.nen):
			node = self.mesh_ref.nodes [self.nodes [i]-1]
			xyz  = np.array (node.xyz)
			aux  = np.dot (xyz,versor)
			min_proj = min ( min_proj,aux)
			max_proj = max ( max_proj,aux)

		return max_proj-min_proj

	#=====================================================
	def get_surface_along_dir (self,vec):
	#=====================================================
		#return area of face with normalclosest to the given vector
		versor = np.array (vec)
		length = sqrt(np.dot (versor,versor))
		versor = versor * (1.0 / length)

		n_faces = self.ref_ele.get_nr_of_faces ()

		if n_faces <= 0:
			return None

		ele_coord = self.get_ele_coord ()
		x = np.zeros(4 * 3).reshape(4,3)

		max_proj = -1.0e38

		face_best = None
		x_best    = np.zeros(4 * 3).reshape(4,3)

		for i in range (n_faces):
			face_type, node_indices = self.ref_ele.get_face_info (i+1)
			if face_type == 'T3':
				face_ele = Element_T3_generic ()
			elif face_type == 'Q4':
				face_ele = Element_Q4_generic()
			elif face_type == 'L2':
				face_ele = Element_L2D_generic()
			else:
				print('error in C_Element:get_surface_along_dir')

			for j, node_index_1 in enumerate(node_indices):
				x [j,0] = ele_coord [node_index_1-1,0]
				x [j,1] = ele_coord [node_index_1-1,1]
				x [j,2] = ele_coord [node_index_1-1,2]

			n = face_ele.get_normal (x)

			proj = np.dot (n,versor)
			if proj > max_proj:
				face_best = face_ele
				max_proj  = proj
				x_best [:,:] = x [:,:]
				if max_proj > 1.0-1.0e-6:
					break

		S = face_best.get_surface (x_best,Element_Surf_generic.QUADR_STD)
		return S


	#=====================================================
	def get_surface (self,quadr_enum):
	#=====================================================
		pass

	#=====================================================
	def get_volume (self,quadr_enum):
	#=====================================================
		pass

	#=====================================================
	def get_length (self,quadr_enum):
	#=====================================================
		pass

	# #=====================================================
	# def get_intersection_with_plane (self):
	# #=====================================================
	#     return None

	#=====================================================
	def get_ref_ele (self):
	#=====================================================
		return self.ref_ele
