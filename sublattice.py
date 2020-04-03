from __future__ import division, print_function
from numpy import array
from matplotlib import pyplot as plt

'''
This class will help streamline the multi-sub-lattice problem
'''

class sublattice(object):
	
	def __init__(self, a_1, a_2, offset, along, ashort):
		self.a1 = a_1
		self.a2 = a_2
		self.offset = offset
		self.pol_long_file  = along 
		self.pol_short_file = ashort
		self.ipol_long  = None 
		self.ipol_short = None 
		self.pvec = None
		self.E0   = None
	#def

	def info(self):
		print("")
		print("SUB-LATTICE:")
		print("\ta1 = < {0:.2f}, {1:.2f}, {2:.2f} >".format(self.a1[0], self.a1[1], self.a1[2]))
		print("\ta2 = < {0:.2f}, {1:.2f}, {2:.2f} >".format(self.a2[0], self.a2[1], self.a2[2]))
		print("\tr0 = < {0:.2f}, {1:.2f}, {2:.2f} >".format(self.offset[0], self.offset[1], self.offset[2]))
		print("\tpol file (long)  = {0:s}".format(self.pol_long_file))
		print("\tpol file (short) = {0:s}".format(self.pol_short_file))
		print("")
	#def

#class


class finite_sublattice(object):
	
	def __init__(self, desc, offset, along, ashort, *args):
		self.type = desc
		self.offset = offset
		self.pos = self.gen_lattice(*args)
		self.pol_long_file  = along 
		self.pol_short_file = ashort
		self.ipol_long  = None 
		self.ipol_short = None 
		self.pvec = None
		self.E0   = None
	#def


	def gen_lattice(self, *args):

		if self.type == "rect":
			'''
			In this case, *args should be: px, py, nx, ny
				- px: lattice constant along the x-direction
				- py: lattice constant along the y-direction
				- nx: number of sites along the x-direction
				- ny: number of sites along the y-direction
			'''

			px = args[0]
			py = args[1]
			nx = args[2]
			ny = args[3]
			plot = False


			pos = []
			CMx = (nx-1)*px/2
			CMy = (ny-1)*py/2
			for i in range(nx):
				for j in range(ny):
					pos.append([px*i - CMx, py*j - CMy, 0.0])
				#
			#
			pos = array(pos)

			if plot:
				x_vals = [pos[i][0] for i in range(pos.shape[0])]
				y_vals = [pos[i][1] for i in range(pos.shape[0])]
				plt.figure(1)
				plt.scatter(x_vals, y_vals)
				plt.xlabel("x (nm)")
				plt.ylabel("y nm")
				plt.show()
			#
			return pos

		else:
			print("\n\tCurrently only rectangular finite lattices supported.\n")
			print("Aborting.")
			exit()
	#def


	def info(self):
		print("")
		print("SUB-LATTICE (finite):")
		print("\tType = {0:s}".format(self.type))
		print("Number of sites: {0:d}".format(self.pos.shape[0]))
		print("\tr0 = < {0:.2f}, {1:.2f}, {2:.2f} >".format(self.offset[0], self.offset[1], self.offset[2]))
		print("\tpol file (long)  = {0:s}".format(self.pol_long_file))
		print("\tpol file (short) = {0:s}".format(self.pol_short_file))
		print("")
	#def

#class