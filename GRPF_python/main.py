from __future__ import print_function, division
from numpy import array, linspace, pi
import matlab.engine
import sublattice as sl
from sys import exit


if __name__ == "__main__":

	# background refractive index:
	nind = 1.5
	
	# Nanoparticle radius [nm]:
	rnp = 20

	# Lattice parameters (currently assumes 2 NP unit cell):
	# Note -- see Pocock et al. ACS Photonics 2018 for details.
	d    = 330.0
	beta = 1.0

	# Band structure parameters:
	kmin = 0.0
	kmax = pi/d
	knum = 15


	# Specify file names:
	template_file_name = "analysis_parameters_TEMPLATE.txt"
	pole_file_name     = "poles.out"
	root_file_name     = "roots.out"   
	 

	#########################	
	#	CODE:		#
	#########################

	
	# start Matlab engine:
	eng = matlab.engine.start_matlab()

	# Prepare root and pole output files:
	pole_file = open(pole_file_name, 'w')
	root_file = open(root_file_name, 'w')

	# Loop over in-plane Bloch vectors:
	k_vals = linspace(kmin, kmax, num=int(knum))
	for kp in k_vals:

		# prepare analysis parameter file:
		try:
			tmp_file = open(template_file_name, 'r')
			tmp_text = tmp_file.read()
			tmp_file.close()
		except IOError:
			print("\nProblem reading in analysis parameter template file.\n") 
			exit()
		#
		tmp_text_filled = tmp_text.format(kp, rnp, nind, d, beta*d/2)	
		tmp_filled_file = open("analysis_parameters.m", 'w')
		tmp_filled_file.write(tmp_text_filled)
		tmp_filled_file.close()		


		# Run Generalized Root and Pole Finder (GRPF):
		# Note: The function (named f) who's poles are desired is defined in pyfunc.py 
		print("Beginning GRPF ...")
		z_root, root_mult, z_pole, pole_mult = eng.GRPF_python(nargout=4)	
		print("GRPF finished.\n")


		# convert output from GRPF to Numpy arrays and convert all to shape [Nx1]:
		z_root    = array(z_root)
		root_mult = array(root_mult)
		z_pole    = array(z_pole)
		pole_mult = array(pole_mult)
		try:
			check = z_root.shape[0]
		except IndexError:
			z_root = array([z_root]).reshape((1,1)) 
		try:
			check = root_mult.shape[0]
		except IndexError:
			root_mult = array([root_mult]).reshape((1,1)) 
		try:
			check = z_pole.shape[0]
		except IndexError:
			z_pole = array([z_pole]).reshape((1,1)) 
		try:
			check = pole_mult.shape[0]
		except IndexError:
			pole_mult = array([pole_mult]).reshape((1,1)) 
		#print("Pole List:")
		#print(z_pole)
		#print("")
		#print(pole_mult)
		#print("")
		#print("Root List:")
		#print(z_root)
		#print("")
		#print(root_mult)
		#print("")
		#print(type(z_pole))
		#print(type(z_root))

		# log poles:
		try:
			for izpole, zpole in enumerate(z_pole[:,0]):
				pole_file.write("{0:f}\t{1:f}\t{2:f}\t{3:.2f}\n"
						.format(kp, zpole.real, zpole.imag, pole_mult[izpole,0]))
			#
		except IndexError:
			print("\n\t** No roots found in search domain**\n")


		# log roots:
		try:
			for iroot, zroot in enumerate(z_root[:,0]): 
				root_file.write("{0:f}\t{1:f}\t{2:f}\t{3:.2f}\n"
						.format(kp, zroot.real, zroot.imag, root_mult[iroot,0]))
			#
		except IndexError:
			print("\n\t** No poles in search domain**\n")
	#end loop over k_vals


	# close pole and root output files:
	pole_file.close()
	root_file.close() 

	
	# stop Matlab engine:
	eng.quit()
	print("Matlab engine stopped.")

	print("main.py finished.\n\n") 

#if









