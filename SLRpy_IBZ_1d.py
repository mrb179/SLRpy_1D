from __future__ import division, print_function
from numpy import pi, sqrt, sin, cos, array, dot, linspace, matmul, cross, exp, conj, zeros, identity, loadtxt, arcsin, arccos, eye, trapz
from scipy import interpolate
from numpy.linalg import inv, solve
from pyfunc import *
from sublattice import *
import pylab as plt
import time
from sys import exit, argv

'''
CREATED: MRB Spring 2020

NOTES:
	- This was adapted from SLRpy_intBZ.py on Xenon (For Jaydeep Paper)
'''

if __name__ == "__main__":

	# Hard-coded input params #

	run_type  = 'site_NN' 			# BZint and site_NN are the only options tested
	NN_order  = 1  
	wmin      = 1240/2.6 			# minimum wavelength 
	wmax      = 1240/3.8 			# maximum wavelength 
	wnum      = 151 					# number of wavelengths 
	nind      = 1.0 				# background refractive index 
	d         = 400 				# chain periodicity [nm]
	t         = 50 					# dipole displacement along x-direction
	r_dip     = [0,1,0] 			# unit vector specifying dipole orientation
	pol_state = None 		
	kx_min    = 0.0 			# minimum k-value 
	kx_max    = pi/d 			# maximum k-value 
	kx_num    = 51
	BZ_fac    = 2  				# multiplicative factor if BZ symmetries are invoked
	long_polarizability_path   = "polarizability_Ag_80nmD_n1p0.out"
	short_polarizability_path  = "polarizability_Ag_80nmD_n1p0.out"
	out_file_name = "slrpy_1d_ibz.out"


	'''
	# Use this block when reading input params from SLRpy input file #	
	
	if len(argv) < 2:
		print("\n\tMust specify an input file as command line argument.\n")
		exit()

	run_settings = parse_input(argv[1])

	if run_settings == -1:
		print("\nError encountered parsing input file: {0:s}".format(argv[1]))
		print("Aborting.\n")
		exit()
	#if

	run_type = run_settings["type"]
	wmin = run_settings["wave_min"]
	wmax = run_settings["wave_max"]
	wnum = run_settings["wave_num"]
	nind = run_settings["nind"]
	lattice_list = run_settings["lattice_list"]
	source_type = run_settings["source_type"]
	pol_state   = run_settings["pol_state"]
	k_start = run_settings["k_start"]
	k_stop  = run_settings["k_stop"]
	k_num   = run_settings["k_num"]

	out_file_name_base = argv[1].split(".in")[0]
	out_file_name = out_file_name_base + ".out"
	'''


	#############
	#	CODE 	#
	#############

	# print summary of sub-lattice settings:
	#for sl in lattice_list:
	#	sl.info()


	# create lattice list (really just to hold polarizability data):
	lattice_list = [sublattice(None, None, [0,0,0], long_polarizability_path, short_polarizability_path)]


	plot_pol  = False
	plot_LS   = False
	plot_ext  = True
	save_data = False


	'''
	Note:

		- For simplicity, theta angle inputs are positive, but in calculations below, 
			theta -> -1*theta to ensure that as theta increases, so too does +k_{//}.
		- pol_state must be equal to either 'TE' or 'TM'
		- currently assumes refractive index = 1.0

		CHECK:
		- what happens if polaizability files have different wavelength spacing?

	'''

	###################################################
	###################################################

	# deal with sub-lattice polarizabilities:
	for sl in lattice_list:

		#################
		# 	POL_LONG	#
		#################
		'''
		polarizability along the long axis
		'''

		# path to file with polarizability data:
		polar_file = sl.pol_long_file

		# read-in polarizability data:
		try:
			pol_data = loadtxt(polar_file)
			# Make sure it is sorted by wavelength:
			pol_data = pol_data[pol_data[:,0].argsort()]
		except IOError:
			print("\nCould not open polarizability file\n")
			exit()
		#try

		# ensure wavelength bounds are contained within domain of polarizabiliy data:
		# generate wavelength array:
		wmin = float(wmin)
		wmax = float(wmax)
		wnum = int(wnum)
		wavemin = max(pol_data[0,0], wmin)
		wavemax = min(pol_data[-1,0],wmax)
		waves = linspace(wmin, wmax, wnum)

		print("")
		print("Min. wavelength [nm]: {0:f}".format(wmin))
		print("Max. wavelength [nm]: {0:f}".format(wmax))
		print("Num. wavelength: {0:d}".format(wnum))
		print("")

		# interpolate polarizability data:
		# create complex polarizability array:
		# 
		# apply scale factor to convert polarizability file data from [um^3] to [nm^3]:
		scale_pol = 10**9
		pol_arr = array([complex(pol_data[i,2]*scale_pol, pol_data[i,3]*scale_pol) for i in range(len(pol_data))])
		# invert:
		ipol_arr = 1./pol_arr
		ipol_data = zeros((len(ipol_arr),3))
		ipol_data[:,0] = pol_data[:,0]
		ipol_data[:,1] = ipol_arr.real
		ipol_data[:,2] = ipol_arr.imag
		#pol_real_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,1])
		#pol_imag_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,2])
		pol_real_func = interpolate.InterpolatedUnivariateSpline(ipol_data[:,0], ipol_data[:,1])
		pol_imag_func = interpolate.InterpolatedUnivariateSpline(ipol_data[:,0], ipol_data[:,2])
		pol_real = pol_real_func(waves)
		pol_imag = pol_imag_func(waves)
		inv_alpha_long = array([complex(pol_real[w], pol_imag[w]) for w in range(len(waves))])
		sl.ipol_long = inv_alpha_long


		#################
		# 	POL_SHORT	#
		#################
		'''
		polarizability along the short axis
		'''

		# path to file with polarizability data:
		polar_file = sl.pol_short_file

		# read-in polarizability data:
		try:
			pol_data = loadtxt(polar_file)
			# Make sure it is sorted by wavelength:
			pol_data = pol_data[pol_data[:,0].argsort()]
		except IOError:
			print("\nCould not open polarizability file\n")
			exit()
		#try

		# ensure wavelength bounds are contained within domain of polarizabiliy data:
		# generate wavelength array:
		wmin = float(wmin)
		wmax = float(wmax)
		wnum = int(wnum)
		wavemin = max(pol_data[0,0], wmin)
		wavemax = min(pol_data[-1,0],wmax)
		waves = linspace(wmin, wmax, wnum)

		print("")
		print("Min. wavelength [nm]: {0:f}".format(wmin))
		print("Max. wavelength [nm]: {0:f}".format(wmax))
		print("Num. wavelength: {0:d}".format(wnum))
		print("")

		# interpolate polarizability data:
		# create complex polarizability array:
		# 
		# apply scale factor to convert polarizability file data from [um^3] to [nm^3]:
		scale_pol = 10**9
		pol_arr = array([complex(pol_data[i,2]*scale_pol, pol_data[i,3]*scale_pol) for i in range(len(pol_data))])
		# invert:
		ipol_arr = 1./pol_arr
		ipol_data = zeros((len(ipol_arr),3))
		ipol_data[:,0] = pol_data[:,0]
		ipol_data[:,1] = ipol_arr.real
		ipol_data[:,2] = ipol_arr.imag
		#pol_real_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,1])
		#pol_imag_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,2])
		pol_real_func = interpolate.InterpolatedUnivariateSpline(ipol_data[:,0], ipol_data[:,1])
		pol_imag_func = interpolate.InterpolatedUnivariateSpline(ipol_data[:,0], ipol_data[:,2])
		pol_real = pol_real_func(waves)
		pol_imag = pol_imag_func(waves)
		inv_alpha_short = array([complex(pol_real[w], pol_imag[w]) for w in range(len(waves))])
		sl.ipol_short = inv_alpha_short
	#for


	# ready output files:
	out_file = open(out_file_name, 'w')

	if run_type == "dispersion":
		print("\n\t** Begining Electric Field Calculation **")
		# Construct k_parr:
		kx = linspace(k_start[0], k_stop[0], k_num) 
		ky = linspace(k_start[1], k_stop[1], k_num) 

		# loop over k-parallel:
		print("Percent k-values completed:")
		for i in range(k_num):
			print("\t{0:.2f}%\n".format(float(i/k_num)*100))
			# loop over wavelengths:
			for iw, w in enumerate(waves):

				kparr = array([kx[i], ky[i], 0.0])

				# determine incident angles (needs to be in deg. now since this is what build_GDyad expects):
				kparr_mag = sqrt(kx[i]**2 + ky[i]**2)
				k         = 2*pi*nind/w 
				theta     = (180./pi)*arcsin(kparr_mag/k) 
				if theta < 0.001:
					print("\n**** Note: theta < 0.001 deg. Defining phi = 0. ****\n")
					phi = 0.
				else:
					phi   = (180./pi)*arccos(kx[i]/kparr_mag)
				#if

				if kparr_mag > k:
					print("WARNING: kx = {0:f}, ky = {1:f} outside of light cone!".format(kx[i], ky[i]))

				# ENABLE following line to print theta, phi, kx, ky (DEBUGGING)
				#print("{0:f}\t{1:f} || {2:f}\t{3:f}\n".format(theta, phi, kx[i], ky[i]))
				

				# Build G-tensor:
				num_sub_latt = len(lattice_list)
				E_vec = zeros(3*num_sub_latt, dtype=complex)
				GD = zeros((3*num_sub_latt, 3*num_sub_latt), dtype=complex)
				I     = eye(3)
		
				for m in range(len(lattice_list)):
					for n in range(len(lattice_list)):
						if m == n :
							# build 3x3 diagonal block:
							sl = lattice_list[m]
							inv_alpha = array([sl.ipol_long[iw], sl.ipol_long[iw], sl.ipol_short[iw]])*I
							A = inv_alpha - build_GDyad_final(w, nind, [0,0,0], theta, phi, sl.a1, sl.a2, 5, include_origin=False)
							# fill suitable block of GD:
							GD[3*m,3*m]     = A[0,0]
							GD[3*m,3*m+1]   = A[0,1]
							GD[3*m,3*m+2]   = A[0,2]
							GD[3*m+1,3*m]   = A[1,0]
							GD[3*m+1,3*m+1] = A[1,1]
							GD[3*m+1,3*m+2] = A[1,2]
							GD[3*m+2,3*m]   = A[2,0]
							GD[3*m+2,3*m+1] = A[2,1]
							GD[3*m+2,3*m+2] = A[2,2]

							# calc E field incident on m_th sub-lattice:
							E_sl = E_field(sl.offset, kparr, theta, phi, pol_state)
							E_vec[3*m:3*m+3] = E_sl
							sl.E0 = E_sl

						else:
							# build 3x3 off-diagonal block:
							sl_m = lattice_list[m]	# row
							sl_n = lattice_list[n]	# column
							r_eval = sl_m.offset - sl_n.offset
							A = -1*build_GDyad_final(w, nind, r_eval, theta, phi, sl_n.a1, sl_n.a2, 5, include_origin=True)
							# fill suitable block of GD:
							GD[3*m,3*n]     = A[0,0]
							GD[3*m,3*n+1]   = A[0,1]
							GD[3*m,3*n+2]   = A[0,2]
							GD[3*m+1,3*n]   = A[1,0]
							GD[3*m+1,3*n+1] = A[1,1]
							GD[3*m+1,3*n+2] = A[1,2]
							GD[3*m+2,3*n]   = A[2,0]
							GD[3*m+2,3*n+1] = A[2,1]
							GD[3*m+2,3*n+2] = A[2,2]
						#
					#
				#
				
				P = 4*pi*solve(GD, E_vec)

				for j in range(num_sub_latt):
					lattice_list[j].pvec = P[j*3:3*j+3]
				#for

				# calculate extinction:
				ext = k*dot(conj(E_vec), P).imag

				# write output to file:
				# FORMAT:   0        1        2          3           4           5
				#			kx		ky 		theta 		phi 		wave 		ext
				out_line = "{0:f}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\n".format(kx[i], ky[i], theta, phi, w, ext)
				out_file.write(out_line)
			#for
		#for
	#if

	
	if run_type == "field":
		print("\n\t** Begining Electric Field Calculation **")
		if run_settings["field_calc_settings"] == None:
			print("Settings for electric field calculation not properly set. Aborting.")
		#
		fcs = run_settings["field_calc_settings"]
		field_file_name = out_file_name_base + "_Efield_{0:d}nm.out".format(int(fcs["wavelength"]))
		Efield_cross_sec(fcs, run_settings, lattice_list, save_name=field_file_name)
		print("Field calculation finished")
		print("")
	#if




	if run_type == "BZint":

		Ek_file  = open("Ek.out", 'w')
		alpha_k_file = open("alpha_k.out", 'w')

		# Note: 1d equivalents of these next 2 functions (build_T_block, and build_S_block) are defined in pyfunc.py
		'''
		def Tk(w_, sl_, re_, kx_, ky_, nind_ ):
			kparr_mag = sqrt(kx_**2 + ky_**2)
			k         = 2*pi*nind_/w_ 
			theta     = (180./pi)*arcsin(kparr_mag/k) 
			if theta < 0.001:
				print("\n**** Note: theta < 0.001 deg. Defining phi = 0. ****\n")
				phi = 0.
			else:
				phi   = (180./pi)*arccos(kx_/kparr_mag)
			#if

			T = build_GDyad_final(w, nind_, re_, theta, phi, sl_.a1, sl_.a2, 10, include_origin=True) 

			return T

		#def


		def Sk(w_, sl_, kx_, ky_, nind_ ):
			kparr_mag = sqrt(kx_**2 + ky_**2)
			k         = 2*pi*nind_/w_ 
			theta     = (180./pi)*arcsin(kparr_mag/k) 
			if theta < 0.001:
				print("\n**** Note: theta < 0.001 deg. Defining phi = 0. ****\n")
				phi = 0.
			else:
				phi   = (180./pi)*arccos(kx_/kparr_mag)
			#if

			S = build_GDyad_final(w, nind_, [0,0,0], theta, phi, sl_.a1, sl_.a2, 10, include_origin=False) 

			return S

		#def

		def make3x3(w_, iw_, sl_, re_, kx_, ky_, nind_):
			inv_alpha = array([sl.ipol_long[iw], sl.ipol_long[iw], sl.ipol_short[iw]])*I

			tk = Tk(w_, sl_, re_, kx_, ky_, nind_ )
			sk = Sk(w_, sl_, kx_, ky_, nind_ )
			B  = inv_alpha - sk 
			Binv = inv(B)
			A  = dot(conj(tk), Binv)
			AA = dot(A, tk)

			return AA 

		#def


		def generate_integrand(w_, iw_, sl_, re_, kx_, ky_, nind_):
			ikx = kx_.shape[0]
			iky = ky_.shape[0]
			B = zeros((ikx, iky, 3, 3), dtype=complex)
			Note:
				Because I'm only populating 1/8 of the BZ with non-zero
				values, the integrand will need to be multiplied by 8 at 
				the end! This also only works for the case of the cubc
				lattice and will need to be generalized later.
			for i in range(len(kx_)):
				for j in range(len(ky_)):
					if kx_[i] >= ky_[j]:
			 			B[i,j,0:3,0:3] = make3x3(w_, iw_, sl_, re_, kx_[i], ky_[j], nind_)
			return B

		#def
		'''

		# construct integrand at a given point in reciprocal space:
		def make3x3(w_, iw_, sl_, d_, t_, kx_, nind_):
			'''
			inputs:
				- w_  	: vacuum wavelength [nm]
				- iw_ 	: index of vacuum wavelength in wavelength list
				- sl_ 	: sub-lattice object (holds polarizability info.) 
				- t_  	: displacement of dipole along x-direc. relative to origin
				- kx_ 	: in-plane wavector [nm^-1]
				- nind_ : background refractive index
			returns:
				- AA    : 3x3 array to be integrated over 
			'''
			# lattice/interaction sum functions take frequencies (eV) as input:
			wc = 2*pi*hbar*c/w_

			inv_alpha = array([sl.ipol_long[iw_], sl.ipol_long[iw_], sl.ipol_short[iw_]])*I
			sk = build_S_block(wc, d_, kx_, nind_)
			tk_01 = build_T_block(wc, d_, t_, kx_, nind_)
			tk_10 = build_T_block(wc, d_, -1*t_, kx_, nind_)
			B  = inv_alpha - sk 
			Binv = inv(B)
			A  = dot(tk_10, Binv)
			AA = dot(A, tk_01)

			# write data to respective files:
			abz = (2*pi/d)
			Ek_file.write("{0:.4f}\t{1:.6e}\t{2:.6e}\t{3:.6e}\n".format(wc, kx_, abz*real(tk_01[1,1]), abz*imag(tk_01[1,1]) ))
			alpha_k_file.write("{0:.4f}\t{1:.6e}\t{2:.6e}\t{3:.6e}\n".format(wc, kx_, real(Binv[1,1]), imag(Binv[1,1]) ))

			return AA 

		#def


		def generate_integrand(w_, iw_, sl_, d_, t_, kx_, nind_):
			'''
			inputs:
				- w_  	: vacuum wavelength [nm]
				- iw_ 	: index of vacuum wavelength in wavelength list
				- sl_ 	: sub-lattice object (holds polarizability info.) 
				- t_  	: displacement of dipole along x-direc. relative to origin
				- kx_ 	: in-plane wavector [nm^-1]
				- nind_ : background refractive index
			returns:
				- B    : [num_kx, 3, 3] array of integrand (3x3) evaluations 
			'''
			ikx = kx_.shape[0]
			B = zeros((ikx, 3, 3), dtype=complex)
			'''
			Note:
				Because I'm only populating 1/2 of the BZ with non-zero
				values, the integrand will need to be multiplied by 2 at 
				the end!
			'''
			for i in range(len(kx_)):
			 	B[i,0:3,0:3] = make3x3(w_, iw_, sl_, d_, t_, kx_[i], nind_)

			return B

		#def


		print("### BZ integrate ###")

		# define identity matrix:
		I = eye(3,3)

		# area (length) of Brillouin Zone:
		ABZ = 2*pi/d

		# create normalized unit vector specifying dipole orientation:
		r_dip = array(r_dip)/sqrt(dot(r_dip,r_dip))

		# initialize array of kx values:
		kx = linspace(kx_min, kx_max, num=kx_num) 
		dkx = (kx_max - kx_min)/(kx_num-1)
		print("\t- k_x resolution: {0:f}".format(dkx))

		# temporary holding containers for 0,0 elements of eta:
		tmp_list_r = []
		tmp_list_i = []

		# begin timer:
		t0 = time.time()

		# loop over wavelengths:
		for iw, w in enumerate(waves):

			t_start = time.time()

			print("Wavelength Number {0:d}: {1:.1f}% complete".format(iw, iw/len(waves)*100))

			eta = zeros((3,3), dtype=complex)

			imat = generate_integrand(w, iw, sl, d, t, kx, nind)

			# loop over 3x3 Cartesian components and integrate: 
			print("**WARNING**")
			print("Currently ONLY looking at i=0, j=0 element of eta! See line ~ 440")
			for i in range(1):
				for j in range(1):
					# this is just a trick to get yy component. Comment out for xx component:
					i = 1
					j = 1

					# integrate over first axis:
					#int_over_x = trapz(imat[:,:,i,j], x=kx, dx=dkx, axis=0)
					int_val = trapz(imat[:,i,j], x=kx, dx=dkx, axis=0)

					# integrate over second axis:
					#int_val = trapz(int_over_x, x=ky, dx=dky, axis=-1)

					eta[i,j] = BZ_fac * int_val/ABZ  # note: BZ_fac can be used to account for BZ symmetry
				#
			#

			# This block for xx component:
			tmp_list_r.append(eta[1,1].real)
			tmp_list_i.append(eta[1,1].imag)
                        
			# This block for xx component:
			#tmp_list_r.append(eta[0,0].real)
			#tmp_list_i.append(eta[0,0].imag)

			# calculate Purcell Factor:
			k = 2*pi*nind/w
			prefactor = 6*pi/(nind*k**3)
			PF = 1 + prefactor*tmp_list_i[-1]




			out_file.write("{0:.2f}\t{1:.6f}\t{2:.6f}\t{3:.3f}\n".format(w, eta[0,0].real, eta[0,0].imag, PF))

			t_end = time.time()
			print("Time elapsed for wavelength {0:d}: {1:.1f}".format(iw, t_end - t_start))

		#for_waves

		Ek_file.close()
		alpha_k_file.close()
	
	#if




	if run_type == "site_NN":

		NN_order = int(NN_order)
		print("### Site Space -- Nearest Neighbor ###")
		print("\tNN Order = {0:d}\n".format(NN_order))
		out_file.write("# Site-Space NN_order = {0:d}".format(NN_order))

		# Determine total number of dipoles and degrees of freedom:
		num_dof = 6*(NN_order + 1)
		num_p   = 2*(NN_order + 1)

		# define identity matrix:
		I = eye(3,3)

		# create normalized unit vector specifying dipole orientation: 
		#r_dip = array(r_dip)/sqrt(dot(r_dip,r_dip)) 	# orientation
		#r_d   = array([t, 0, 0]) 						# position
		r_dip = array([0,0,1]) 	# orientation
		r_d   = array([0, 0, t]) 						# position
		print("Dipole Orientation:")
		print(r_dip)
		print("Dipole Position:")
		print(r_d)

		# create list of NP positions (length = 2*NN_order + 1)
		r_NP_list = [ array([n*d, 0, 0]) for n in range(-NN_order,NN_order+1) ]
		print("\nNP positions:")
		print(r_NP_list)

		# create sublattice object to hold emitter dipole 
		# For now using same polarizability as NP (physically, answer should not depend on this ... (should test))
		#beta_obj = sublattice(None, None, [0,0,0], long_polarizability_path, short_polarizability_path)

		# construct E_vec to only drive emitter dipole:
		E_vec = zeros((num_dof,1), dtype=complex)
		E_vec[0:3,0] = r_dip
		print("\nDriving Field")
		print(E_vec)

		# container to hold Purcell factor data:
		PF_NN   = []
		Es_arr  = []

		# loop over wavelengths:
		for iw, w in enumerate(waves):

			# construct M tensor:
			M = zeros((num_dof, num_dof), dtype=complex)
			for i in range(num_p):
				for j in range(num_p):
					if i == j :
						# diagonal (polarizability) block:
						if i == 0:
							# emitter polarizability:
							sl = lattice_list[0]
							A = array([sl.ipol_long[iw], sl.ipol_long[iw], sl.ipol_short[iw]])*I

							# fill suitable block of GD:
							M[3*i,3*i]     = A[0,0]
							M[3*i,3*i+1]   = A[0,1]
							M[3*i,3*i+2]   = A[0,2]
							M[3*i+1,3*i]   = A[1,0]
							M[3*i+1,3*i+1] = A[1,1]
							M[3*i+1,3*i+2] = A[1,2]
							M[3*i+2,3*i]   = A[2,0]
							M[3*i+2,3*i+1] = A[2,1]
							M[3*i+2,3*i+2] = A[2,2]
						else:
							# NP polarizability:
							sl = lattice_list[0]
							A = array([sl.ipol_long[iw], sl.ipol_long[iw], sl.ipol_short[iw]])*I
	
							# fill suitable block of GD:
							M[3*i,3*i]     = A[0,0]
							M[3*i,3*i+1]   = A[0,1]
							M[3*i,3*i+2]   = A[0,2]
							M[3*i+1,3*i]   = A[1,0]
							M[3*i+1,3*i+1] = A[1,1]
							M[3*i+1,3*i+2] = A[1,2]
							M[3*i+2,3*i]   = A[2,0]
							M[3*i+2,3*i+1] = A[2,1]
							M[3*i+2,3*i+2] = A[2,2]
							
					else:
						# off-diagonal (coupling) block:
						if i == 0:
							reval   = r_d
						else:
							reval   = r_NP_list[i-1]
						if j == 0:
							rsource = r_d
						else:
							print(j)
							rsource = r_NP_list[j-1] 

						A  = -1*GD_single(reval, rsource, w, nind)
	
						M[3*i,3*j]     = A[0,0]
						M[3*i,3*j+1]   = A[0,1]
						M[3*i,3*j+2]   = A[0,2]
						M[3*i+1,3*j]   = A[1,0]
						M[3*i+1,3*j+1] = A[1,1]
						M[3*i+1,3*j+2] = A[1,2]
						M[3*i+2,3*j]   = A[2,0]
						M[3*i+2,3*j+1] = A[2,1]
						M[3*i+2,3*j+2] = A[2,2]
				#for
			# end: M_tensor construction loop

			# solve matrix equation for dipole moments:
			P = solve(M, E_vec)

			# pick out the moment of the driven emitter (should have shape (3,1)):
			d_emitter = P[0:3,0]

			# compute scattered field at the location of the emitter:
			k = 2*pi*nind/w
			Es = complex(0,0)
			for n in range(num_p - 1):
				Es += dot(GD_single(r_d, r_NP_list[n], w, nind), P[3*(n+1):3*(n+2),0])
			#for 
			Es_arr.append(Es)

			# compute Purcell factor:
			prefactor = 6*pi/(nind*k**3)
			#PF = 1 + prefactor*dot( conj(d_emitter.T), Es ).imag[0]/(sqrt(dot(conj(d_emitter.T), d_emitter))).real
			PF = 1 + prefactor*dot( conj(d_emitter.T), Es ).imag /((dot(conj(d_emitter.T), d_emitter).real))
			PF_NN.append(PF)

			# write to output file:
			out_file.write("{0:.2f}\t{1:.4f}\n".format(w, PF))

		# end: wavelength loop

		# Optional Plotting:
		if False:
			print(PF_NN)
			plt.figure(1)
			plt.plot(1240/waves, PF_NN, color="red", label="NN = {0:d}".format(NN_order))
			plt.xlabel("Energy (eV)")
			plt.ylabel("Purcell Factor")
			#plt.ylim((0,40))
			plt.show()
		#if

	#if




	# house-keeping:
	out_file.close()

	print("Finished Successfully")
	print("Data written to: {0:s}".format(out_file_name))

#if
