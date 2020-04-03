from __future__ import division, print_function
from numpy import pi, sqrt, sin, cos, array, dot, linspace, matmul, cross, exp, conj, zeros, identity, loadtxt, arcsin, arccos, eye, trapz
from scipy import interpolate
from numpy.linalg import inv, solve
#from sum_dyadic_test import *
from pyfunc import *
#from SLRpy_input import *
from sublattice import *
import pylab as plt
import time
from sys import exit, argv

'''
CREATED: MRB Spring 2020

NOTES:

'''

if __name__ == "__main__":

	run_type = disp
	wmin = 
	wmax = 
	wnum = 
	nind = 
	pol_state   = run_settings["pol_state"]
	k_start = run_settings["k_start"]
	k_stop  = run_settings["k_stop"]
	k_num   = run_settings["k_num"]
	'''

	out_file_name_base = argv[1].split(".in")[0]
	out_file_name = out_file_name_base + ".out"
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
	'''

	out_file_name_base = argv[1].split(".in")[0]
	out_file_name = out_file_name_base + ".out"

	# print summary of sub-lattice settings:
	for sl in lattice_list:
		sl.info()


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

		# Testing:
		if run_type == "field":
			waves = [run_settings["field_calc_settings"]["wavelength"]]
		# End Testing

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
		pol_real = pol_real_func(waves)*4*pi
		pol_imag = pol_imag_func(waves)*4*pi
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
		pol_real = pol_real_func(waves)*4*pi
		pol_imag = pol_imag_func(waves)*4*pi
		inv_alpha_short = array([complex(pol_real[w], pol_imag[w]) for w in range(len(waves))])
		sl.ipol_short = inv_alpha_short
	#for


	# ready output file:
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

		I     = eye(3)

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
			#print(kx_.shape)
			#print(ky_.shape)
			ikx = kx_.shape[0]
			iky = ky_.shape[0]
			B = zeros((ikx, iky, 3, 3), dtype=complex)
			#print("")
			#print("shape of B:")
			#print(B.shape)
			#print("")
			'''
			Note:
				Because I'm only populating 1/8 of the BZ with non-zero
				values, the integrand will need to be multiplied by 8 at 
				the end! This also only works for the case of the cubc
				lattice and will need to be generalized later.
			'''
			for i in range(len(kx_)):
				for j in range(len(ky_)):
					if kx_[i] >= ky_[j]:
			 			B[i,j,0:3,0:3] = make3x3(w_, iw_, sl_, re_, kx_[i], ky_[j], nind_)
			return B

		#def


		print("### BZ integrate ###")

		# emitter position [nm] w.r.t. origin at center of a NP: 
		rd = array([50., 0., 0.])

		# Reciprocal Space:
		b1, b2 = recipLatt(sl.a1, sl.a2)
		ABZ = norm(cross(b1, b2))

		'''
		## Incident Field Polarization ##
		'''
		E0 = array([1,0,0])


		# set integration boundaries in k-space:
		kx_min = 0.0
		kx_max = pi/450.
		kx_num = 75
		ky_min = 0.0
		ky_max = pi/450.
		ky_num = 75

		kx = linspace(kx_min, kx_max, num=kx_num) 
		ky = linspace(ky_min, ky_max, num=ky_num)
		dkx = (kx_max - kx_min)/(kx_num-1)
		print("\t- k_x resolution: {0:f}".format(dkx))
		dky = (ky_max - ky_min)/(ky_num-1)
		print("\t- k_y resolution: {0:f}".format(dky))

		# temporary holding containers for 0,0 elements of eta:
		tmp_list_r = []
		tmp_list_i = []

		# begin timer:
		t0 = time.time()

		# loop over wavelengths:
		for iw, w in enumerate(waves):

			t_start = time.time()

			print("Wavelength Number {0:d}: {1:.1f}% complete".format(iw, iw/len(waves)))

			eta = zeros((3,3), dtype=complex)

			imat = generate_integrand(w, iw, sl, rd, kx, ky, nind)

			# loop over 3x3 Cartesian components and integrate: 
			print("**WARNING**")
			print("Currently ONLY looking at i=0,j=0 element of eta! See line ~ 440")
			for i in range(1):
				for j in range(1):
					# integrate over first axis:
					int_over_x = trapz(imat[:,:,i,j], x=kx, dx=dkx, axis=0)

					# integrate over second axis:
					int_val = trapz(int_over_x, x=ky, dx=dky, axis=-1)

					eta[i,j] = 8 * int_val/ABZ  # note: factor of 8 is because only integrated 1/8 of BZ
				#
			#

			tmp_list_r.append(eta[0,0].real)
			tmp_list_i.append(eta[0,0].imag)

			out_file.write("{0:.2f}\t{1:.6f}\t{2:.6f}\n".format(w, eta[0,0].real, eta[0,0].imag))

			t_end = time.time()
			print("Time elapsed for wavelength {0:d}: {1:.1f}".format(iw, t_end - t_start))

		#for_waves

		#plt.figure()
		#plt.plot(waves, tmp_list_r, color='red', label='real')
		#plt.plot(waves, tmp_list_i, color='blue', label='imag.')
		#plt.xlabel("waves")
		#plt.show()

	
	#if




	if source_type == "dipole_array":
		print("dipole_array source not yet implemented")
	#if

	out_file.close()

	print("Finished Successfully")
	print("Data written to: {0:s}".format(out_file_name))

	'''
	# plot inverse polarizability data [if plot_pol = True]:
	if(plot_pol):
		plt.figure(1)
		plt.plot(waves, pol_real, color="red", label="Real")
		plt.plot(waves, pol_imag, color="blue", label="Imag.")
		if(plot_LS):
			plt.plot(waves, SL_R, color="red", ls="--")
			plt.plot(waves, SL_I, color="blue", ls="--")
		plt.xlabel("$\lambda$ [nm]")
		plt.ylabel(r"$\alpha^{-1}$ [nm$^{-3}$]")
		plt.legend(loc="best")
		plt.show()
	#if

	# calculate extinction efficiency:
	k_array = 2*pi/waves
	LS_array = array(LS_array)
	geo_cross = pi*50.0**2
	ext = 4*pi*k_array*((inv_alpha - LS_array)**(-1)).imag

	if(plot_ext):
		plt.figure(2)
		plt.plot(waves, ext/geo_cross)
		plt.xlabel("$\lambda$ [nm]")
		plt.ylabel("Extinction Efficiency")
		plt.ylim((0,40))
		plt.show()
	#if

	if(save_data):
		out_file_name = "slurpie.out"
		out_file = open(out_file_name, 'w')
		# print header:
		# **__TO_DO__**
		for i in range(len(waves)):
			out_line = "{0:f}\t{1:e}\t{2:e}\t{3:e}\t{4:e}\t{5:f}\n".format(waves[i], pol_real[i], pol_imag[i], SL_R[i], SL_I[i], ext[i])
			out_file.write(out_line)
		#for
		out_file.close()
	#if
	'''


#if
