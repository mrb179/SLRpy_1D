from __future__ import division, print_function
from numpy import pi, sqrt, real, imag,  sin, cos, array, dot, linspace, matmul, cross, exp, conj, zeros, identity, loadtxt, arcsin, arccos, eye, vstack, savetxt
from scipy import interpolate
from numpy.linalg import inv
from pyfunc import GD_single
import pylab as plt
from sys import exit
#import multiprocessing as mp

'''
MODIFICATIONS:

- adding in ability to have n != 1.00 -> Finished: see (quantitative) comparison with FDTD in:
	/Users/marcbourgeois/Desktop/Work/QD-SLR/FDTD_SPECTRA   (run plot_ext.py)
- still need to go back to compare with finite lattice (c++) code 
- later include parallelization
'''

if __name__ == "__main__":

	#################################################
	#			INPUTS:			#
	#################################################
	
	# Vector position of NP [nm]:
	rp = [0,0,0]

	# Vector position of dipole:
	rd = [70,0,0]

	# NP polarizability paths:
	long_polarizability_path   = "polarizability_Ag_80nmD_n1p0.out"
	short_polarizability_path  = "polarizability_Ag_80nmD_n1p0.out"

	# refractive index of surrounding medium:
	nind = 1.00

	# set wavelength scan parameters: [nm]: 
	wmin = 325.
	wmax = 600.
	wnum = 201

	# unit vector specifying dipole orientation:
	nd = [0,1,0]

	# Toggles:
	PF_lattice = False
	plot_FDTD  = False
	plot_A     = False
	plot_IBZ   = False
	IBZ_file   = "slrpy_ibz_450nmP_AgNP_80nmD_40nmH_n1.47_75x75kpoints.out"

	# Enable various toggle options:
	single_file_name = "PF_single_Ag_80nmD_n1p0_30nm_trans_alt.out"
	save_single      = True

	###################################################
	###################################################

	#################
	# 	POL_LONG	#
	#################
	'''
	polarizability along the long axis
	'''

	# path to file with polarizability data:
	polar_file = long_polarizability_path

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
	#pol_arr = array([complex(pol_data[i,2]*scale_pol, pol_data[i,3]*scale_pol) for i in range(len(pol_data))])
	pol_real_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,2]*scale_pol)
	pol_imag_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,3]*scale_pol)
	pol_real = pol_real_func(waves)
	pol_imag = pol_imag_func(waves)
	alpha_long = array([complex(pol_real[w], pol_imag[w]) for w in range(len(waves))])


	#################
	# 	POL_SHORT	#
	#################
	'''
	polarizability along the short axis
	'''

	# path to file with polarizability data:
	polar_file = short_polarizability_path

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

	'''
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
	'''
	# interpolate polarizability data:
	# create complex polarizability array:
	# 
	# apply scale factor to convert polarizability file data from [um^3] to [nm^3]:
	scale_pol = 10**9
	#pol_arr = array([complex(pol_data[i,2]*scale_pol, pol_data[i,3]*scale_pol) for i in range(len(pol_data))])
	pol_real_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,2]*scale_pol)
	pol_imag_func = interpolate.InterpolatedUnivariateSpline(pol_data[:,0], pol_data[:,3]*scale_pol)
	pol_real = pol_real_func(waves)
	pol_imag = pol_imag_func(waves)
	alpha_short = array([complex(pol_real[w], pol_imag[w]) for w in range(len(waves))])






	#########################################
	#			PURCELL FACTOR CALC. 		#
	#########################################

	print("")
	print("NOTE:")
	print("\tUnlike my other codes, this PF script uses SI(ish) units")
	print("\tso polarizability is NOT normalized by 4*pi.")
	print("")

	#####################################
	# load BZI data and interpolate:	#
	#####################################

	rp = array(rp)	# NP position
	rd = array(rd)	# emitter position
	nd = array(nd)	# emitter orientation

	# ready output file:
	#out_file = open(out_file_name, 'w')


	# container for single NP Purcell factors:
	pfc_single    = []
	
	# container to hold self-propagator values:
	A_list = []	

	if PF_lattice:
		# ready files for lattice Purcell factor calculations:
		ibz_dat = loadtxt(IBZ_file)
		print(ibz_dat[:,0])
		if ibz_dat[0,0] > min(waves):
			print("min wave < min ibz wave. aborting")
			exit()
		elif ibz_dat[-1,0] < max(waves):
			print(ibz_dat[0,-1])
			print("max wave > max ibz wave. aborting")
			exit()
		ibz_real_func = interpolate.InterpolatedUnivariateSpline(ibz_dat[:,0], ibz_dat[:,1]/(4*pi))
		ibz_imag_func = interpolate.InterpolatedUnivariateSpline(ibz_dat[:,0], ibz_dat[:,2]/(4*pi))
		ibz_real = ibz_real_func(waves)
		ibz_imag = ibz_imag_func(waves)
		ibz_complex = array([complex(ibz_real[i], ibz_imag[i]) for i in range(wnum)])
		
		# container for lattice Purcell factors:
		pf_lattice = []
	#if

	for iw, w in enumerate(waves):

		print("\n__NOTE__\nThere are still questions related to factors of epsilon (only trust for n = 1 for now.)\n")

		# wave vector in medium
		k = 2*pi*nind/w
	
		# calculate Green Dyad from rd to rp:
		# Gs is the script_G in my Overleaf notes
		Gs = GD_single(rd, rp, w, nind)/(nind**2) # this matches currently the definition form sum_dyadic_test

		# create polarizability tensor at current wavelength:
		alphaT = array([alpha_long[iw], alpha_long[iw], alpha_short[iw]])*eye(3)

		# calculate single-NP self-propagator:
		A = (nind**4)*matmul(nd.T, matmul( Gs, matmul(alphaT ,matmul(Gs,nd))))

		# calculate single-NP Purcell factor:
		prefactor = 6*pi/(nind*k**3)
		PF_single = prefactor*A
		
		# log single NP PF:
		pfc_single.append(PF_single.imag) # <---- When I changed this to PF_single.imag + PF_single.real, seemed to fix sign issue
		# log single-NP self-propagator:
		A_list.append(A)
		
		if PF_lattice:
			PF = prefactor*ibz_complex[iw]
			pf_lattice.append(PF.imag)
		#if
		
	#for
	A_list = array(A_list)

	plt.figure()
	#plt.plot(1240/waves, pfc, color="blue", label="CDA")
	plt.plot(1240/waves, pfc_single, color="black", label="single")
	#plt.xlim((min(waves), max(waves)))
	plt.ylim((0,max(pfc_single)*1.2))

	if plot_FDTD:
		fdtd_dat = loadtxt("purcell_Ag_disk_40nmR_40nmH_Palik_n1.0_50nmX_n1.0_fdtd.txt", skiprows=3, delimiter=',')
		plt.plot(fdtd_dat[:,0], fdtd_dat[:,1], color="black", label="FDTD")
	#if
	plt.legend(loc="best")
	plt.xlabel("Energy (eV)")
	plt.ylabel("Purcell Factor")
	#plt.savefig("PF_Ravindra_CDA.pdf")

	if plot_A:
		plt.figure(2)
		plt.plot(1240/waves, A_list.real, color='red', label="Real")
		plt.plot(1240/waves, A_list.imag, color='blue', label="Imag.")
		plt.xlabel("Wavelength (nm)")
		plt.xlim((min(1240/waves), max(1240/waves)))
		#plt.xlim((1.6,2.4))
		plt.ylabel("A")
		if plot_IBZ:
			ibz_dat = loadtxt(IBZ_file)
			plt.scatter(1240/ibz_dat[:,0], ibz_dat[:,1]/(4*pi), color='grey')
			plt.scatter(1240/ibz_dat[:,0], ibz_dat[:,2]/(4*pi), color='green')
	#
	#plt.savefig("GF_singleNP_SLR.pdf", transparent=True)
	plt.show()

	# save single NP PF data if enabled:
	if save_single:
		head_str = "# 0: wave [nm]\t1: PF\t2: Re(A)\t3: Im(A)"
		out_dat  = vstack((waves, pfc_single,real(A_list), imag(A_list)))
		savetxt(single_file_name, out_dat.T, header=head_str)
		print("Data written to: {0:s}".format(single_file_name))
	#	

	
	'''
	#############################################
	#			CDA PURCELL FACTOR CALC. 		#
	#############################################

	omegas = 1240./waves

	def polarizability_v2(E):
		E0 = 3.3 
		gm = 0.08
		A  = 5.8*10**(-30)
		alpha =  A/(E0**2 - E**2 - complex(0,1)*E*gm )
		return alpha
	#def

	# params above chosen for 10 D dipole moment at 590 nm
	eps_0 = (8.8541878128*10**(-12))
	beta = array([polarizability_v2(oi)*10**(9)/eps_0 for oi in omegas])

	# Scale emitter polarizability beta to have same units as NP polarizability
	#	alpha -> beta is in [SI] units, so first divide by vacuum permittivity
	# 	and then convert m^3 -> nm^3
	#eps_0 = (8.8541878128*10**(-12))
	#beta_scaled = beta *(10**(27))/eps_0
	#print(beta_scaled.real)
	plot_D = False

	if plot_D:
		print(beta.real)
		#scale = 3.33*10**(-30)

		plt.figure()
		plt.plot(waves, beta.real, color="green")
		plt.plot(waves, beta.imag, color="blue")
		plt.xlabel("wavelength (nm)")
		plt.ylabel("$beta$")
		plt.show()
	#if

	ibeta = 1./beta
	print(A_list)
	V = (1./(ibeta - (A_list))).imag #/ (beta.imag)
	print("")
	print(V)

	plt.figure()
	plt.plot(waves, V)
	plt.show()

	if save_data:
		out_dat = vstack((waves, pfc))
		savetxt(out_file_name, out_dat.T)
		print("Data written to: {0:s}".format(out_file_name))
	#
	'''

	#out_file.close()

	print("Finished Successfully")

#if
