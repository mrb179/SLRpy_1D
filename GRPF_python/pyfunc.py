from __future__ import print_function, division
from numpy import linspace, sqrt, pi, array, abs, meshgrid, zeros, exp, log, real, imag
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import spherical_jn, spherical_yn
from scipy.optimize import minimize
from mpmath import fp, polylog, re, im
from SFP_dielectrics import Ag_JC
from sys import exit

######################
##    CONSTANTS:	##
######################

# complex unit:
ic = complex(0,1)

# speed of light [nm/s]:
c = 2.99792458*10**(17)

# reduced Plank constant [eV*s]:
hbar = 6.582119569*10**(-16)



##########################
##    FUNCTION DEFs:	##
##########################

def psi_n(n,z):
	return z*spherical_jn(n,z, derivative=False)

def psi_n_deriv(n,z):
	return spherical_jn(n,z) + z*spherical_jn(n,z, derivative=True)

def chi_n(n,z):
	return z*(spherical_jn(n,z) + complex(0,1)*spherical_yn(n,z))

def chi_n_deriv(n,z):
	return ( (spherical_jn(n,z) + complex(0,1)*spherical_yn(n,z)) + z*(spherical_jn(n,z, derivative=True) + complex(0,1)*spherical_yn(n,z, derivative=True)) )

def eval_a1(w, R, eps_in, eps_out):
	# returns a1 Mie scattering coefficient
	# w: complex-valued frequency (eV)
	# R: sphere radius [nm]
	# eps_in: complex valued dielectric function inside NP
	hbar = 6.582119*10**(-16) # [eV*s]
	c    = 2.99792458*10**(8) #[m/s]
	x = (R*10**(-9))*w*sqrt(eps_out)/(hbar*c)
	m = sqrt(eps_in)/sqrt(eps_out)
	num   = m*psi_n(1,m*x)*psi_n_deriv(1,x) - psi_n(1,x)*psi_n_deriv(1,m*x)
	denom = m*psi_n(1,m*x)*chi_n_deriv(1,x) - chi_n(1,x)*psi_n_deriv(1,m*x)
	return (num/denom) 
# def

def Ag_drude(w):
	# parameters taken (approx.) from Yang et al. PRB 91, 235137 (2015).
	# w is allowed to be complex-valued
	einf = 4
	wp   = 8.9  # [eV]
	t    = 0.25 # [eV]
	eps  = einf - (wp**2)/(w**2 +complex(0,1)*w*t)
	return eps
#def

def calc_pol(wr, wi, r=10, nout=1.0):
	'''
	- wr : real part of frequency [eV]
	- wi : imag part of frequency [eV]
	- r  : NP radius [nm]
	- nout : background refractive index

	NOTE: 
	- Left out eps_0 in the numerator 
	- defined k as just omega*n, NOT as omega*n/c, since this just scales the result 

	'''
	wc = complex(wr, wi)
	k = wc*nout
	#k = wr*nout
	eps_Ag = Ag_drude(wc)
	a1 = eval_a1(wc, r, eps_Ag, nout**2)
	alpha = 6*pi*complex(0,1)*(nout**2)*a1/(k**3)
	return abs(alpha**(-1))
#def

def calc_pol_vecinput(ww, r, nout):
	'''
	- ww : this is a 1d vector with shape (2,), where ww[0] and ww[1] are the real and imag. parts of frequency, respectively. 
	- r  : NP radius [nm]
	- nout : background refractive index

	NOTE: 
	- Left out eps_0 in the numerator 
	'''
	wc = complex(ww[0], ww[1])
	k = wc*nout/(hbar*c)
	#k = wr*nout
	eps_Ag = Ag_drude(wc)
	a1 = eval_a1(wc, r, eps_Ag, nout**2)
	alpha = 6*pi*complex(0,1)*(nout**2)*a1/(k**3)
	#return abs(alpha**(-1))	# this was from when I was minimizing inv alpha
	return alpha**(-1)
#def

def mypolylog(s,z):
	PL = polylog(s,z)
	PL_real = fp.re(PL)
	PL_imag = fp.im(PL)
	return complex(PL_real, PL_imag)
#def

def S_1d_trans(ww, a, kp, nout):
	'''
	Input Arguments:
	- ww : this is a 1d vector with shape (2,), where ww[0] and ww[1] are the real and imag. parts of frequency [eV], respectively.
	- a  : site periodicity [nm]
	- kp : in-plane wave vector along chain axis [nm^(-1)]

	Notes:
	- form of S_trans taken from [eq. 8] in Cherqui, C., Bourgeois, M., et al. Accounts of Chemical Research (2019) 
	'''
	wc = complex(ww[0], ww[1])

	k = wc*nout/(hbar*c)
	#k = ww[0]*nout/(hbar*c)

	#print("check 1")
	term_1 = (a**2)*(k**2)*log( (1 - exp(ic*a*(k-kp))) * (1 - exp(ic*a*(k+kp))) )

	#print("check 2")
	term_2 = ic*a*k*mypolylog(2, exp(ic*a*(k-kp)))

	#print("check 3")
	term_3 = ic*a*k*mypolylog(2, exp(ic*a*(k+kp)))

	#print("check 4")
	term_4 = mypolylog(3, exp(ic*a*(k-kp)) )
	#print(term_4)

	#print("check 5")
	term_5 = mypolylog(3, exp(ic*a*(k+kp)) )
	#print(term_5)

	#print("check 6")
	St = ((4*pi*(a**3))**(-1))*( term_1 - term_2 - term_3 - term_4 - term_5 )

	return St

#def 

def S_1d_trans_ALT(ww, a, kp, nout):
	'''
	Input Arguments:
	- ww : this is a 1d vector with shape (2,), where ww[0] and ww[1] are the real and imag. parts of frequency [eV], respectively.
	- a  : site periodicity [nm]
	- kp : in-plane wave vector along chain axis [nm^(-1)]

	Notes:
	- form of S_trans taken from [eq. 8] in Cherqui, C., Bourgeois, M., et al. Accounts of Chemical Research (2019) 
	'''
	wc = complex(ww[0], ww[1])

	k = wc*nout/(hbar*c)

	term_1 = (a**2)*(k**2)*mypolylog(1, exp(ic*a*(k-kp)) )

	term_2 = (a**2)*(k**2)*mypolylog(1, exp(ic*a*(k+kp)) )

	term_3 = ic*a*k*mypolylog(2, exp(ic*a*(k-kp)))

	term_4 = ic*a*k*mypolylog(2, exp(ic*a*(k+kp)))

	term_5 = mypolylog(3, exp(ic*a*(k-kp)) )

	term_6 = mypolylog(3, exp(ic*a*(k+kp)) )

	St = -1*((4*pi*(a**3))**(-1))*( term_6 + term_5 - term_3 - term_4 - term_1 - term_2 )

	return St

#def 

def empty_latt_1d(kp, aa, nout, bi):
	'''
	Input Arguments:
	- kp   : in-plane wave vector [nm^(-1)] (real-valued)
	- aa   : lattice periodicity [nm]
	- nout : background refractive index (real-valued)
	- bi   : Band Index

	Returns:
	Ek : energy of band [eV] at in-plane wave vector kp.  
	'''
	Ek = (hbar*c/nout)*(kp + (2*bi*pi/aa))
	return Ek
#def 

def dispersion_1d(kvals, aa, nout, mbi):
	'''
	Input Arguments:
	- kvals: (array-like) vector of  
	- aa   : lattice periodicity [nm]
	- nout : background refractive index (real-valued)
	- mbi  : max. band index -- calculates empty lattice dispersion from -bi to bi 

	Returns:
	Ek :   
	'''
	# create container to hold dispersion of bands [-bi, bi]
	band_dispersion_list = []
	# loop over bands:
	for bi in range(-mbi, mbi+1):
		band_dispersion = array([empty_latt_1d(kp, aa, nout, bi) for kp in kvals])
		band_dispersion_list.append(bi*band_dispersion)
	#for
	return band_dispersion_list
#def

def g_SLR_inv(ww, a, kp, r, nout):
	'''
	- ww : this is a 1d vector with shape (2,), where ww[0] and ww[1] are the real and imag. parts of frequency, respectively. 
	- a  : site periodicity [nm]
	- kp : in-plane wave vector along chain axis [nm^(-1)]
	- r  : NP radius [nm]
	- nout : background refractive index

	NOTE:
	- Currently *** returning absolute value *** of g_SLR_inv in order to perform minimization
	'''

	# this returns the (complex-valued) inverse of the NP polarizability
	ainv = calc_pol_vecinput(ww, r, nout)
	
	# this is the 1d lattice sum (single sub-lattice)
	lsum = S_1d_trans_ALT(ww, a, kp, nout)

	gSLR = ainv - lsum  

	return abs(gSLR)
	'''
	NOTE: Not sure what is going on here... grid search plot doesn't look correct
	'''
#def





##########################
##    		CODE:		##
##########################

if __name__ == '__main__':

	'''
	Run Types: 
	- "grid_search": brute force evaluation over complex frequency grid 
	- "minimize_NM": Nelder-Mead minimization 
	- "lattice_1d" : 1d periodic lattice
	- "dispersion_1d_test" : testing empty lattice dispersion
	'''

	nout = 1.00
	R_NP = 25
	run_type = "grid_search"




	## Real-Frequency Spectrum ##

	waves = linspace(300, 500, num=151)
	k     = 2*pi*nout/waves
	waves_eV = 1240/waves

	#eps_Ag = Ag_JC(waves)
	eps_Ag = Ag_drude(waves_eV)

	a1 = array([eval_a1(waves_eV[i], R_NP, eps_Ag[i], nout) for i in range(len(waves))])
	ext = (6*pi/(k*k))*a1.real/(pi*R_NP**2)



	if run_type == "grid_search":

		print("\nrun_type: grid_search")

		w_real = linspace(3.0, 4.0, 101)
		w_imag = linspace(-0.4, 0.2, 101)
		#w_real = linspace(2.8, 3.4, 101)
		#w_imag = linspace(-0.3, 0, 101)
		Wr, Wi = meshgrid(w_real, w_imag)

		print("")
		print("Array Dimensions:")
		print(Wr.shape)
		print(Wi.shape)
		print("")

		# Commenting out since this had issues, brute forcing for now:
		#alpha_mag = calc_pol(Wr, Wi)

		holding_pen = zeros((len(w_imag), len(w_real)), dtype=float)
		for i, ww_r in enumerate(w_real):
			for j, ww_i in enumerate(w_imag):
				holding_pen[j,i] = calc_pol(ww_r, ww_i, r=R_NP, nout=nout )
			#for
		#for
		print(holding_pen)
		print("")
		print(holding_pen.shape)


		## PLOTTING ##
		font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

		matplotlib.rc('font', **font)

		# plot surface plot of |alpha(omega_complex)|:
		fig = plt.figure(1)
		ax = plt.axes(projection="3d")
		ax.plot_surface(Wr, Wi, holding_pen, cmap='viridis', edgecolor='none')
		#ax.set_xlabel("$Re(\Omega_{m})$")
		#ax.set_ylabel('$Im(\Omega_{m})$')



		# plot single LSP spectrum (all real frequencies):
		evplot = 1240./waves
		plt.figure(2)
		plt.plot(evplot, ext, color=[0,0.2,0.6])
		plt.ylim((0,1.1*max(ext)))
		plt.xlim((min(evplot), max(evplot)))
		plt.xlabel("Energy (eV)")
		#plt.ylabel("")
		plt.show()

	# End_grid_search 



	elif run_type == "minimize_NM":

		## PARAMETERS ##

		# initial guess:
		w0 = array([3.4, -0.02])

		# solution tolerance:
		my_tol = 10**(-3)

		max_iterations = 200

		# print convergence details:
		display_conv = True



		print("\nrun_type: minimize_NM\n")

		min_result = minimize(calc_pol_vecinput, x0=w0, args=(R_NP, nout), method='Nelder-Mead', tol=my_tol, options={'maxiter': max_iterations, 'disp': display_conv, 'return_all': False, 'initial_simplex': None, 'xatol': 0.0001, 'fatol': 0.0001})
		#min_result = minimize(calc_pol_vecinput, x0=w0, method='Nelder-Mead', tol=my_tol, options={'maxiter': max_iterations, 'disp': display_conv, 'return_all': False, 'initial_simplex': None, 'xatol': 0.0001, 'fatol': 0.0001})

		print("")
		print("__Minimization Results:__\n")
		if min_result.success:
			print("Minimization finished sucessfully:")
			print("\tomega (real) [eV] = {0:.3f}".format(min_result.x[0]))
			print("\tomega (imag) [eV] = {0:.3f}".format(min_result.x[1]))
		else:
			print("** minimization failed **\n")
			print(min_result.message)
		print("")
		print("Number of iterations performed: {0:d}".format(min_result.nit))
		print("")

	# end minimize_NM


	elif run_type == "lattice_1d":

		# lattice periodicity [nm]:
		aa = 470
		kpp = pi/(2*aa)
		# R_NP and nout variables set above at beginning of Code block

		E_num  = 151
		Er_min = 2.63
		Er_max = 2.64
		Ei_min = -0.015
		Ei_max = 0.005
		Evals_real = linspace(Er_min, Er_max, num=E_num)
		Evals_imag = linspace(Ei_min, Ei_max, num=E_num)

		w_vec = array([array([Evals_real[i], Evals_imag[i]]) for i in range(E_num)])

		#print("Real parts:")
		##print(real(w_vec))
		#print("")
		#print("Imag. parts:")
		#print(imag(w_vec))

		print("\nrun_type: lattice_1d\n")

		# Calculate polarizability and lattice spectra (Real Frequencies):
		LS_trans = array([S_1d_trans(wi, aa, kpp, nout) for wi in w_vec])
		alpha_inv = array([(calc_pol_vecinput(wi, R_NP, nout)) for wi in w_vec])

		# Brute force grid search of domain of complex frequencies (Gamma point only)
		w_real = linspace(2.8, 3.5, num=101)
		w_imag = linspace(-1, 0.005, num=101)
		#w_real = linspace(2.6, 2.65, num=201)
		#w_imag = linspace(-0.05, 0.01, num=201)
		Wr, Wi = meshgrid(w_real, w_imag)
		holding_pen = zeros((len(w_imag), len(w_real)), dtype=float)
		for i, ww_r in enumerate(w_real):
			for j, ww_i in enumerate(w_imag):
				omega_c = array([ww_r, ww_i])
				holding_pen[j,i] = g_SLR_inv(omega_c, aa, kpp, R_NP, nout)
				#holding_pen[j,i] = -1*S_1d_trans(omega_c, aa, kpp, nout)
			#for
		#for

		# Plot real frequency spectra:
		fig = plt.figure(1)
		plt.plot(Evals_real, -1*real(LS_trans), color="red")
		plt.plot(Evals_real, -1*imag(LS_trans), color="blue")
		plt.plot(Evals_real, real(alpha_inv), color="black")
		#plt.xlim((1.0,3.5))
		#plt.ylim((-0.5*10**(-7),7*10**(-7)))

		# plot grid search data:
		fig = plt.figure(2)
		ax = plt.axes(projection="3d")
		ax.plot_surface(Wr, Wi, holding_pen, cmap='inferno', edgecolor='none')

		enable_num_min = False
		if enable_num_min:
			# initial value:
			#w0 = array([3.3, -0.5])
			w0 = array([2.6, -0.01])
			# solution tolerance:
			my_tol = 10**(-4)
			max_iterations = 200
			# print convergence details:
			display_conv = True
			min_result = minimize(g_SLR_inv, x0=w0, args=(aa, kpp, R_NP, nout), method='Nelder-Mead', tol=my_tol, options={'maxiter': max_iterations, 'disp': display_conv, 'return_all': False, 'initial_simplex': None, 'xatol': 0.0001, 'fatol': 0.0001})

			print("")
			print("__Minimization Results:__\n")
			if min_result.success:
				print("Minimization finished sucessfully:")
				print("\tomega (real) [eV] = {0:.5f}".format(min_result.x[0]))
				print("\tomega (imag) [eV] = {0:.5f}".format(min_result.x[1]))
			else:
				print("** minimization failed **\n")
				print(min_result.message)
			print("")
			print("Number of iterations performed: {0:d}".format(min_result.nit))
			print("")
			# add initial guess w0 to grid search surface plot:
			ax.scatter(w0[0], w0[1], zs=g_SLR_inv([w0[0], w0[1]], aa, kpp, R_NP, nout), zdir='z', s=20, c=[0,0,1], depthshade=True)
			# add minimized result to grid search surface plot:
			ax.scatter(min_result.x[0], min_result.x[1], zs=g_SLR_inv([min_result.x[0], min_result.x[1]], aa, kpp, R_NP, nout), zdir='z', s=20, c=[1,0,0], depthshade=True)
		#if

		plt.xlabel("$\Re[{\omega}]$")
		plt.ylabel("$\Im[{\omega}$]")
		#plt.zlabel("det($G_{SLR}(\omega)$)")
		plt.show()

	# end lattice_1d



	elif run_type == "dispersion_1d_test":

		aa = 470 

		print("\n** run_type = dispersion_1d_test **\n")	

		kvals = linspace(0,pi/aa, num=41)

		# Generate empty-lattice dispersion:
		band_dispersion = dispersion_1d(kvals, aa, nout, 1)

		# Generate complex CDM dispersion bands:
		w0_photon = array([2.64, -0.006])
		w0_LSP    = array([3.3, -0.7]) 
		# solution tolerance:
		my_tol = 10**(-5)
		max_iterations = 200
		# print convergence details:
		display_conv = False

		# loop over k_parallel to build dispersion diagram:
		omega_opt_photon = []
		omega_opt_LSP = []
		for kpp in kvals:
			if len(omega_opt_LSP) > 0:
				# use previous optimized value as initial guess for current iteration
				w0_photon = omega_opt_photon[-1]
				w0_LSP = omega_opt_LSP[-1]
			min_photon = minimize(g_SLR_inv, x0=w0_photon, args=(aa, kpp, R_NP, nout), method='Nelder-Mead', tol=my_tol, options={'maxiter': max_iterations, 'disp': display_conv, 'return_all': False, 'initial_simplex': None, 'xatol': 0.0001, 'fatol': 0.0001})
			min_LSP = minimize(g_SLR_inv, x0=w0_LSP, args=(aa, kpp, R_NP, nout), method='Nelder-Mead', tol=my_tol, options={'maxiter': max_iterations, 'disp': display_conv, 'return_all': False, 'initial_simplex': None, 'xatol': 0.0001, 'fatol': 0.0001})
			if min_photon.success:
				omega_opt_photon.append(min_photon.x)
			else:
				omega_opt_photon.append(array([0,0]))
			if min_LSP.success:
				omega_opt_LSP.append(min_LSP.x)
			else:
				omega_opt_LSP.append(array([0,0]))

		# seperate real and imaginary frequency components:
		omega_photon_real = [omega_opt_photon[i][0] for i in range(len(omega_opt_photon))]
		omega_photon_imag = [omega_opt_photon[i][1] for i in range(len(omega_opt_photon))]
		omega_LSP_real = [omega_opt_LSP[i][0] for i in range(len(omega_opt_LSP))]
		omega_LSP_imag = [omega_opt_LSP[i][1] for i in range(len(omega_opt_LSP))]

		plt.figure(1)
		for i in range(len(band_dispersion)):
			plt.plot(kvals, band_dispersion[i], label="Band Index = {0:d}".format(i-1))
		plt.scatter(kvals, omega_photon_real, color=[0,0,1])
		plt.scatter(kvals, omega_LSP_real, color=[1,0,0])
		plt.xlabel("k$_{||}$")
		plt.ylabel("Energy [eV]")
		plt.ylim((2,3.8))
		plt.legend(loc='best')

		plt.figure(2)
		ax = plt.axes(projection="3d")
		ax.plot(kvals, omega_photon_real, omega_photon_imag, color=[0,0,1])
		ax.plot(kvals, omega_LSP_real, omega_LSP_imag, color=[1,0,0])
		plt.xlabel("k$_{||}$")
		plt.ylabel("$\Re[\omega]$ [eV]")
		plt.ylim((2,3.8))

		plt.show()

	# end dispersion_1d_test


	else:

		print("\n** run_type not recognized **\n")





