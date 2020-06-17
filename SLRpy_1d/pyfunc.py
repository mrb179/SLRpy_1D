from __future__ import print_function, division
from numpy import linspace, sqrt, pi, array, abs, meshgrid, zeros, exp, log, real, imag, eye, dot, outer
import numpy.linalg as nl
import matplotlib
#from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy.special import spherical_jn, spherical_yn
from scipy.optimize import minimize
from mpmath import fp, polylog, lerchphi,  re, im
from SFP_dielectrics import Ag_JC
from sys import exit

def f_analytic(x, y, k_parallel, rnp, nind, d, t):
	ic = complex(0,1)
	z = complex(x,y)
	w = (z-1)*((z-ic)**2)*((z+1)**3)/(z+ic)
	wr = w.real
	wi = w.imag
	return w
#

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


def Ag_drude(w, params='charles'):
	'''
	w is complex-valued frequency [eV]
	'''
	if params == "charles":
		# these are the values that Charles always uses for Ag
		einf = 5.77
		wp   = 9.18	# eV
		t    = 0.018 	# [eV]

	elif params == "PRB": 
		#parameters taken (approx.) from Yang et al. PRB 91, 235137 (2015).
		einf = 4
		wp   = 8.9  	# [eV]
		t    = 0.25 	# [eV]

	elif params == 'topolog':
		# parameters taken from Pocock et al. ACS Photonics (2018) -- Topological Chain Paper
		print("\n*** WARNING *** The current parameters for epsilon(w) are for nind = 1.5\n")
		einf = 9.1
		wp   = 9.08332  # [eV]
		t    = 0.071 	# [eV]
	else:
		print("\nvalue of param passed to Ag_drude are not recognized\n")
		exit()
	#if

	eps  = einf - (wp**2)/(w**2 +complex(0,1)*w*t)
	return eps
#def

def calc_pol(wc, r=10, nout=1.0):
	'''
	- wc : complex frequency [eV]
	- r  : NP radius [nm]
	- nout : background refractive index

	NOTE: 
	- Left out eps_0 in the numerator 
	- defined k as just omega*n, NOT as omega*n/c, since this just scales the result 

	'''
	#wc = complex(wr, wi)
	k = wc*nout/(hbar*c)
	eps_Ag = Ag_drude(wc)
	a1 = eval_a1(wc, r, eps_Ag, nout**2)
	alpha = 6*pi*complex(0,1)*(nout**2)*a1/(k**3)
	return alpha
#def



# DEFINITION OF pyfunc.f() for single-NP example. This function is essentially identical to calc_pol() 
# 	defined immediately above.
def f(w_real, w_imag, r_NP, nout):
	wc = complex(w_real, w_imag)
	return calc_pol(wc, r=r_NP, nout=nout)
#def


def calc_pol_vecinput(ww, r, nout):
	'''
	- ww : this is a 1d vector with shape (2,), where ww[0] and ww[1] are the real and imag. parts of frequency, respectively. 
	- r  : NP radius [nm]
	- nout : background refractive index

	NOTE: 
	- Left out eps_0 in the numerator 
	'''
	nout = 1.0
	r = 25.0 
	wc = complex(wr, wi)
	k = wc*nout/(hbar*c)
	#k = wr*nout
	eps_Ag = Ag_drude(wc)
	a1 = eval_a1(wc, r, eps_Ag, nout**2)
	alpha = 6*pi*complex(0,1)*(nout**2)*a1/(k**3)
	return alpha**(-1)
#def

def calc_invpol(wc, r, nout):
	'''
	- ww : this is a 1d vector with shape (2,), where ww[0] and ww[1] are the real and imag. parts of frequency, respectively. 
	- r  : NP radius [nm]
	- nout : background refractive index

	NOTE: 
	- Left out eps_0 in the numerator 
	'''
	k = wc*nout/(hbar*c)
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

def mylerch(z, s, a):
	lt = lerchphi(z, s, a)
	lt_real = fp.re(lt)
	lt_imag = fp.im(lt)
	return complex(lt_real, lt_imag)
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


#################################
#	Interaction Sum		#
#################################

def T_1d_long_pos(wc, d, t, kx, nind):
	'''
	Calculates 1,2  off-diagonal component of the longitudical lattice sum tensor as defined in 
	Eq. 7 (positive sign convention) of the Supporting Information of Pocock et al. ACS Photonics (2018). 
	
	Inputs:
		- wc: complex angular frequency [eV]
		- d : unit cell length [nm]
		- t : gap distance between two NPs within each unit cell [nm]
	'''
	k = wc*nind/(hbar*c)
	prefactor = (2/(4*pi*d**3))
	term_1 = exp(ic*k*t)*mylerch(exp(ic*(k+kx)*d), 3, t/d) + exp(-ic*k*t)*exp(ic*(k-kx)*d)*mylerch(exp(ic*(k-kx)*d), 3, 1-(t/d)) - \
		ic*k*d*( exp(ic*k*t)*mylerch(exp(ic*(k+kx)*d), 2, t/d) + exp(-ic*k*t)*exp(ic*(k-kx)*d)*
		mylerch(exp(ic*(k-kx)*d), 2, 1-(t/d))  )
	return prefactor*term_1
#def 

def T_1d_long_neg(wc, d, t, kx, nind):
	'''
	Calculates 1,2  off-diagonal component of the longitudical lattice sum tensor as defined in 
	Eq. 7 (negitive sign convention) of the Supporting Information of Pocock et al. ACS Photonics (2018). 
	
	Inputs:
		- wc: complex angular frequency [eV]
		- d : unit cell length [nm]
		- t : gap distance between two NPs within each unit cell [nm]
	'''
	k = wc*nind/(hbar*c)
	prefactor = (2/(4*pi*d**3))
	term_1 = exp(ic*k*t)*mylerch(exp(ic*(k-kx)*d), 3, t/d) + exp(-ic*k*t)*exp(ic*(k+kx)*d)*mylerch(exp(ic*(k+kx)*d), 3, 1-(t/d)) \
		- ic*k*d*( exp(ic*k*t)*mylerch(exp(ic*(k-kx)*d), 2, t/d) + exp(-ic*k*t)*exp(ic*(k+kx)*d) \
		*mylerch(exp(ic*(k+kx)*d), 2, 1-(t/d))  )
	return prefactor*term_1
#def 

def T_1d_long(wc, d, t, kx, nind):
	'''
	- Uses the sign of t to decide whethere to use +/- version of interaction sum T.
	- Note that the positive upper(+)/lower(-) sign convention in Eq. 7 of the Pocock
		paper correspond to the evaluation point to the neg./pos. x-direction w.r.t.
		the offset of the sublattice sourcing the interaction in my notes. This is 
		the reason for the conditions on t below. 
	- In the general case for the interaction between sub-lattices i and j, t will be
		the difference between the coordinates of the sublattices in the 0th unit 
		cell.   
	'''
	if t > 0:
		return T_1d_long_neg(wc, d, t, kx, nind)
	else:
		return T_1d_long_pos(wc, d, -1*t, kx, nind)
#def 

def T_1d_trans_pos(wc, d, t, kx, nind):
	'''
	Calculates 1,2  off-diagonal component of the longitudical lattice sum tensor as defined in 
	Eq. 7 (positive sign convention) of the Supporting Information of Pocock et al. ACS Photonics (2018). 
	
	Inputs:
		- wc: complex angular frequency [eV]
		- d : unit cell length [nm]
		- t : gap distance between two NPs within each unit cell [nm]
	'''
	k = wc*nind/(hbar*c)
	prefactor = (-1/(4*pi*d**3))
	term_1 = exp(ic*k*t)*mylerch(exp(ic*(k+kx)*d), 3, t/d) + exp(-ic*k*t)*exp(ic*(k-kx)*d)*mylerch(exp(ic*(k-kx)*d), 3, 1-(t/d)) - \
		ic*k*d*( exp(ic*k*t)*mylerch(exp(ic*(k+kx)*d), 2, t/d) + exp(-ic*k*t)*exp(ic*(k-kx)*d)* mylerch(exp(ic*(k-kx)*d), 2, 1-(t/d))  ) + \
		-((k*d)**2)*(exp(ic*k*t)*mylerch(exp(ic*(k+kx)*d),1,t/d) + exp(-ic*k*t)*exp(ic*(k-kx)*d)*mylerch(exp(ic*(k-kx)*d), 1, 1-(t/d))   ) 
	return prefactor*term_1
#def 

def T_1d_trans_neg(wc, d, t, kx, nind):
	'''
	Calculates 1,2  off-diagonal component of the longitudical lattice sum tensor as defined in 
	Eq. 7 (negitive sign convention) of the Supporting Information of Pocock et al. ACS Photonics (2018). 
	
	Inputs:
		- wc: complex angular frequency [eV]
		- d : unit cell length [nm]
		- t : gap distance between two NPs within each unit cell [nm]
	'''
	k = wc*nind/(hbar*c)
	prefactor = (-1/(4*pi*d**3))
	term_1 = exp(ic*k*t)*mylerch(exp(ic*(k-kx)*d), 3, t/d) + exp(-ic*k*t)*exp(ic*(k+kx)*d)*mylerch(exp(ic*(k+kx)*d), 3, 1-(t/d)) - \
		ic*k*d*( exp(ic*k*t)*mylerch(exp(ic*(k-kx)*d), 2, t/d) + exp(-ic*k*t)*exp(ic*(k+kx)*d)* mylerch(exp(ic*(k+kx)*d), 2, 1-(t/d))  ) + \
		-((k*d)**2)*(exp(ic*k*t)*mylerch(exp(ic*(k-kx)*d),1,t/d) + exp(-ic*k*t)*exp(ic*(k+kx)*d)*mylerch(exp(ic*(k+kx)*d), 1, 1-(t/d))   ) 
	return prefactor*term_1
#def 

def T_1d_trans(wc, d, t, kx, nind):
	'''
	- Uses the sign of t to decide whethere to use +/- version of interaction sum T.
	- see notes above for T_1d_long()
	'''
	if t > 0:
		return T_1d_trans_neg(wc, d, t, kx, nind)
	else:
		return T_1d_trans_pos(wc, d, -1*t, kx, nind)
#def 


def build_T_block(wc, d, t, kx, nind):
	'''
	Input Arguments:
	- wc : complex frequency [eV].
	- d  : intra-sub-lattice site periodicity [nm]
	- t  : gap distance between two NPs within each unit cell [nm]
	- kx : in-plane wave vector along chain (x) axis [nm^(-1)]
	- nind: background refractive index.
	
	NOTE: it is assumed the chain axis is along x-direction
	'''
	T = zeros((3,3), dtype=complex)
	T[0,0] = T_1d_long(wc, d, t, kx, nind) 
	T[1,1] = T_1d_trans(wc, d, t, kx, nind)
	T[2,2] = T[1,1]
	return T
#def


#################################
#	Lattice Sum		#
#################################

def S_diag_1d_long(wc, d, kx, nind):
	'''
	Input Arguments:
	- wc : complex frequency [eV].
	- d  : intra-sub-lattice site periodicity [nm]
	- kx : in-plane wave vector along chain (x) axis [nm^(-1)]
	- nind: background refractive index.

	Notes:
	- form of S_trans taken from [eq. 8] in Cherqui, C., Bourgeois, M., et al. Accounts of Chemical Research (2019) 
	'''
	k = wc*nind/(hbar*c)

	S_long = 2*((4*pi*(d**3))**(-1))*( mypolylog(3, exp(ic*(k-kx)*d)) + mypolylog(3, exp(ic*(k+kx)*d)) -ic*k*d \
		*(mypolylog(2, exp(ic*(k-kx)*d)) + mypolylog(2, exp(ic*(k + kx)*d)) ) )
	return S_long
#def

def S_diag_1d_trans(wc, d, kx, nind):
	'''
	Input Arguments:
	- wc : complex frequency [eV].
	- d  : intra-sub-lattice site periodicity [nm]
	- kx : in-plane wave vector along chain (x) axis [nm^(-1)]
	- nind: background refractive index.

	Notes:
	- form of S_trans taken from [eq. 8] in Cherqui, C., Bourgeois, M., et al. Accounts of Chemical Research (2019) 
	'''
	k = wc*nind/(hbar*c)

	S_trans = -1*((4*pi*(d**3))**(-1))*( mypolylog(3, exp(ic*(k-kx)*d)) + mypolylog(3, exp(ic*(k+kx)*d)) -ic*k*d \
		*(mypolylog(2, exp(ic*(k-kx)*d)) + mypolylog(2, exp(ic*(k + kx)*d)) ) - ((k*d)**2) \
		*(mypolylog(1, exp(ic*(k-kx)*d)) + mypolylog(1, exp(ic*(k + kx)*d)) ) )
	return S_trans
#def


def build_S_block(wc, d, kx, nind):
	'''
	Input Arguments:
	- wc : complex frequency [eV].
	- d  : intra-sub-lattice site periodicity [nm]
	- kx : in-plane wave vector along chain (x) axis [nm^(-1)]
	- nind: background refractive index.
	
	NOTE: it is assumed the chain axis is along x-direction
	'''
	S = zeros((3,3), dtype=complex)
	S[0,0] = S_diag_1d_long(wc, d, kx, nind) 
	S[1,1] = S_diag_1d_trans(wc, d, kx, nind) 
	S[2,2] = S[1,1]
	return S
#def

 
def GD_single(reval, rsource, wc, nind):
	I = eye(3)
	r_vec = reval - rsource
	R = sqrt(dot(r_vec, r_vec))
	rhat_rhat = (R**(-2))*outer(r_vec, r_vec)
	#kc = wc*nind/(hbar*c)
	kc = 2*pi*nind/wc
	prefactor = exp(ic*kc*R)/(4*pi*R)
	term_1 = (1 + ic/(kc*R) - 1/(kc*R)**2)	
	term_2 = (-1 - 3*ic/(kc*R) + 3/(kc*R)**2)
	G = ((kc**2)/(nind**2))*prefactor*( term_1*I + term_2*rhat_rhat)
	return G*nind**2
	# Note: extra nind**2 needed to reproduce Pocock results (understood)
#def


# THIS WILL BE THE FUNCYION TO REPRODUCE THE ACS PHOTONICS PAPER:
#
# NOTE:
# 	- GRPF_python/main_dispersion.py is set up to run with this as pyfunc.f() with the analysis_parameters_TEMPLATE.m templated input files. 
def f_ACS_Photon_dispersion(w_real, w_imag, k_parallel, rnp, nind, d, t):
	'''
	This function constructs the two-NP unit cell G_SLR matrix and returns its determinant 
	evaluated at the complex angular frequency w = w_real + i*w_imag.
	'''

	mode_dir = "long"
	
	wc = complex(w_real, w_imag)

	if mode_dir == "long":
		print("** Longitudinal Polarized Mode **")
		S = zeros((2,2), dtype=complex)
		S[0,0] = S_diag_1d_long(wc, d, k_parallel, nind)
		S[1,1] = S[0,0]
		S[0,1] = T_1d_long_pos(wc, d, t, k_parallel, nind) 
		S[1,0] = T_1d_long_neg(wc, d, t, k_parallel, nind) 
	else:
		print("mode direction not yet implemented")
		exit()
	
	# calculate the NP polarizability tensor:
	inv_alpha_mie = calc_invpol(wc, rnp, nind)
	
	I = eye(2, dtype=complex)
	A = inv_alpha_mie*I - S
	
	R = A[0,0]*A[1,1] - A[0,1]*A[1,0]

	return R
#


