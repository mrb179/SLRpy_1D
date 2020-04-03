from __future__ import print_function, division
from numpy import linspace, sqrt, pi, array, abs, meshgrid, zeros, exp, log, real, imag, eye
from numpy.linalg import det
import matplotlib
from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy.special import spherical_jn, spherical_yn
from scipy.optimize import minimize
from mpmath import fp, polylog, lerchphi,  re, im
from SFP_dielectrics import Ag_JC
from sys import exit
from pyfunc import *
from numpy import array, zeros, eye, sqrt, outer, exp, pi, dot 

hbar = 6.582119569*10**(-16) 	# [eV s]
c    = 2.99792458*10**(17)	# [nm/s]
ic   = complex(0,1) 

def sg(reval, iN, a, wc, kx, nind):
	I = eye(3)
	r_vec = reval - array([iN*a, 0, 0])
	R = sqrt(dot(r_vec, r_vec))
	rhat_rhat = (R**(-2))*outer(r_vec, r_vec)
	kc = wc*nind/(hbar*c)
	prefactor = exp(ic*kc*R)/(4*pi*R)
	term_1 = (1 + ic/(kc*R) - 1/(kc*R)**2)	
	term_2 = (-1 - 3*ic/(kc*R) + 3/(kc*R)**2)
	G = ((kc**2)/(nind**2))*prefactor*( term_1*I + term_2*rhat_rhat)*exp(ic*kx*(iN*a))
	return G*nind**2
	# Note: extra nind**2 needed to reproduce Pocock results (understood)
#def

def S_1D_x_brute(N, a, wc, kx, nind):
	S = zeros((3,3), dtype=complex)
	for iN in range(-N, N+1):
		if iN != 0:
			S += sg(array([0,0,0]), iN, a, wc, kx, nind)
	return S
#def

def T_1D_x_brute(reval, N, a, wc, kx, nind):
	T = zeros((3,3), dtype=complex)
	for iN in range(-N,N+1):
		T += sg(reval, iN, a, wc, kx, nind)
	#for
	return T
#def
	

if __name__ == "__main__":

	# Units [eV]:
	wave_min = 1.5
	wave_max = 3.5
	wave_num = 101
	
	# Lattice Params:
	test_type = "T_trans"
	r_np = 20
	nind = 1.3 
	d = 470
	t = 230
	kx = (pi/d)*0.7
	'''
	Notes:
	- d is the unit cell periodicity [nm]
	- t is the offset along the x-direction between sub-lattices [nm]
	- test_type can be "S_trans" or "T_trans" 
	'''

	# number of terms in brute force eval = 2N + 1:
	N = 201

	# omega values [eV]:
	omegas_real = linspace(wave_min, wave_max, num=wave_num)
	omegas_comp = array([complex(wr, 0) for wr in omegas_real]) 

	if test_type == "S_trans":
		print("\nTesting S (trans)\n")	
		#LS = array([S_diag_1d_long(wc, d, kx, nind) for wc in omegas_comp])
		LS = array([S_diag_1d_trans(wc, d, kx, nind) for wc in omegas_comp])
		LS_BF = array([S_1D_x_brute(N, d, wc, kx, nind) for wc in omegas_comp])
	
		inv_alpha = array([calc_invpol(wc, r_np, nind) for wc in omegas_comp])

		#	ev 1/(s*eV) *nm^3* s/nm * 1/(nm**2)
		# NOTE: extra factor of nind might be to compensate definition of alpha
		#ext = (waves/(nind*c*hbar))*imag(inv_alpha**(-1))/(pi*r_np**2)	
		#plt.figure()
		#plt.plot(1240/waves, ext, color="blue")

		waves = 1240/omegas_real
		print(LS_BF.shape)	

		plt.figure()
		plt.plot(waves, real(LS), color="red", label="Re. Polylog")
		plt.plot(waves, imag(LS), color="blue", label="Im. Polylog")
		plt.plot(waves, real(LS_BF[:,1,1]), color="black", ls="--", label="Re. Brute")
		plt.plot(waves, imag(LS_BF[:,1,1]), color="green", ls="--", label="Im. Brute")
		#plt.plot(waves, real(inv_alpha), color="orange")
		#plt.plot(waves, imag(inv_alpha), color="purple")
		plt.legend(loc="best")
		plt.show()
	#END: if test_type == "S_trans"

	elif test_type == "S_long":
		print("\nTesting S (long)\n")	
		#LS = array([S_diag_1d_long(wc, d, kx, nind) for wc in omegas_comp])
		LS = array([S_diag_1d_long(wc, d, kx, nind) for wc in omegas_comp])
		LS_BF = array([S_1D_x_brute(N, d, wc, kx, nind) for wc in omegas_comp])
	
		inv_alpha = array([calc_invpol(wc, r_np, nind) for wc in omegas_comp])

		#	ev 1/(s*eV) *nm^3* s/nm * 1/(nm**2)
		# NOTE: extra factor of nind might be to compensate definition of alpha
		#ext = (waves/(nind*c*hbar))*imag(inv_alpha**(-1))/(pi*r_np**2)	
		#plt.figure()
		#plt.plot(1240/waves, ext, color="blue")

		waves = 1240/omegas_real
		print(LS_BF.shape)	

		plt.figure()
		plt.plot(waves, real(LS), color="red", label="Re. Polylog")
		plt.plot(waves, imag(LS), color="blue", label="Im. Polylog")
		plt.plot(waves, real(LS_BF[:,0,0]), color="black", ls="--", label="Re. Brute")
		plt.plot(waves, imag(LS_BF[:,0,0]), color="green", ls="--", label="Im. Brute")
		#plt.plot(waves, real(inv_alpha), color="orange")
		#plt.plot(waves, imag(inv_alpha), color="purple")
		plt.legend(loc="best")
		plt.show()
	#END: if test_type == "S_long"
	
	elif test_type == "T_long":
		r_eval = array([t, 0,0])
		print("\nTesting T (long)\n")	
		# Lerch transcendent:
		Tarr = array([T_1d_long(wc, d, t, kx, nind) for wc in omegas_comp])
		# brute force eval:
		T_BF = array([T_1D_x_brute(r_eval, N, d, wc, kx, nind) for wc in omegas_comp])
	
		print(type(Tarr))	
		print(type(Tarr[0]))		
	
		plt.figure()
		waves = 1240/omegas_real
		plt.plot(waves, real(Tarr), color="red", label="Re. Lerch")
		plt.plot(waves, imag(Tarr), color="blue", label="Im. Lerch")
		plt.plot(waves, real(T_BF[:,0,0]), color="black", ls="--", label="Re. Brute")
		plt.plot(waves, imag(T_BF[:,0,0]), color="green", ls="--", label="Im. Brute")
		#plt.plot(waves, real(inv_alpha), color="orange")
		#plt.plot(waves, imag(inv_alpha), color="purple")
		plt.legend(loc="best")
		plt.show()
	#END: if test_type == "T_long"
	
	elif test_type == "T_trans":
		r_eval = array([t, 0,0])
		print("\nTesting T (trans)\n")	
		# Lerch transcendent:
		Tarr = array([T_1d_trans(wc, d, t, kx, nind) for wc in omegas_comp])
		# brute force eval:
		T_BF = array([T_1D_x_brute(r_eval, N, d, wc, kx, nind) for wc in omegas_comp])
	
		plt.figure()
		waves = 1240/omegas_real
		plt.plot(waves, real(Tarr), color="red", label="Re. Lerch")
		plt.plot(waves, imag(Tarr), color="blue", label="Im. Lerch")
		plt.plot(waves, real(T_BF[:,1,1]), color="black", ls="--", label="Re. Brute")
		plt.plot(waves, imag(T_BF[:,1,1]), color="green", ls="--", label="Im. Brute")
		#plt.plot(waves, real(inv_alpha), color="orange")
		#plt.plot(waves, imag(inv_alpha), color="purple")
		plt.legend(loc="best")
		plt.show()
	#END: if test_type == "T_trans"
	
	else:
		print("\nrun_type not found.\n")
		exit()
	#END: else
#fi
