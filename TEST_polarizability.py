from __future__ import division, print_function
from numpy import array, linspace, loadtxt, pi, real, imag
from matplotlib import pyplot as plt
from pyfunc import calc_pol, ic, hbar, c
import pathlib


dir_names = ["20nmD_n1p0", "40nmD_n1p0", "80nmD_n1p0", "20nmD_n1p33", "40nmD_n1p33", "80nmD_n1p33"]
base_name = "/home/mrb179/Projects/Mie/single_NP/Ag/DRUDE_MODEL_CHARLES/"



f1, ax1 = plt.subplots()
f2, ax2 = plt.subplots()

for dir_base in dir_names:

	# generated polarizability files:
	file_path = base_name + dir_base + "/polarizability.out"
	qdat      = loadtxt(file_path)
	scale     = 10**(9)		# converts um^3 to nm^3
	waves_vac = qdat[:,0] 		# [nm]
	k_med 	  = 2*pi/qdat[:,1] 	# [nm^-1]
	pol_real  = qdat[:,2]*scale 	# [nm^3]
	pol_imag  = qdat[:,3]*scale 	# [nm^3]
	np_diam   = float(dir_base.split("nm")[0])
	ax1.plot(waves_vac, pol_real, ls="-", label="{0:s}".format(dir_base))
	ax2.plot(waves_vac, pol_imag, ls="-", label="{0:s}".format(dir_base))

	# construct background refractive index:
	tmp = dir_base.split("_n")[-1]
	whole = float(tmp.split("p")[0])
	decim = float(tmp.split("p")[-1])/100
	nind = whole + decim
	print("")
	print("{0:s}".format(dir_base))
	print("{0:.2f}".format(nind))
	print("") 

	# Analytical polarizability
	omegas_eV = hbar*2*pi*c/waves_vac
	omega_complex = array([complex(ww, 0) for ww in omegas_eV])
	#calc_pol()		
	pol_analyt = array([calc_pol(wc, r=(np_diam/2), nout=nind) for wc in omega_complex ])
	ax1.plot(waves_vac, real(pol_analyt)/nind**2, ls="--")
	ax2.plot(waves_vac, imag(pol_analyt)/nind**2, ls="--")

#end
ax1.set_ylabel("$\Re{\alpha}$")
ax1.set_xlabel("Wavelength (nm)")
ax1.set_xlim((350,650))
ax1.legend(loc="best")
#ax1.savefig("polariz_test_Ag_ext.pdf")

ax2.set_ylabel("$\Im{ \alpha }$")
ax2.set_xlabel("Wavelength (nm)")
ax2.set_xlim((350,650))
ax2.legend(loc="best")
#ax2.savefig("polariz_test_Ag_ext.pdf")

plt.show() 


