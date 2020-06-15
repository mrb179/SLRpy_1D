from __future__ import division, print_function
from numpy import pi, sqrt, real, imag,  sin, cos, array, dot, linspace, matmul, cross, exp, conj, zeros, identity, loadtxt, arcsin, arccos, eye, vstack, savetxt
from numpy.linalg import inv
import pylab as plt

f1_name = "PF_single_Ag_80nmD_n1p0_10nm_trans.out" 
f2_name = "PF_single_Ag_80nmD_n1p0_20nm_trans.out"
f3_name = "PF_single_Ag_80nmD_n1p0_30nm_trans.out"
f4_name = "PF_single_Ag_80nmD_n1p0_30nm_trans_alt.out"
#f4_name = "PF_single_Ag_80nmD_n1p0_50nm_trans.out"
#f5_name = ""
#f6_name = ""

file_list = [f1_name, f2_name, f3_name, f4_name]

plt.figure()
for file_name in file_list:
	dat = loadtxt(file_name, skiprows=1)
	plt.plot(1240/dat[:,0], dat[:,1], label="")
	#plt.xlim((min(waves), max(waves)))
#
plt.xlabel("Energy (eV)")
plt.ylabel("PF")
plt.show()
