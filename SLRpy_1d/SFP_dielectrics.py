#!/usr/bin/python

from scipy import *
from cmath import *

import inspect

#Custom
from SFP import hp,evnm,generate_waves

################################################################################
#   Copyright M. G. Blaber 2011
#   This file is part of SFP = Street Fighting Plasmonics
#   Street Fighting Plasmonics is free software: you can redistribute it 
#   and/or modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 3 of the 
#   License, or (at your option) any later version.
#
#   SFP is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with SFP.  If not, see <http://www.gnu.org/licenses/>.
################################################################################
# 
# SFP_dielectrics.py
#  version 0.2
# 
# This file gives an interpolated fit to dielectric constants of silver and gold
# The original data is that of Johnson & Christy[1] and Palik[2]. The fit is to
# a function by Etchegoin and Leru[3].
# 
#     See the SFP_dielectric_examples.py for examples of how to use this file.
# 
# References
#     [1]  P. B. Johnson and R. W. Christy, 
#             Phys. Rev. B 6, 4370 (1972).
#     [2]  D. W. Lynch, W. R. Hunter in "Handbook of Optical Constants of 
#             Solids"; E. D. Palik, Ed.; Academic Press: 1985; Vol. 1.
#     [3]  P. G. Etchegoin, E. C. Le Ru, and M. Meyer, 
#             J. Chem. Phys. 125, 164705 (2006).
# 
# Definitions: 
#     Ag         = Silver
#     JC         = Johnson and Christy [1]
#     LH         = Palik [2]
#     SS         = Surface scattering
#     GNR        = Gamma_{Non-Radiative} - = 2 eps_2 ( d{eps_1}/d{omega} )
#                  Non-radiative resonance linewidth
#     wavelength = wavelength in nm
#     energy_ev  = energy in ev.
#     gscat      = drude phenomenological scattering rate that is due to
#                  surface scattering
#
# ------------------- PROVIDES THE FOLLOWING FUNCTIONS -------------------------
#
# internal_dielectrics_test()       = Test function to ensure some things are
#                                     working. Not Comprehensive.
# dielectric_deriv(f,x)             = Takes a function f as the argument (see
#                                     below). x is the argument to f at which 
#                                     the derivative is calculated.
# dielectric_deriv_gamma(f,x,gscat) = Same as above but for functions that requi
#                                     re the scattering rate as an argument.
# generate_waves(min,max,dwave)     = returns an array of wavelengths in bet
#                                     ween min and max in steps of dwave. Can of
#                                     course be used for energies.
# generate_dielectric_tab(filename,
#              min,max,dwave,m)     = Builds a ddscat style tab file. from wave
#                                     =min to max in dwave steps. using refract
#                                     ive index m.
# generate_tab(fname,min,max,dwave,
#              function,[m],[gscat])= Builds a ddscat style tab file. from wave
#                                     =min to max in dwave steps. function is 
#                                     one of the funcs below. m is the diel of
#                                     the medium (optional) and gscat is the op
#                                     tional scattering rate.
# Ag_path_to_g(P)                   = Converts the mean free path (P in nm) to 
#                                     a scattering rate (g=gamma) in eV.

#                                      Metal  Author input  EPS  SS   GNR
# Ag_JC(wavelength)                     Ag    JC     nm      1    0    0
# Ag_JC_SS(wavelength,gscat)            Ag    JC     nm,eV   1    1    0
# Ag_JC_ev(energy_ev)                   Ag    JC     eV      1    0    0
# Ag_JC_SS_ev(energy_ev,gscat)          Ag    JC     eV,eV   1    1    0
# Ag_JC_GNR(wavelength)                 Ag    JC     nm      0    0    1
# Ag_JC_SS_GNR(wavelength,gscat)        Ag    JC     nm,eV   0    1    1
# Ag_LH(wavelength)                     Ag    LH     nm      1    0    0
# Ag_LH_SS(wavelength,gscat)            Ag    LH     nm,eV   1    1    0
# Ag_LH_ev(energy_ev)                   Ag    LH     eV      1    0    0
# Ag_LH_SS_ev(energy_ev,gscat)          Ag    LH     eV,eV   1    1    0
# Ag_LH_GNR(wavelength)                 Ag    LH     nm      0    0    1
# Ag_LH_SS_GNR(wavelength,gscat)        Ag    LH     nm,eV   0    1    1
#
#                                      Metal  Author  input  EPS  SS   GNR
# Au_JC(wavelength)                     Au    JC     nm      1    0    0
# Au_JC_SS(wavelength)                  Au    JC     nm      1    1    0
# Au_JC_ev(energy_ev)                   Au    JC     eV      1    0    0
# Au_JC_SS_ev(energy_ev,gscat)          Ag    JC     eV,eV   1    1    0
# Au_JC_GNR(wavelength)                 Ag    JC     nm      0    0    1
# Au_JC_SS_GNR(wavelength,gscat)        Ag    JC     nm,eV   0    1    1
# Au_LH(wavelength)                     Au    LH     nm      1    0    0
# Au_LH_SS(wavelength)                  Au    LH     nm      1    1    0
# Au_LH_ev(energy_ev)                   Au    LH     eV      1    0    0
# Au_LH_SS_ev(energy_ev,gscat)          Ag    LH     eV,eV   1    1    0
# Au_LH_GNR(wavelength)                 Ag    LH     nm      0    0    1
# Au_LH_SS_GNR(wavelength,gscat)        Ag    LH     nm,eV   0    1    1

###########################################################
######################## Gen Tab ##########################
###########################################################

################################################################################
# generate_dielectric_tab(filename,
#              min,max,dwave,m)       = Builds a ddscat style tab file. from wave
#                                     =min to max in dwave steps. using refract
#                                     ive index m.
def generate_dielectric_tab(filename,minwave,maxwave,dwave,m_in):
    eps_in=m_in**2
    waves=generate_waves(minwave,maxwave,dwave)
    nwaves=len(waves)
    
    file = open(filename,'w')
    file.write("# Dielectric with m=%f\n"%(m_in))
    file.write(" 1 0 0 2 3 = columns for wave, Re(n), Im(n), eps1, eps2    \n")
    file.write("    #LAMBDA(um)     eps1(Real)     eps2(Imag)\n")
    for iwave in range(nwaves):
        wave=waves[iwave]
        file.write("%15.8f%15.8f%15.8f\n"%((wave)/1000.0,real(eps_in),imag(eps_in)))
    #for
    file.close()
#def

################################################################################
# generate_tab(fname,min,max,dwave,
#              function,[m],[gscat])= Builds a ddscat style tab file. from wave
#                                     =min to max in dwave steps. function is 
#                                     one of the funcs below. m is the diel of
#                                     the medium (optional) and gscat is the op
#                                     tional scattering rate.
def generate_tab(filename,minwave,maxwave,dwave,func_pointer,m=1.00,gscat=0.00):
    
    nargs=len(inspect.getargspec(func_pointer)[0])
    funcname=func_pointer.__name__
    epsm=sqrt(m)
    
    
    waves=generate_waves(minwave,maxwave,dwave)
    nwaves=len(waves)
    
    file = open(filename,'w')
    file.write("# %s, m=%f, gscat=%f\n"%(funcname,m,gscat))
    file.write(" 1 0 0 2 3 = columns for wave, Re(n), Im(n), eps1, eps2    \n")
    file.write("    #LAMBDA(um)     eps1(Real)     eps2(Imag)\n")
    for iwave in range(nwaves):
        wave=waves[iwave]
        if(nargs==1):
            eps=func_pointer(wave)/epsm
        elif(nargs==2):
            eps=func_pointer(wave,gscat)/epsm
        #endif
        file.write("%15.8f%15.8f%15.8f\n"%((wave/m)/1000.0,real(eps),imag(eps)))
    #for
    file.close()
        
    
#def 


###########################################################
######################### SILVER ##########################
###########################################################

def Ag_path_to_g(P,A):
    # Uses the Fermi velocity of silver to calculate the scattering
    # rate based on the mean free path P.
    # 1e9 
    gamma_scat=(A*(1.39e6*1.e9)/P)*hp
    return gamma_scat
#def

################################################################################
# JOHNSON AND CHRISTY
#
# Ag_JC(wavelength)
# Ag_JC_SS(wavelength,gscat)
#
# Ag_JC_ev(energy_ev)
# Ag_JC_SS_ev(energy_ev,gscat)
# 
# Ag_JC_GNR(wavelegnth)
# Ag_JC_SS_GNR(wavelength,gscat)

###################################
# DIELECTRIC FUNCTIONS (WAVELENGTH)

# Silver - JC - Leru - Wavelength nm
def Ag_JC(wavelength):
    #L is the wavelength and should be in nm.
    l=wavelength
    einf= 1.46932164354844;   einfi= 0.245032510494354
    lp=  137.964321665857;    gp=  74484.2077687798
    a1=  0.163221158200948;   l1=  307.954950886876; g1=  4277.72063494036; p1=  -1.70922134985866
    a2=  0.967277507909697;   l2=  260.730135593913; g2=  1051.4738085748 ; p2=  0.154666427887571
    a3=  0.615093779453933;   l3=  236.212760050413; g3=  1008.62208663784; p3=  -1.78230663012765
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;      ib3=czero;
    cone=complex(0.0,1.0);
    #The Three interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    iba=exp( cone*p3)/ ( 1.0/l3 - 1.0/l  - cone/g3 )
    ibb=exp(-cone*p3)/ ( 1.0/l3 + 1.0/l  + cone/g3 )
    ib3=(a3/l3)*(iba+ibb) #Interband 3
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ib1 + ib2 + ib3
    
    return eps
#def Ag_JC

# Silver - JC - McMahon Parameters - Wavelength nm
def Ag_JC_MM(wavelength):
    # wavelength and should be in nm.
    # speed of light in [nm/s]
    c = 2.99792458*10**(17)
    # convert wavelength to angular frequency:
    w = 2.*pi*c/wavelength
    einf= 1.17152;	wd= 1.39604*10**(16);	gd= (1./(2*pi))*12.6126*10**(13);
    Le1= 2.23994;	Lw1= 8.25718*10**(15); 	Ld1= 1.95614*10**(14);
    Le2= 0.222651;	Lw2= 3.05707*10**(15); 	Ld2= 8.52675*10**(14);
    czero=complex(0.0,0.0)
    cone=complex(0.0,1.0);
    ib1=czero;	ib2=czero;	drude=czero;
    # interband transition 1:
    ib1 = (Le1*Lw1**2)/( Lw1**2 - 2.*cone*w*Ld1 - w**2 )    
    # interband transition 2:
    ib2 = (Le2*Lw2**2)/( Lw2**2 - 2.*cone*w*Ld2 - w**2 ) 
    # Drude component:
    drude = einf -  (wd**2)/( w**2 + cone*w*gd )
    # Dielectric function:
    eps = drude #+ ib1 + ib2
    return eps 
#def Ag_JC_MM

# Silver - JC - Eduardo Coronado Parameters - Wavelength nm
def Ag_JC_EC(wavelength):
    # wavelength and should be in nm.
    # speed of light in [nm/s]
    c = 2.99792458*10**(17)
    # convert wavelength to angular frequency:
    w = 2.*pi*c/wavelength
    einf= 1.00;	wd= 1.38*10**(16);	gd= 2.7*10**(13);
    czero=complex(0.0,0.0)
    cone=complex(0.0,1.0);
    ib1=czero;	ib2=czero;	drude=czero;
    # Drude component:
    drude = einf -  (wd**2)/( w**2 + cone*w*gd )
    # Dielectric function:
    eps = drude + ib1 + ib2
    return eps 
#def Ag_JC_EC

# Silver - JC - Leru - DAMPING - Wavelength nm
def Ag_JC_SS(wavelength,gscat):
    #L is the wavelength and should be in nm.
    #g2 is an additional damping parameter and should be in eV
    l=wavelength
    einf= 1.46932164354844;   einfi= 0.245032510494354
    lp=  137.964321665857;    gp=  evnm/(evnm/74484.2077687798+gscat)
    a1=  0.163221158200948;   l1=  307.954950886876; g1=  4277.72063494036; p1=  -1.70922134985866
    a2=  0.967277507909697;   l2=  260.730135593913; g2=  1051.4738085748 ; p2=  0.154666427887571
    a3=  0.615093779453933;   l3=  236.212760050413; g3=  1008.62208663784; p3=  -1.78230663012765
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;      ib3=czero;
    cone=complex(0.0,1.0);
    #The Three interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    iba=exp( cone*p3)/ ( 1.0/l3 - 1.0/l  - cone/g3 )
    ibb=exp(-cone*p3)/ ( 1.0/l3 + 1.0/l  + cone/g3 )
    ib3=(a3/l3)*(iba+ibb) #Interband 3
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ib1 + ib2 + ib3
    
    return eps
#def Ag_JC

################################################################################
# DIELECTRIC FUNCTIONS (ENERGY)

# Silver - JC - Leru - Energy ev
def Ag_JC_ev(energy_ev):
    return Ag_JC(1239.854/energy_ev)
#def Ag_JC_eV

# Silver - JC - Leru - Energy ev
def Ag_JC_SS_ev(energy_ev,gscat):
    return Ag_JC_SS(1239.854/energy_ev,gscat)
#def Ag_JC_eV

################################################################################
# Nonradiative linewidth (Wavelength)
# 
# Ag_JC_GNR(wavelength)
# Ag_JC_SS_GNR(wavelength,gscat)

def Ag_JC_GNR(wavelength):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Ag_JC(wavelength)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv(Ag_JC_ev, ev)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2.0*eps2/tmp
    
#def

def Ag_JC_SS_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Ag_JC_SS(wavelength,gscat)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Ag_JC_SS_ev, ev,gscat)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2.0*eps2/tmp

#def Ag_JC_eV

################################################################################
# LYNCH AND HUNTER (Palik)
#
# Ag_LH(wavelength)               #No Surf Scattering
# Ag_LH_SS(wavelength,gscat)      #With Surf Scattering
#
# Ag_LH_ev(energy_ev)             #No Surf Scattering
# Ag_LH_SS_ev(energy_ev,gscat)    #With Surf Scattering
# 
# Ag_LH_GNR(wavelength)        #No Surf Scattering
# Ag_LH_SS_GNR(wavelength,gscat)  #With Surf Scattering


# Silver - LH - Leru - Wavelength nm
def Ag_LH(wavelength):
    #L is the wavelength and should be in nm.
    l=wavelength
    einf= 0.983048524950333;    einfi= 0.337143329306706
    lp=  145.777666579854;        gp=  17192.158229537
    a1=  0.153920567777116;        l1=  306.624437448325;    g1=  4647.05951655734;    p1=  -1.60205817994344
    a2=  0.852173299670293;        l2=  265.13703382702;    g2=  793.646031672091;    p2=  -0.479756609595443
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    cone=complex(0.0,1.0);
    #The Three interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ib1 + ib2
    
    return eps
#def Ag_LH

# Silver - LH - Leru - Wavelength nm
def Ag_LH_SS(wavelength,gscat):
    #L is the wavelength and should be in nm.
    #g2 is an additional damping parameter and should be in eV
    l=wavelength
    einf= 0.983048524950333;    einfi= 0.337143329306706
    lp=  145.777666579854;        gp=  evnm/(evnm/17192.158229537+gscat)
    a1=  0.153920567777116;        l1=  306.624437448325;    g1=  4647.05951655734;    p1=  -1.60205817994344
    a2=  0.852173299670293;        l2=  265.13703382702;    g2=  793.646031672091;    p2=  -0.479756609595443
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    cone=complex(0.0,1.0);
    #The Three interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ib1 + ib2
    
    return eps
#def Ag_LH

################################################################################
# DIELECTRIC FUNCTIONS (ENERGY)

# Silver - LH - Leru - Energy ev
def Ag_LH_ev(energy_ev):
    return Ag_LH(1239.854/energy_ev)
#def Ag_LH_eV

# Silver - LH - Leru SURFACE SCAT - Energy ev
def Ag_LH_SS_ev(energy_ev,gscat):
    return Ag_LH_SS(1239.854/energy_ev,gscat)
#def Ag_LH_eV

################################################################################
# Nonradiative linewidth (Wavelength)
# 
# Ag_LH_GNR(wavelength)
# Ag_LH_SS_GNR(wavelength,gscat)

def Ag_LH_GNR(wavelength):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Ag_LH(wavelength)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv(Ag_LH_ev, ev)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2*eps2/tmp
    
#def

def Ag_LH_SS_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Ag_LH_SS(wavelength,gscat)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Ag_LH_SS_ev, ev,gscat)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2*eps2/tmp

#def Ag_LH_eV




###########################################################
########################## GOLD ###########################
###########################################################

def Au_path_to_g(P,A):
    # Uses the Fermi velocity of gold to calculate the scattering
    # rate based on the mean free path P.
    # 1e9 
    gamma_scat=(A*(1.40e6*1.e9)/P)*hp
    return gamma_scat
#def

################################################################################
# JOHNSON AND CHRISTY
#
# Au_JC(wavelength)
# Au_JC_SS(wavelength,gscat)



def Au_JC(wavelength):
    #L is the wavelength and should be in nm.
    l=wavelength
    einf= 1.53;    einfi= 0.0
    lp=145.0;    gp=14500.0
    a1=0.94;    l1=468.0;    g1=2300.0;    p1=-0.785398
    a2=1.36;    l2=331.0;    g2=940.0;    p2=-0.785398
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    cone=complex(0.0,1.0);
    #The Three interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ib1 + ib2
    
    return eps
#def Au_JC

def Au_JC_SS(wavelength,gscat):
    #L is the wavelength and should be in nm.
    l=wavelength
    einf= 1.53;    einfi= 0.0
    lp=145.0;    gp=evnm/(evnm/14500.0+gscat)
    a1=0.94;    l1=468.0;    g1=2300.0;    p1=-0.785398
    a2=1.36;    l2=331.0;    g2=940.0;    p2=-0.785398
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    cone=complex(0.0,1.0);
    #The Three interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ib1 + ib2
    
    return eps
#def Au_JC_SS

# Gold - JC - Leru - Energy ev
def Au_JC_ev(energy_ev):
    return Au_JC(1239.854/energy_ev)
#def Au_JC_eV

# Gold - JC - Leru - Energy ev - Surface Scatterinfd
def Au_JC_SS_ev(energy_ev,gscat):
    return Au_JC_SS(1239.854/energy_ev,gscat)
#def Au_JC_eV


def Au_JC_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Au_JC(wavelength)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Au_JC_ev, ev)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2*eps2/tmp

#def Au_JC_SS_GNR


def Au_JC_SS_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Au_JC_SS(wavelength,gscat)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Au_JC_SS_ev, ev,gscat)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2*eps2/tmp

#def Au_JC_SS_GNR


################################################################################
# Lynch and Hunter
#
# Au_LH(wavelength)
# Au_LH_SS(wavelength,gscat)

def Au_LH(wavelength):
    #L is the wavelength and should be in nm.
    l=wavelength
    einf= 1.786232666464; einfi= 0.474161581212655
    lp=  135.865937688189; gp=  16088.906016254
    a1=  -0.789833578871083; l1=  314.245873751911; g1=  1469.6864723325; p1=  -9.84007589746256
    a2=  1.0499806855811;    l2=  463.507318325189; g2=  2290.02237298882;p2=  -0.951356810003336
    a3=  0.0973282595247586; l3=  57.9550132461995; g3=  497.421718909885;p3=  6.44681517896198
    a4=  -0.789447146394695; l4=  156.777953986622; g4=  283.587000505755;p4=  9.29750037444383
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    ib3=czero;        ib4=czero;
    cone=complex(0.0,1.0);
    #The Four interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    iba=exp( cone*p3)/ ( 1.0/l3 - 1.0/l  - cone/g3 )
    ibb=exp(-cone*p3)/ ( 1.0/l3 + 1.0/l  + cone/g3 )
    ib3=(a3/l3)*(iba+ibb) #Interband 3
    iba=exp( cone*p4)/ ( 1.0/l4 - 1.0/l  - cone/g4 )
    ibb=exp(-cone*p4)/ ( 1.0/l4 + 1.0/l  + cone/g4 )
    ib4=(a4/l4)*(iba+ibb) #Interband 4
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ ib1 + ib2 + ib3 + ib4
    
    return eps
#def Au_LH

def Au_LH_SS(wavelength,gscat):
    #L is the wavelength and should be in nm.
    l=wavelength
    einf= 1.786232666464; einfi= 0.474161581212655
    lp=  135.865937688189; gp=  gp=evnm/(evnm/16088.906016254+gscat)
    a1=  -0.789833578871083; l1=  314.245873751911; g1=  1469.6864723325; p1=  -9.84007589746256
    a2=  1.0499806855811;    l2=  463.507318325189; g2=  2290.02237298882;p2=  -0.951356810003336
    a3=  0.0973282595247586; l3=  57.9550132461995; g3=  497.421718909885;p3=  6.44681517896198
    a4=  -0.789447146394695; l4=  156.777953986622; g4=  283.587000505755;p4=  9.29750037444383
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    ib3=czero;        ib4=czero;
    cone=complex(0.0,1.0);
    #The Four interband transitions
    iba=exp( cone*p1)/ ( 1.0/l1 - 1.0/l  - cone/g1 )
    ibb=exp(-cone*p1)/ ( 1.0/l1 + 1.0/l  + cone/g1 )
    ib1=(a1/l1)*(iba+ibb) #Interband 1
    iba=exp( cone*p2)/ ( 1.0/l2 - 1.0/l  - cone/g2 )
    ibb=exp(-cone*p2)/ ( 1.0/l2 + 1.0/l  + cone/g2 )
    ib2=(a2/l2)*(iba+ibb) #Interband 2
    iba=exp( cone*p3)/ ( 1.0/l3 - 1.0/l  - cone/g3 )
    ibb=exp(-cone*p3)/ ( 1.0/l3 + 1.0/l  + cone/g3 )
    ib3=(a3/l3)*(iba+ibb) #Interband 3
    iba=exp( cone*p4)/ ( 1.0/l4 - 1.0/l  - cone/g4 )
    ibb=exp(-cone*p4)/ ( 1.0/l4 + 1.0/l  + cone/g4 )
    ib4=(a4/l4)*(iba+ibb) #Interband 4
    #Drude Part and Total
    drude=einf +cone*einfi - 1.00 / ( (lp**2.0)*( l**-2.0 +cone/(gp*l) ) )
    eps=drude+ ib1 + ib2 + ib3 + ib4
    return eps
#def Au_LH_SS

# Gold - LH - Energy ev
def Au_LH_ev(energy_ev):
    return Au_LH(1239.854/energy_ev)
#def Au_LH_eV

# Gold - LH -  Energy ev - Surface Scattering
def Au_LH_SS_ev(energy_ev,gscat):
    return Au_LH_SS(1239.854/energy_ev,gscat)
#def Au_LH_eV


def Au_LH_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Au_LH(wavelength)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Au_LH_ev, ev)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2*eps2/tmp

#def Au_LH_SS_GNR


def Au_LH_SS_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Au_LH_SS(wavelength,gscat)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Au_LH_SS_ev, ev,gscat)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2*eps2/tmp

#def Au_LH_SS_GNR


###########################################################
########################## Aluminum #######################
###########################################################

def Al_path_to_g(P,A):
    # Uses the Fermi velocity of aluminum to calculate the scattering
    # rate based on the mean free path P.
    # 1e9 
    gamma_scat=(A*(2.03e6*1.e9)/P)*hp
    return gamma_scat
#def


###########################################################
# Palik
#
def Al_Palik(wavelength):
    #L is the wavelength and should be in nm.
    ev=1239.854/wavelength
    einf= 1.00;    einfi= 0.0
    lp=14.344;    gp=0.651
    a1=4.283;    l1=1.645;    g1=0.203;
    a2=0.069;    l2=3.933;    g2=0.036;
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    cone=complex(0.0,1.0);
    #The Two interband transitions
    iba=a1*l1**2
    ibb=ev*(ev+2*cone*g1)-l1
    ib1=(iba/ibb) #Interband 1
    iba=a2*l2**2
    ibb=ev*(ev+2*cone*g2)-l2
    ib2=(iba/ibb) #Interband 2
    #Drude Part and Total
    drude=einf +cone*einfi - (lp**2)/(ev*(ev+cone*gp))
    eps=drude - ib1 - ib2
    return eps
#def Al_Palik_ev

# Aluminum - Palik - Wavelength
def Al_Palik_ev(energy_ev):
    return Al_Palik(1239.854/energy_ev)
#def Au_JC_eV



def Al_Palik_SS(wavelength,gscat):
    #L is the wavelength and should be in nm.
    ev=1239.854/wavelength
    einf= 1.00;    einfi= 0.0
    lp=14.344;    gp=0.651+gscat
    a1=4.283;    l1=1.645;    g1=0.203;
    a2=0.069;    l2=3.933;    g2=0.036;
    czero=complex(0.0,0.0)
    iba=czero;        ibb=czero;
    ib1=czero;        ib2=czero;
    cone=complex(0.0,1.0);
    #The Two interband transitions
    iba=a1*l1**2
    ibb=ev*(ev+2*cone*g1)-l1
    ib1=(iba/ibb) #Interband 1
    iba=a2*l2**2
    ibb=ev*(ev+2*cone*g2)-l2
    ib2=(iba/ibb) #Interband 2
    #Drude Part and Total
    drude=einf +cone*einfi - (lp**2)/(ev*(ev+cone*gp))
    eps=drude - ib1 - ib2
    return eps
#def Al_Palik_ev

def Al_Palik_SS_ev(energy_ev,gscat):
    return Al_Palik_SS(1239.854/energy_ev,gscat)
#def Au_JC_eV




def Al_Palik_GNR(wavelength):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Al_Palik(wavelength)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv(Al_Palik_ev, ev)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2.0*eps2/tmp
#def

def Al_Palik_SS_GNR(wavelength,gscat):
    # 2 * eps2 / sqrt( (de1/dw)**2 + (de2/dw)**2 )
    eps=Al_Palik_SS(wavelength,gscat)
    eps1=real(eps)
    eps2=imag(eps)
    
    ev=evnm/wavelength
    #derivative with respect to frequency
    deps_dw=dielectric_deriv_gamma(Al_Palik_SS_ev, ev,gscat)
    tmp=real(sqrt( real(deps_dw)**2 + imag(deps_dw)**2 ))
    
    return 2.0*eps2/tmp
#def Al_Palik_eV


#def Cu_JC_SS(wavelength, radius):
#	eps = 
#def


###########################################################
#################### SUPPORT FUNCTIONS ####################
###########################################################


def internal_dielectrics_test():
    print("Ag_JC(600):    ",Ag_JC(600), "Actual: (-15.3712092445+0.303070673894j)")
    print("Ag_LH(600):    ",Ag_LH(600), "   Actual: (-14.037966287+0.9407517842j)")
    print("Au_JC(600):    ",Au_JC(600), " Actual: (-8.92636318208+1.08246173221j)")
    
#def

def dielectric_deriv(f,x):
    #F is the function we're taking the derivative of, x is the parameter to give to f.
    # A ... X ... B
    dyb=0
    dya=5
    dx=0.002
    iloop=0
    while(abs(abs(dyb)-abs(dya))>0.001):
        dx=dx/2.0
        
        xa=x-dx
        xb=x+dx
    
        ya=f(xa)
        y=f(x)
        yb=f(xb)
    
        dya=y-ya
        dyb=yb-y
        
        iloop=iloop+1
        
        if (iloop>10):
            print("Derivative does not converge")
            return nan
        #if
        
    #while
    
    return (dya/dx+dyb/dx)/2.0
#def dielectric_deriv


def dielectric_deriv_gamma(f,x,gscat):
    #F is the function we're taking the derivative of, x and gscat is the parameter to give to f.
    # A ... X ... B
    dyb=0
    dya=5
    dx=0.002
    iloop=0
    while(abs(abs(dyb)-abs(dya))>0.001):
        dx=dx/2.0
        
        xa=x-dx
        xb=x+dx
        
        ya=f(xa,gscat)
        y=f(x,gscat)
        yb=f(xb,gscat)
        
        dya=y-ya
        dyb=yb-y
        
        iloop=iloop+1
        
        if (iloop>10):
            print("Derivative does not converge")
        #if
        
    #while
    return (dya/dx+dyb/dx)/2.0
#def dielectric_deriv


