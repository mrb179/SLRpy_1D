#!/usr/bin/python
################################################################################
#   Copyright M. G. Blaber 2011
#   This file is part of SFP = Street Fighting Plasmonics
#   Street Fighting Plasmonics is free software: you can redistribute it 
#   and/or modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 3 of the 
#   License, or (at your option) any later version.
#
#   Foobar is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

from numpy import zeros
from numpy import concatenate
#from scipy.interpolate import interp1d
debug=0

################################################################################
# 
# SFP.py
#  version 0.1
# 
# This file provides constants and some basic functions.
#
#     See the SFP_examples.py for examples of how to use this file.
# 
# References
#     [1]  
# 
# Definitions: 
#     wavelength = wavelength in nm
#     energy_ev  = energy in ev.
#     gscat      = drude phenomenological scattering rate that is due to
#                  surface scattering
#
# ------------------- PROVIDES THE FOLLOWING FUNCTIONS -------------------------
#
# generate_waves(minwave,maxwave,dwave) 
# 
# readfile(file, colx, coly)        = returns an array of the file "file". Col x
#                                     is the column number in the file that corre
#                                     corresponds to x. etc.
#
# readAllfile(file)                  = returns an array of the file "file". 
#
# interpolate(file, colx, coly, 
#             minx, maxx, dx )      = returns an interpolated version of the file
#                                     Col x is the column number in the file that
#                                     corresponds to x. minx, maxx are the 
#                                     boundaries of the interpolated data. dx is 
#                                     the step size.
#                                    
# smooth(x,window_len=10,
#        window='hanning')          = x is 1D array of data to smooth.
#                                     window_len is the number of points to 
#                                     average the data over.
# 
################################################################################
#Constants

#Convert ev to nm and vice versa
evnm=1239.854

# Plancks constant
hp=4.135667516E-15;


###########################################################
#################### Generate _waves ######################
###########################################################

def generate_waves(minwave,maxwave,dwave):
	# Wave stetup
	nwaves=int((maxwave-minwave)/dwave)+1
	#print "nwave = ", nwaves
	waves=zeros(nwaves)
	l=minwave
	i=0
	while l <= maxwave:
		waves[i]=l
		l=l+dwave
		i=i+1
	#while
	return waves
#done

###########################################################
###################### readfile ###########################
###########################################################

def readAllFile(filename):

    in_data=None

#This is pretty much the fastest way possible to read a file.
    file = open(filename)
    iline=0
    idataline=0
    for line in file:
        iline=iline+1
        s=line.strip().split()
        if(len(s)!=0):
            char1=s[0][0]
            if ( char1 == "#" or char1 == "!" or char1=="%" ):
                #print("%s"%(line))
                donothing=1 #This is a placeholder. It  does nothing..
            else:
                if (in_data == None):
                    in_data = zeros((0,len(s)))
                #endif
                #Make sure there are enough columns on this line
                if( len(in_data) >0 and len(in_data[0,:]) != len(s) ):
                    print("Error on line %d. Wrong number of columns"%(iline))
                else:
                    idataline=idataline+1
                    #Convert text to numbers.
                    for x in range(len(s)):
                        s[x]=float(s[x])
                    #for
                    
                    in_data=concatenate( ( in_data,[s] ) )
                #if
                
            #if #char
        #if len 0
    #for

    ndata=len(in_data)
    if(debug):
        print("in_data",in_data)
        print("Length of the file in lines = %d"%(iline))
        print("%7d:%14.8f %14.8f"%(0,in_data[0,0],in_data[0,1]))
        print("%7d:%14.8f %14.8f"%(ndata,in_data[ndata-1][0],in_data[ndata-1][1]))
    #if debug
    
    return in_data
    
#end function


###########################################################
###################### readfile ###########################
###########################################################

def readfile(filename, colx, coly):

    in_data=zeros((0,2))

#This is pretty much the fastest way possible to read a file.
    file = open(filename)
    iline=0
    idataline=0
    for line in file:
        iline=iline+1
        s=line.strip().split()
        if(len(s)!=0):
            char1=s[0][0]
            if ( char1 == "#" or char1 == "!" or char1=="%" ):
        		#print("%s"%(line))
                donothing=1 #This is a placeholder. It  does nothing..
            else:
                #Make sure there are enough columns on this line
                if (len(s) < max(colx,coly) ):
                    print("Error on line %d. Not enough columns"%(iline))
                else:
                    idataline=idataline+1
                    in_data=concatenate( ( in_data, [[float(s[colx-1]),float(s[coly-1])]] ) )
                #if
            #if
        #if
    #for

    ndata=len(in_data)
    if(debug):
        print("Length of the file in lines = %d"%(iline))
        print("%7d:%14.8f %14.8f"%(0,in_data[0][0],in_data[0][1]))
        print("%7d:%14.8f %14.8f"%(ndata,in_data[ndata-1][0],in_data[ndata-1][1]))
    #if
    
    return in_data
    
#end function

###########################################################
###################### readfile ###########################
###########################################################

def savefile(filename,data):

#Check that the dimensionality of the data is 2.
    assert  len(data.shape) == 2 , "Data must have only two dimensions"
    
    #Which dimension is longest?
    if(data.shape[0]>=data.shape[1]):
        longdim=0
        shortdim=1
    else:
        longdim=1
        shortdim=0
    #if

    file = open(filename,'w')
    
    for idata in range(data.shape[longdim]):
    
        for icol in range(data.shape[shortdim]):
            
            if (longdim==0):
                file.write(" %f "%(data[idata,icol]))
            else:
                file.write(" %f "%(data[icol,idata]))
            #fi
        #for
        file.write("\n")
    #for

#end function

###########################################################
######################## Interpolate ######################
###########################################################

def interpolate(filename, colx, coly, minx,maxx,dx):

    in_data=readfile(filename,colx,coly)
    
    ndata=len(in_data)
    
    
    #print "Begin Interpolation."
    
    #Check data boundaries
    if( minx < in_data[0][0] or maxx > in_data[ndata-1][0] ):
        print("Min X or Max X is outside the domain of the data in the file")
        sys.exit()
    #if
    
    newx=generate_waves(minx,maxx,dx)
    
    idmin=0
    #Reduce the necessary range in the data.
    for idat in range (len(in_data)):
        if(in_data[idat,0]<minx):
            idmin=idat
        else:
            break
        #if
    #for
    
    idmax=len(in_data)
    #Reduce the necessary range in the data. (Count down)
    for idat in range (len(in_data)-1,-1,-1):
        if(in_data[idat,0]>maxx):
            idmax=idat
        else:
            break
        #if
    #for
    
    #New Range
    #print "New Range = (%d, %d) = (%f, %f)"%(idmin,idmax,in_data[idmin,0],in_data[idmax,0])
    
    interped = interp1d(in_data[idmin:idmax+1,0], in_data[idmin:idmax+1,1], kind='cubic')
    
    newy=interped(newx)

    out_data=zeros((len(newx),2))    
    out_data[:,0]=newx
    out_data[:,1]=newy
    
    return out_data
    
#end function





###########################################################
##################### Smooth ##############################
###########################################################

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    #Remove the additional data due to the windowing
    if window_len%2 == 0:
      cutoff_bottom=(window_len)/2
      cutoff_top=(window_len)/2
    else:
      cutoff_bottom=(window_len-1)/2
      cutoff_top=(window_len-1)/2
    #if
    y=y[cutoff_bottom:len(y)-cutoff_top]
    return y

#end function





