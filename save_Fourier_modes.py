#----------------------------------------------------------------------------#
#   Fully spectral space - Fourier Modes                          16-02-13
#   transformer                                             fourier_modes.py
#   Last modified: Mon 04 Mar 2013 12:32:39 GMT
#----------------------------------------------------------------------------#
"""Convert spectral space to Fourier modes and then save these to files for 
plotting"""


# MODULES
import sys
from scipy import *
import cPickle as pickle
import ConfigParser

# FUNCTIONS

def cheb(yIndex, chebIndex):
    """Take a yIndex(cies) in the array, change it into a y value in the system,
    then calculate the Chebyshev polynomial."""
    return cos(chebIndex*arccos(-1. + (2./(yDataPts-1.))*yIndex))

def y_point(yIndex):
    return -1. + (2./(yDataPts-1.))*yIndex
#
# MAIN
#
config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
N = config.getint('settings', 'N')
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')

fp.close()

yDataPts = 80 # Number of y points after transformation to real space

filename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

(U,V,W,conxx,conyy,conzz,conxy,conxz,conyz) = pickle.load(open('pf'+filename, 'r'))

transformed = zeros((yDataPts,2*N+1), dtype=complex)

for n in range(2*N+1):
    for y in range(yDataPts):
        for m in range(M):
            transformed[y, n] += conxx[n*M+m]*cheb(y,m)
del n, m, y

yPtsArray = y_point(r_[0:yDataPts])
for n in range(2*N+1):
    savingarray = vstack((yPtsArray, real(transformed[:,n]),\
                          imag(transformed[:,n]))).T
    savetxt('fou'+str(n-N)+filename[:-7]+'.dat', savingarray)
del n 
