#----------------------------------------------------------------------------#
#   Fully spectral space - Chebychev Modes                          16-02-13
#   transformer                                             fourier_modes.py
#   Last modified: Mon 04 Mar 2013 12:29:03 GMT
#----------------------------------------------------------------------------#
"""Convert spectral space to Chebyschev modes and then save these to files for 
plotting"""


# MODULES
import sys
from scipy import *
import cPickle as pickle
import ConfigParser

# FUNCTIONS

def z_point(zIndex):
    return zLength*((zIndex/(1.*(zDataPts-1.))) - 0.5)

def fou(zIndex, fouIndex):
    return exp(1j*(2*pi)*(fouIndex-N)*(((1.*zIndex)/(zDataPts-1.))-0.5))

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

zDataPts = 100 # Number of y points after transformation to real space
gamma = pi / 2.
zLength = 2.*pi/gamma

filename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

(U,V,W,conxx,conyy,conzz,conxy,conxz,conyz) = pickle.load(open('pf'+filename, 'r'))

transformed = zeros((zDataPts,M), dtype=complex)

for m in range(M):
    for z in range(zDataPts):
        for n in range(2*N+1):
            transformed[z, m] += conxx[n*M+m]*fou(z, n)
del n, m, z

zPtsArray = z_point(r_[0:zDataPts])
for m in range(M):
    filearray = vstack((zPtsArray, real(transformed[:,m]),\
                          imag(transformed[:,m]))).T
    savetxt('cheb'+str(m)+filename[:-7]+'.dat', filearray)
del m 
