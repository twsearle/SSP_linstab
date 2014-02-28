
from scipy import *
import cPickle as pickle
import ConfigParser
import sys
import RStransform

#FUNCTIONS

def to_real(velA):
    """convert velocity to a 2D real space velocity array"""
    print 'converting velocity to real space...'
    realVel = zeros((yDataPts, zDataPts), dtype='complex')
    for zind in xrange(zDataPts):
        for yind in xrange(yDataPts):
            for n in xrange(2*N+1):
                for m in xrange(M):
                    realVel[yind, zind] += velA[n*M+m]*cheb(yind, m)*fou(zind, n)
    del m, n, yind, zind

    return real(realVel)

def faster_FC_transform(vecin) :
    """
    This uses f2py. much faster but uses 64 bit floats instead of doubles so it
    is a little less accurate.
    """

    vecout = RStransform.rstransform(vecin, N, M, numZs, numYs, gamma)

    return real(vecout)


def cheb(yIndex, chebIndex):
    """Take a yIndex(cies) in the array, change it into a y value in the system,
    then calculate the Chebyshev polynomial."""
    return cos(chebIndex*arccos(-1. + (2./(yDataPts-1.))*yIndex))

def y_point(yIndex):
    return -1. + (2./(yDataPts-1.))*yIndex

def z_point(zIndex):
    zLength = 2.*pi/gamma
    return zLength*((zIndex/(1.*(zDataPts-1.))) - 0.5)

def x_point(xIndex):
    xLength = 2.*pi/k 
    return xIndex/(xDataPts-1.)*xLength

def fou(zIndex, fouIndex):
    return exp(1j*(2*pi)*(fouIndex-N)*(((1.*zIndex)/(zDataPts-1.))-0.5))

def save_field(mat, name):
    """takes a matrix of real values in the y and z planes and saves them 
    in a file such that they are readable by gnuplot."""

    settingsline = '#Reynolds: '+str(Re)+' beta: '+str(beta)+' Weissenburg: '\
                +str(Weiss)+ ' Amp: '+str(Amp)

    #Open file, write data, close file
    f = open(name+'.dat', 'w')
    delim = ', \t\t'
    f.write(settingsline)
    f.write('\n')
    for n in range(len(mat[0,:])):
        for m in range(len(mat[:,0])):
            f.write(str(z_point(n)))
            f.write(delim)
            f.write(str(y_point(m)))
            f.write(delim)
            f.write(str(mat[m,n]))
            f.write('\n')
        f.write('\n')
    f.close()

#MAIN
config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
N = config.getint('settings', 'N')
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')
k = config.getfloat('settings', 'k')
gamma = pi/2.


fp.close()

filename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

outFileName = 'real-pf{0}'.format(filename)


element_number = r_[0:M]
YPOINTS = cos(pi*element_number/(M-1))
zDataPts = 25
yDataPts = 50
numZs = zDataPts
numYs = yDataPts


(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(open('pf'+filename, 'r'))

Nselect =1 
# Set all components to zero apart from the selected Fourier mode.
#U[:(N-Nselect)*M]                     = 0
#U[(N-Nselect+1)*M:(N+Nselect)*M]      = 0
#U[(N+Nselect+1)*M:]                   = 0 
#V[:(N-Nselect)*M]                     = 0
#V[(N-Nselect+1)*M:(N+Nselect)*M]      = 0
#V[(N+Nselect+1)*M:]                   = 0 
#W[:(N-Nselect)*M]                     = 0
#W[(N-Nselect+1)*M:(N+Nselect)*M]      = 0
#W[(N+Nselect+1)*M:]                   = 0 
#Cxx[:(N-Nselect)*M]                 = 0
#Cxx[(N-Nselect+1)*M:(N+Nselect)*M]  = 0  
#Cxx[(N+Nselect+1)*M:]               = 0  
#Cyy[:(N-Nselect)*M]                 = 0
#Cyy[(N-Nselect+1)*M:(N+Nselect)*M]  = 0  
#Cyy[(N+Nselect+1)*M:]               = 0  
#Czz[:(N-Nselect)*M]                 = 0
#Czz[(N-Nselect+1)*M:(N+Nselect)*M]  = 0
#Czz[(N+Nselect+1)*M:]               = 0  
#Cxy[:(N-Nselect)*M]                 = 0
#Cxy[(N-Nselect+1)*M:(N+Nselect)*M]  = 0
#Cxy[(N+Nselect+1)*M:]               = 0  
#Cxz[:(N-Nselect)*M]                 = 0
#Cxz[(N-Nselect+1)*M:(N+Nselect)*M]  = 0
#Cxz[(N+Nselect+1)*M:]               = 0  
#Cyz[:(N-Nselect)*M]                 = 0
#Cyz[(N-Nselect+1)*M:(N+Nselect)*M]  = 0
#Cyz[(N+Nselect+1)*M:]               = 0  

#V     = zeros((2*N+1)*M, dtype='complex')
#W     = zeros((2*N+1)*M, dtype='complex')
Cyy   = zeros((2*N+1)*M, dtype='complex')
Czz   = zeros((2*N+1)*M, dtype='complex')

Cxy   = zeros((2*N+1)*M, dtype='complex')
Cxz   = zeros((2*N+1)*M, dtype='complex')
Cyz   = zeros((2*N+1)*M, dtype='complex')

#(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(open('real_base_profile.pickle', 'r'))

print 'finished reading in files'

#U2 = to_real(U)
U = faster_FC_transform(U)
#print shape(U), shape(U2)
#print allclose(U,U2,atol=1e-6)
#print U
#print U2
#exit(1)

#save_field(U, 'U')
V = faster_FC_transform(V)
W = faster_FC_transform(W)
Cxx = faster_FC_transform(Cxx)
#save_field(Cxx, 'Cxx')
Cyy = faster_FC_transform(Cyy)
#save_field(Cyy, 'Cyy')
Czz = faster_FC_transform(Czz)
#save_field(Czz, 'Czz')
Cxy = faster_FC_transform(Cxy)
#save_field(Cxy, 'Cxy')
Cxz = faster_FC_transform(Cxz)
#save_field(Cxz, 'Cxz')
Cyz = faster_FC_transform(Cyz)
#save_field(Cyz, 'Cyz')

pickle.dump((U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz), open(outFileName, 'w'))

#f = open('VW'+'.dat', 'w')
#delim = ', \t\t'
#for n in r_[0:zDataPts:5]:
#    for m in r_[0:yDataPts:5]:
#        f.write(str(z_point(n)))
#        f.write(delim)
#        f.write(str(y_point(m)))
#        f.write(delim)
#        f.write(str(W[m,n]))
#        f.write(delim)
#        f.write(str(V[m,n]))
#        f.write('\n')
#    f.write('\n')
#f.close()

#xz plane (top down) with array[xindx, zindx]
xDataPts = 100

topdownU = zeros((xDataPts, zDataPts), dtype='complex')

#for xindx in xrange(xDataPts):
#    topdownU[xindx, :] = U[yDataPts/2-1, :]
#del xindx
#
#fs = open('U_xz_plane.dat', 'w')
#delim = ', '
#for zindx in xrange(zDataPts):
#    for xindx in xrange(xDataPts):
#        fs.write(str(z_point(zindx)))
#        fs.write(delim)
#        fs.write(str(x_point(xindx)))
#        fs.write(delim)
#        fs.write(str(real(topdownU[xindx, zindx])))
#        fs.write('\n')
#    fs.write('\n')
#fs.flush()
#fs.close()

