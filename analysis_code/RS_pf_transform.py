
from scipy import *
import cPickle as pickle
import ConfigParser
import sys
import RStransform
import scipy.weave as weave
from scipy.weave.converters import blitz
import argparse

# SETTINGS -------------------------------------------------------------------

config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
cfgN = config.getint('settings', 'N')
cfgM = config.getint('settings', 'M')
cfgRe = config.getfloat('settings', 'Re')
cfgbeta = config.getfloat('settings','beta')
cfgWeiss = config.getfloat('settings','Weiss')
cfgAmp = config.getfloat('settings', 'Amp')
gamma = pi/2.

fp.close()

argparser = argparse.ArgumentParser()

argparser.add_argument("-N", type=int, default=cfgN, 
                help='Override Number of Fourier modes given in the config file')
argparser.add_argument("-M", type=int, default=cfgM, 
                help='Override Number of Chebyshev modes in the config file')
argparser.add_argument("-Re", type=float, default=cfgRe, 
                help="Override Reynold's number in the config file") 
argparser.add_argument("-b", type=float, default=cfgbeta, 
                help='Override beta of the config file')
argparser.add_argument("-Wi", type=float, default=cfgWeiss, 
                help='Override Weissenberg number of the config file')
argparser.add_argument("-amp", type=float, default=cfgAmp,
                help='Override amplitude of the streamwise vortices from the config file')

args = argparser.parse_args()
N = args.N 
M = args.M
Re = args.Re
beta = args.b
Weiss = args.Wi
Amp = args.amp


filename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

outFileName = 'real-pf{0}'.format(filename)

# ----------------------------------------------------------------------------

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

def c_FC_transform(vecIn):
    """
    use scipy weave to do the fourier transform. Blitz arrays ought to be
    accurate and fast
    """

    weaveCode = r"""
    #include <cmath>
    #include <complex>
    int n, m, zind, yind;
    double cheb, tmp, tmp2;
    double PI;
    PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
    for (zind=0; zind<zDataPts; zind++) {
        for (yind=0; yind<yDataPts; yind++) {
            for (n=0; n<2*N+1; n++) {
                for (m=0; m<M; m++) {
                tmp = acos(-1. + (2./((double)yDataPts-1.))*(double)yind);
                cheb = cos((double)m*tmp);
                tmp2=2.*PI*((double)n-(double)N)*(((double)zind/((double)zDataPts-1.))-0.5);
                fou.imag() = sin(tmp2);
                fou.real() = cos(tmp2);
                vecOut(yind, zind) += vecIn(n*M+m)*cheb*fou;
                }
            }
        }
    }
    """
    vecOut = zeros((yDataPts, zDataPts), dtype='D', order='C')
    EI = 1.j
    fou = 1. + 1.j
    weave.inline(weaveCode,['vecOut', 'vecIn', 'N', 'M', 'zDataPts',
                            'yDataPts', 'fou'],type_converters = blitz, compiler='gcc',
                 headers=["<cmath>", "<complex>" ] )

    return real(vecOut)

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

def mk_single_diffy():
    """Makes a matrix to differentiate a single vector of Chebyshev's, 
    for use in constructing large differentiation matrix for whole system"""
    # make matrix:
    mat = zeros((M, M), dtype='d')
    for m in range(M):
        for p in range(m+1, M, 2):
            mat[m,p] = 2*p*oneOverC[m]

    return mat

def mk_diff_y():
    """Make the matrix to differentiate a velocity vector wrt y."""
    D = mk_single_diffy()
    MDY = zeros( (vecLen,  vecLen) )
     
    for cheb in range(0,vecLen,M):
        MDY[cheb:cheb+M, cheb:cheb+M] = D
    del cheb
    return MDY

def mk_diff_z():
    """Make matrix to do fourier differentiation wrt z."""
    MDZ = zeros( (vecLen, vecLen), dtype='complex')

    n = -N
    for i in range(0, vecLen, M):
        MDZ[i:i+M, i:i+M] = eye(M, M, dtype='complex')*n*gamma*1.j
        n += 1
    del n, i
    return MDZ


#MAIN

# Set the oneOverC function: 1/2 for m=0, 1 elsewhere:
oneOverC = ones(M)
oneOverC[0] = 1. / 2.
#set up the CFunc function: 2 for m=0, 1 elsewhere:
CFunc = ones(M)
CFunc[0] = 2.
zLength = 2.*pi/gamma
vecLen = (2*N+1)*M

element_number = r_[0:M]
YPOINTS = cos(pi*element_number/(M-1))
zDataPts = 25
yDataPts = 50
numZs = zDataPts
numYs = yDataPts


(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(open('pf'+filename, 'r'))

#Nselect =1 
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
#Cyy   = zeros((2*N+1)*M, dtype='complex')
#Czz   = zeros((2*N+1)*M, dtype='complex')

#Cxy   = zeros((2*N+1)*M, dtype='complex')
#Cxz   = zeros((2*N+1)*M, dtype='complex')
#Cyz   = zeros((2*N+1)*M, dtype='complex')

#(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(open('real_base_profile.pickle', 'r'))

print 'finished reading in files'

# Vorticity
MDY = mk_diff_y()
MDZ = mk_diff_z()

omegaX = dot(MDY, W) - dot(MDZ, V)
# MDX = 0 because base profile no x dependence 
omegaY = dot(MDZ, U) 
omegaZ = - dot(MDY, U)

omegaX = c_FC_transform(omegaX)
omegaY = c_FC_transform(omegaY)
omegaZ = c_FC_transform(omegaZ)

#U2 = to_real(U)
U = c_FC_transform(U)
#print shape(U1), shape(U2)
#print allclose(U1,U2)
#print U1
#print U2
#U = faster_FC_transform(U)

#save_field(U, 'U')
V = c_FC_transform(V)
W = c_FC_transform(W)
Cxx = c_FC_transform(Cxx)
#save_field(Cxx, 'Cxx')
Cyy = c_FC_transform(Cyy)
#save_field(Cyy, 'Cyy')
Czz = c_FC_transform(Czz)
#save_field(Czz, 'Czz')
Cxy = c_FC_transform(Cxy)
#save_field(Cxy, 'Cxy')
Cxz = c_FC_transform(Cxz)
#save_field(Cxz, 'Cxz')
Cyz = c_FC_transform(Cyz)
#save_field(Cyz, 'Cyz')


pickle.dump((U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz), open(outFileName, 'w'))

pickle.dump((omegaX, omegaY, omegaZ), open('vorticity-pf'+filename, 'w'))

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

