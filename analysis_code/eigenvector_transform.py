
from scipy import *
import cPickle as pickle
import ConfigParser
import argparse
import sys
import scipy.weave as weave
from scipy.weave.converters import blitz

# SETTINGS --------------------------------------------------------------------

config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
cfgN = config.getint('settings', 'N')
cfgM = config.getint('settings', 'M')
cfgRe = config.getfloat('settings', 'Re')
cfgbeta = config.getfloat('settings','beta')
cfgWeiss = config.getfloat('settings','Weiss')
cfgAmp = config.getfloat('settings', 'Amp')
cfgk = config.getfloat('settings', 'k')

fp.close()

gamma = pi/2. 


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
argparser.add_argument("-kx", type=float, default=cfgk, 
                help='Override kx of the config file')
argparser.add_argument("-dim2", 
                help='specify just yz eigenvector transform', action='store_true')

args = argparser.parse_args()
N = args.N 
M = args.M
Re = args.Re
beta = args.b
Weiss = args.Wi
Amp = args.amp
k = args.kx

baseFileName = '-k{k}-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    k=k, N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

inFilename = 'full-evecs{0}'.format(baseFileName)
outFilename = 'real-evec{0}'.format(baseFileName)

element_number = r_[0:M]
YPOINTS = cos(pi*element_number/(M-1))
zDataPts = 25   # half the number of points we end up with
zDataPts_2D = 50 # 2D output uses double the number of points per period 
yDataPts = 50
xDataPts = 50

# -----------------------------------------------------------------------------


#FUNCTIONS

def to_3Dreal(velA):
    """convert velocity to a 2D real space velocity array"""
    print 'converting velocity to real space...'
    Vel2D = zeros((yDataPts, zDataPts), dtype='complex')
    for zind in xrange(zDataPts):
        for yind in xrange(yDataPts):
            for n in xrange(2*N+1):
                for m in xrange(M):
                    Vel2D[yind, zind] += velA[n*M+m]*cheb(yind, m)*fou(zind, n)
    del m, n, yind, zind

    #double the size in the z dimension to make it easier to see.
    realVel = zeros((xDataPts,yDataPts, 2*zDataPts))
    realVel[0,:,:] = hstack((Vel2D,Vel2D))
    del Vel2D

    print 'filling out x dimension...'
    for xindx in xrange(1,xDataPts):
        realVel[xindx,:,:] = realVel[0,:,:]*exp(1.j*k*x_point(xindx)) +  \
                            conjugate(realVel[0,:,:])*exp(-1.j*k*x_point(xindx))
    del xindx
    realVel[0,:,:] = realVel[0,:,:] + conjugate(realVel[0,:,:])

    return realVel

def c_to_3Dreal(vecIn):

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

    #double the size in the z dimension to make it easier to see.
    realVel = zeros((xDataPts,yDataPts, 2*zDataPts))
    realVel[0,:,:] = hstack((vecOut,vecOut))
    del vecOut

    print 'filling out x dimension...'
    for xindx in xrange(1,xDataPts):
        realVel[xindx,:,:] = realVel[0,:,:]*exp(1.j*k*x_point(xindx)) +  \
                            conjugate(realVel[0,:,:])*exp(-1.j*k*x_point(xindx))
    del xindx
    realVel[0,:,:] = realVel[0,:,:] + conjugate(realVel[0,:,:])

    return realVel

def c_to_2Dreal(vecIn):
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
    zLength = 4 #2Pi/gamma where gamma = pi/2
    return zLength*((zIndex/(1.*(zDataPts-1.))) - 0.5)

def x_point(xIndex):
    if k==0:
        xLength = 10
    else:
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
    delim = ', '
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

def save_csv_for_paraview((xVec, yVec, zVec), filename):

    Nst = 50
    xspace = 4.*pi/(Nst)
    yspace = 2./(Nst-1)
    zspace = 8./(Nst-1)

    xPoints = r_[0 : 4.*pi : xspace] 
    yPoints = r_[-1 : 1 : yspace]
    zPoints = r_[-4 : 4 : zspace]

    fp = open(filename, 'w')
    fp.write('x,    y,   z,   u,   v,   w\n')
    for xIndx in range(xDataPts):
        for yIndx in range(yDataPts):
            for zIndx in range(2*zDataPts):
                fp.write('{0:10.5f}, '.format(xPoints[xIndx]))
                fp.write('{0:10.5f}, '.format(yPoints[yIndx]))
                fp.write('{0:10.5f}, '.format(zPoints[zIndx]))
                fp.write('{0:10.5f}, '.format(xVec[xIndx, yIndx, zIndx]))
                fp.write('{0:10.5f}, '.format(yVec[xIndx, yIndx, zIndx]))
                fp.write('{0:10.5f}'.format(zVec[xIndx, yIndx, zIndx]))
                fp.write('\n')
        fp.flush()
    del xIndx, yIndx, zIndx

    fp.close()



#MAIN


(du,dv,dw,_,dcxx,dcyy,dczz,dcxy,dcxz,dcyz) = pickle.load(open(inFilename, 'r'))

print 'finished reading in files'


if not args.dim2:
    # Vorticity
    oneOverC = ones(M)
    oneOverC[0] = 1. / 2.
    #set up the CFunc function: 2 for m=0, 1 elsewhere:
    CFunc = ones(M)
    CFunc[0] = 2.
    zLength = 2.*pi/gamma
    vecLen = (2*N+1)*M
    MDY = mk_diff_y()
    MDZ = mk_diff_z()

    print 'vorticity'
    dOmegaX = dot(MDY,dw) - dot(MDZ,dv)
    dOmegaY = dot(MDZ,du) - k*dw
    dOmegaZ = k*dv - dot(MDY,du)

    dOmegaX = c_to_3Dreal(dOmegaX)
    dOmegaY = c_to_3Dreal(dOmegaY)
    dOmegaZ = c_to_3Dreal(dOmegaZ)

    print 'du'
    du = c_to_3Dreal(du)
    print 'dv'
    dv = c_to_3Dreal(dv)
    print 'dw'
    dw = c_to_3Dreal(dw)
    print 'dcxx'
    dcxx = c_to_3Dreal(dcxx)
    print 'dcyy'
    dcyy = c_to_3Dreal(dcyy)
    print 'dczz'
    dczz = c_to_3Dreal(dczz)
    print 'dcxy'
    dcxy = c_to_3Dreal(dcxy)
    print 'dcxz'
    dcxz = c_to_3Dreal(dcxz)
    print 'dcyz'
    dcyz = c_to_3Dreal(dcyz)

    save_csv_for_paraview((du,dv,dw), outFilename[:-7]+'.txt')

    outVortFilename = 'vorticity-evec{0}.txt'.format(baseFileName[:-7])
    save_csv_for_paraview((dOmegaX, dOmegaY, dOmegaZ), outVortFilename)

    pickle.dump((du,dv,dw,dcxx,dcyy,dczz,dcxy,dcxz,dcyz), open(outFilename, 'w'))

    pickle.dump((dOmegaX, dOmegaY, dOmegaZ),
                open('vorticity-evec{0}'.format(baseFileName), 'w'))

if args.dim2:
    print "performing only the 2D transform"
    print "setting zDataPts to zDataPts_2D"
    zDataPts = zDataPts_2D

    # Vorticity
    oneOverC = ones(M)
    oneOverC[0] = 1. / 2.
    #set up the CFunc function: 2 for m=0, 1 elsewhere:
    CFunc = ones(M)
    CFunc[0] = 2.
    zLength = 2.*pi/gamma
    vecLen = (2*N+1)*M
    MDY = mk_diff_y()
    MDZ = mk_diff_z()

    print 'vorticity'
    dOmegaX = dot(MDY,dw) - dot(MDZ,dv)
    dOmegaY = dot(MDZ,du) - k*dw
    dOmegaZ = k*dv - dot(MDY,du)

    dOmegaX = c_to_2Dreal(dOmegaX)
    dOmegaY = c_to_2Dreal(dOmegaY)
    dOmegaZ = c_to_2Dreal(dOmegaZ)

    print 'du'
    du = c_to_2Dreal(du)
    print 'dv'
    dv = c_to_2Dreal(dv)
    print 'dw'
    dw = c_to_2Dreal(dw)
    print 'dcxx'
    dcxx = c_to_2Dreal(dcxx)
    print 'dcyy'
    dcyy = c_to_2Dreal(dcyy)
    print 'dczz'
    dczz = c_to_2Dreal(dczz)
    print 'dcxy'
    dcxy = c_to_2Dreal(dcxy)
    print 'dcxz'
    dcxz = c_to_2Dreal(dcxz)
    print 'dcyz'
    dcyz = c_to_2Dreal(dcyz)

    outFilename = 'real-evec2D{0}'.format(baseFileName)
    pickle.dump((du,dv,dw,dcxx,dcyy,dczz,dcxy,dcxz,dcyz), open(outFilename, 'w'))

    pickle.dump((dOmegaX, dOmegaY, dOmegaZ),
                open('vorticity-evec2D{0}'.format(baseFileName), 'w'))
