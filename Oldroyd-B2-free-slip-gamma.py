#----------------------------------------------------------------------------#
#   Fully Spectral Linear stability analysis Oldroyd B Model
#   Last modified: Tue  4 Mar 10:59:12 2014
#----------------------------------------------------------------------------#
""" Perform Linear stability analysis to find eigenvalues for the stability 
of the streaky flow"""

# MODULES
import sys
import time
from scipy import *
from scipy import linalg
from scipy.sparse import linalg as sparselinalg
import cPickle as pickle
import ConfigParser
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
cfgk = config.getfloat('settings', 'k')
cfgPiDivide = config.getfloat('settings', 'pi divide')

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
argparser.add_argument("-kx", type=float, default=cfgk,
                help='Override kx from the config file')
argparser.add_argument("-piDivide", type=float, default=cfgPiDivide,
                help='Override piDivide from the config file')
argparser.add_argument("-evecs", 
                help = 'output eigenvectors instead of eigenvalues',
                       action="store_true")

args = argparser.parse_args()
N = args.N 
M = args.M
Re = args.Re
beta = args.b
Weiss = args.Wi
Amp = args.amp
k = args.kx
piDivide = args.piDivide


filename = '-N{N}-M{M}-Re{Re}-b{beta}-Wi{Weiss}-amp{Amp}-gdiv{gdiv}.pickle'.format(\
            N=N,M=M,Re=Re,beta=beta,Weiss=Weiss,Amp=Amp, gdiv=piDivide)

#FUNCTIONS

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

def cheb_prod_mat(velA):
    """Function to return a matrix for left-multiplying two matrices
    of velocities."""

    D = zeros((M, M), dtype='complex')

    #failcount = 0
    for n in range(M):
        for m in range(-M+1,M):     # Bottom of range is inclusive
            itr = abs(n-m)
            if (itr < M):
                D[n, abs(m)] += 0.5*oneOverC[n]*CFunc[itr]*CFunc[abs(m)]*velA[itr]
    del m, n, itr
    return D

def prod_mat(velA):
    """Function to return a matrix ready for the left dot product with another
    velocity vector"""
    MM = zeros((vecLen, vecLen), dtype='complex')

    #First make the middle row
    midMat = zeros((M, vecLen), dtype='complex')
    for n in range(2*N+1):       # Fourier Matrix is 2*N+1 cheb matricies
        yprodmat = cheb_prod_mat(velA[n*M:(n+1)*M])
        endind = 2*N+1-n
        midMat[:, (endind-1)*M:endind*M] = yprodmat
    del n

    #copy matrix into MM, according to the matrix for spectral space
    # top part first
    for i in range(0, N):
        MM[i*M:(i+1)*M, :] = column_stack((midMat[:, (N-i)*M:], zeros((M, (N-i)*M))) )
    del i
    # middle
    MM[N*M:(N+1)*M, :] = midMat
    # bottom 
    for i in range(0, N):
        MM[(i+N+1)*M:(i+2+N)*M, :] = column_stack((zeros((M, (i+1)*M)), midMat[:, :(2*N-i)*M] ))
    del i

    return MM

def mk_bigM():
    """ makes the matrix to appear on the left hand side of the generalised
    eigenvalue problem"""
    bigM = zeros((10*vLen, 10*vLen), dtype=complex)

    ################Navier Stokes x direction:###############
    #*u
    bigM[0:vLen, 0:vLen] = - Re*GRAD + beta*LAPLACIAN
    #*v
    bigM[0:vLen, vLen:2*vLen] = - Re*MMDYU
    #*w
    bigM[0:vLen, 2*vLen:3*vLen] = - Re*MMDZU
    #*p
    bigM[0:vLen, 3*vLen:4*vLen] = - 1.j*k*eye(vLen,vLen)
    #cxx
    bigM[0:vLen, 4*vLen:5*vLen] = (1-beta)*oneOverWeiss*1.j*k*eye(vLen,vLen)
    #cyy
    #czz
    #cxy
    bigM[0:vLen, 7*vLen:8*vLen] = (1-beta)*oneOverWeiss*MDY
    #cxz
    bigM[0:vLen, 8*vLen:9*vLen] = (1-beta)*oneOverWeiss*MDZ
    #cyz

    ################Navier Stokes y direction:###############
    #*u
    #*v
    bigM[vLen:2*vLen, vLen:2*vLen] = - Re*GRAD - Re*MMDYV \
                                     + beta*LAPLACIAN
                                     
    #*w
    bigM[vLen:2*vLen, 2*vLen:3*vLen] = - Re*MMDZV
    #*p
    bigM[vLen:2*vLen, 3*vLen:4*vLen] = - MDY
    #cxx
    #cyy
    bigM[vLen:2*vLen, 5*vLen:6*vLen] = (1-beta)*oneOverWeiss*MDY
    #czz
    #cxy
    bigM[vLen:2*vLen, 7*vLen:8*vLen] = (1-beta)*oneOverWeiss*1.j*k*eye(vLen,vLen)
    #cxz
    #cyz
    bigM[vLen:2*vLen, 9*vLen:10*vLen] = (1-beta)*oneOverWeiss*MDZ

    ################Navier Stokes z direction:###############
    #*u
    #*v
    bigM[2*vLen:3*vLen, vLen:2*vLen] = -Re*MMDYW
    #*w
    bigM[2*vLen:3*vLen, 2*vLen:3*vLen] = - Re*GRAD - Re*MMDZW \
                                         + beta*LAPLACIAN
    #*p
    bigM[2*vLen:3*vLen, 3*vLen:4*vLen] = - MDZ
    #cxx
    #cyy
    #czz
    bigM[2*vLen:3*vLen, 6*vLen:7*vLen] = (1-beta)*oneOverWeiss*MDZ
    #cxy
    #cxz
    bigM[2*vLen:3*vLen, 8*vLen:9*vLen] = (1-beta)*oneOverWeiss*1.j*k*eye(vLen,vLen)
    #cyz
    bigM[2*vLen:3*vLen, 9*vLen:10*vLen] = (1-beta)*oneOverWeiss*MDY

    ################Incompressability equation:###############
    #*u
    bigM[3*vLen:4*vLen, 0:vLen] = 1.j*k*eye(vLen,vLen)
    #*v
    bigM[3*vLen:4*vLen, vLen:2*vLen] = MDY
    #*w
    bigM[3*vLen:4*vLen, 2*vLen:3*vLen] = MDZ
    #*p
    #cxx
    #cyy
    #czz
    #cxy
    #cxz
    #cyz

    ################cxx equation:####################
    #*u
    bigM[4*vLen:5*vLen, 0:vLen] = 2.j*k*MMCXX + 2*dot(MMCXY,MDY) \
                                + 2*dot(MMCXZ,MDZ)
    #*v
    bigM[4*vLen:5*vLen, vLen:2*vLen] = -prod_mat(dot(MDY,conxx)) 
    #*w
    bigM[4*vLen:5*vLen, 2*vLen:3*vLen] = -prod_mat(dot(MDZ,conxx)) 
    #*p
    #cxx
    bigM[4*vLen:5*vLen, 4*vLen:5*vLen] = -oneOverWeiss*eye(vLen,vLen) - GRAD
    #cyy
    #czz
    #cxy
    bigM[4*vLen:5*vLen, 7*vLen:8*vLen] = 2*MMDYU
    #cxz
    bigM[4*vLen:5*vLen, 8*vLen:9*vLen] = 2*MMDZU
    #cyz
    
    ################cyy equation:####################
    #*u
    #*v
    bigM[5*vLen:6*vLen, vLen:2*vLen] = +2j*k*MMCXY +2*dot(MMCYY,MDY) \
                                     + 2*dot(MMCYZ,MDZ) \
                                     - prod_mat(dot(MDY,conyy))
    #*w
    bigM[5*vLen:6*vLen, 2*vLen:3*vLen] =  -prod_mat(dot(MDZ,conyy))
    #*p
    #cxx
    #cyy
    bigM[5*vLen:6*vLen, 5*vLen:6*vLen] = -oneOverWeiss*eye(vLen,vLen) - GRAD \
                                       + 2*MMDYV
    #czz
    #cxy
    #cxz
    #cyz
    bigM[5*vLen:6*vLen, 9*vLen:10*vLen] = 2*MMDZV

    ################czz equation:####################
    #*u
    #*v
    bigM[6*vLen:7*vLen, vLen:2*vLen] = -prod_mat(dot(MDY,conzz))
    #*w
    bigM[6*vLen:7*vLen, 2*vLen:3*vLen] =  -prod_mat(dot(MDZ,conzz)) + 2.j*k*MMCXZ \
                                       +  2*dot(MMCYZ,MDY) + 2*dot(MMCZZ,MDZ)
    #*p
    #cxx
    #cyy
    #czz
    bigM[6*vLen:7*vLen, 6*vLen:7*vLen] = - oneOverWeiss*eye(vLen,vLen) - GRAD \
                                       + 2*MMDZW 
    #cxy
    #cxz
    #cyz
    bigM[6*vLen:7*vLen, 9*vLen:10*vLen] = 2*MMDYW

    
    ################cxy equation:####################
    #*u
    bigM[7*vLen:8*vLen, 0:vLen] = + dot(MMCYY,MDY) + dot(MMCYZ,MDZ)
    #*v
    bigM[7*vLen:8*vLen, vLen:2*vLen] = -prod_mat(dot(MDY,conxy)) + 1.j*k*MMCXX \
                                       + dot(MMCXZ,MDZ)
    #*w
    bigM[7*vLen:8*vLen, 2*vLen:3*vLen] = - prod_mat(dot(MDZ,conxy)) \
                                        - dot(MMCXY,MDZ)
    #*p
    #cxx
    #cyy
    bigM[7*vLen:8*vLen, 5*vLen:6*vLen] =  MMDYU
    #czz
    #cxy
    bigM[7*vLen:8*vLen, 7*vLen:8*vLen] = -oneOverWeiss*eye(vLen,vLen) - GRAD\
                                         + MMDYV
    #cxz
    bigM[7*vLen:8*vLen, 8*vLen:9*vLen] =  MMDZV
    #cyz
    bigM[7*vLen:8*vLen, 9*vLen:10*vLen] =  MMDZU
    

    ################cxz equation:####################
    #*u
    bigM[8*vLen:9*vLen, 0:vLen] =  + dot(MMCYZ,MDY) + dot(MMCZZ,MDZ)
    #*v
    bigM[8*vLen:9*vLen, vLen:2*vLen] = - prod_mat(dot(MDY,conxz)) \
                                      - dot(MMCXZ,MDY)
    #*w
    bigM[8*vLen:9*vLen, 2*vLen:3*vLen] = - prod_mat(dot(MDZ,conxz)) + 1.j*k*MMCXX\
                                         + dot(MMCXY,MDY) 
    #*p
    #cxx
    #cyy
    #czz
    bigM[8*vLen:9*vLen, 6*vLen:7*vLen] = MMDZU
    #cxy
    bigM[8*vLen:9*vLen, 7*vLen:8*vLen] = MMDYW
    #cxz
    bigM[8*vLen:9*vLen, 8*vLen:9*vLen] = -oneOverWeiss*eye(vLen,vLen) - GRAD\
                                         + MMDZW
    #cyz
    bigM[8*vLen:9*vLen, 9*vLen:10*vLen] =  MMDYU

    ###############cyz equation:####################
    #*u
    bigM[9*vLen:10*vLen, 0:vLen] = -1.j*k*MMCYZ
    #*v
    bigM[9*vLen:10*vLen, vLen:2*vLen] = - prod_mat(dot(MDY,conyz)) + 1.j*k*MMCXZ \
                                        + dot(MMCZZ,MDZ)
    #*w
    bigM[9*vLen:10*vLen, 2*vLen:3*vLen] = - prod_mat(dot(MDZ,conyz)) + 1.j*k*MMCXY \
                                          + dot(MMCYY,MDY) 
    #*p
    #cxx
    #cyy
    bigM[9*vLen:10*vLen, 5*vLen:6*vLen] =  MMDYW
    #czz
    bigM[9*vLen:10*vLen, 6*vLen:7*vLen] =  MMDZV
    #cxy
    #cxz
    #cyz
    bigM[9*vLen:10*vLen, 9*vLen:10*vLen] = - oneOverWeiss*eye(vLen,vLen) - GRAD

    #Apply Boundary Conditions for u, v, w (dudy=v=dwdy=0 at y=+-1):
    #u
    for i in range((2*N+1)):
        bigM[M*(i+1)-2,:] = hstack((zeros(M*i), DERIVTOP, zeros(10*vLen-M*(i+1))))
        bigM[M*(i+1)-1,:] = hstack((zeros(M*i), DERIVBOT, zeros(10*vLen-M*(i+1))))
    del i
    #v
    for i in range((2*N+1), 2*(2*N+1)):
        bigM[M*(i+1)-2,:] = hstack((zeros(M*i), BTOP, zeros(10*vLen-M*(i+1))))
        bigM[M*(i+1)-1,:] = hstack((zeros(M*i), BBOT, zeros(10*vLen-M*(i+1))))
    del i
    #w
    for i in range(2*(2*N+1), 3*(2*N+1)):
        bigM[M*(i+1)-2,:] = hstack((zeros(M*i), DERIVTOP, zeros(10*vLen-M*(i+1))))
        bigM[M*(i+1)-1,:] = hstack((zeros(M*i), DERIVBOT, zeros(10*vLen-M*(i+1))))
    del i

    return bigM

#
# MAIN
#
#Start the clock:
startTime = time.time()

print """
----------------------------------------
N     = {0}
M     = {1}
Re    = {2}
beta  = {3}
Weiss = {4}
amp   = {5}
k     = {6}
----------------------------------------
""". format(N, M, Re, beta, Weiss, Amp, k)

# Unpickle the answer from part1 and the V and W vectors

fpickle = open('pf'+filename, 'r')
(U,V,W,conxx,conyy,conzz,conxy,conxz,conyz) = pickle.load(fpickle)
fpickle.close()

# Setup variables:
gamma = pi / piDivide 
zLength = 2.*pi/gamma
vecLen = M*(2*N+1)
vLen = vecLen
print vLen, type(vLen)
oneOverWeiss = 1./Weiss
# Set the oneOverC function: 1/2 for m=0, 1 elsewhere:
oneOverC = ones(M)
oneOverC[0] = 1. / 2.
#set up the CFunc function: 2 for m=0, 1 elsewhere:
CFunc = ones(M)
CFunc[0] = 2.

#Boundary arrays:
BTOP = ones(M)
BBOT = ones(M)
BBOT[1:M:2] = -1

singleDY = mk_single_diffy()
DERIVTOP = zeros((M), dtype='complex')
DERIVBOT = zeros((M), dtype='complex')
for j in range(M):
    DERIVTOP[j] = dot(BTOP, singleDY[:,j]) 
    DERIVBOT[j] = dot(BBOT, singleDY[:,j])
del j

# make some useful matrices
MDY = mk_diff_y()
MDZ = mk_diff_z()

MMU = prod_mat(U)
MMV = prod_mat(V)
MMW = prod_mat(W)
MMCXX = prod_mat(conxx)
MMCYY = prod_mat(conyy)
MMCZZ = prod_mat(conzz)
MMCXY = prod_mat(conxy)
MMCXZ = prod_mat(conxz)
MMCYZ = prod_mat(conyz)

MMDYU = prod_mat(dot(MDY, U))
MMDZU = prod_mat(dot(MDZ, U))
MMDYV = prod_mat(dot(MDY, V))
MMDZV = prod_mat(dot(MDZ, V))
MMDYW = prod_mat(dot(MDY, W))
MMDZW = prod_mat(dot(MDZ, W))


LAPLACIAN = -(k**2)*eye(vLen,vLen) + dot(MDY,MDY) + dot(MDZ,MDZ)
GRAD      = 1.j*k*MMU + dot(MMV,MDY) + dot(MMW,MDZ)

# Make the matrix for the generalised eigenvalue problem

equations_matrix= mk_bigM()

# Make the scaling matrix for RHS of equation
RHS = eye(10*vLen,10*vLen) 
RHS[:3*vLen, :3*vLen] = Re*eye(3*vLen,3*vLen)
# Zero all elements corresponding to p equation
RHS[3*vLen:4*vLen, :] = zeros((vLen,10*vLen))

# Apply boundary conditions to RHS
for i in range(3*(2*N+1)):
    RHS[M*(i+1)-1, M*(i+1)-1] = 0
    RHS[M*(i+1)-2, M*(i+1)-2] = 0
del i

if args.evecs:
    print 'finding eigenvectors and eigenvalues'
    # Use library function to solve for eigenvalues/vectors
    print 'in linalg.eig time=', (time.time() - startTime)
    eigenvals, evecs = linalg.eig(equations_matrix, RHS, overwrite_a=True)

    # Save output

    #make large_evs, as large as eigenvals, but only contains real part of large, 
    #physical eigenvalues. Rest are zeros. Index of large_evs same as that of 
    #eigenvalues 
    large_evs = zeros(len(eigenvals))
    for i in xrange(10*(2*N+1)*M):
        if (real(eigenvals[i]) > 0) and (real(eigenvals[i]) < 50):
            large_evs[i] = real(eigenvals[i])
    del i

    lead_index = argmax(large_evs)

    du =   evecs[           :(2*N+1)*M,   lead_index]
    dv =   evecs[  (2*N+1)*M:2*(2*N+1)*M, lead_index]
    dw =   evecs[2*(2*N+1)*M:3*(2*N+1)*M, lead_index]
    dp =   evecs[3*(2*N+1)*M:4*(2*N+1)*M, lead_index]
    dcxx = evecs[4*(2*N+1)*M:5*(2*N+1)*M, lead_index]
    dcyy = evecs[5*(2*N+1)*M:6*(2*N+1)*M, lead_index]
    dczz = evecs[6*(2*N+1)*M:7*(2*N+1)*M, lead_index]
    dcxy = evecs[7*(2*N+1)*M:8*(2*N+1)*M, lead_index]
    dcxz = evecs[8*(2*N+1)*M:9*(2*N+1)*M, lead_index]
    dcyz = evecs[9*(2*N+1)*M:10*(2*N+1)*M, lead_index]

    eigarray = vstack((real(eigenvals), imag(eigenvals))).T
    #remove nans and infs from eigenvalues
    #eigarray = eigarray[~isnan(eigarray).any(1), :]
    #eigarray = eigarray[~isinf(eigarray).any(1), :]

    print 'chosen eig: {e}'.format(e=lead_index)
    savetxt('ev-k{k}{fn}.dat'.format(k=k, fn=filename[:-7],
                                                     ), eigarray)

    pickle.dump((du,dv,dw,dp,dcxx,dcyy,dczz,dcxy,dcxz,dcyz), open('full-evecs-k{k}{fn}'.format(k=k,fn=filename), 'w'))

else:
    # eigenvalues only
    print 'finding eigenvalues only'
    # Use library function to solve for eigenvalues/vectors
    print 'in linalg.eig time=', (time.time() - startTime)
    eigenvals = linalg.eigvals(equations_matrix, RHS, overwrite_a=True)

    # Save output

    eigarray = vstack((real(eigenvals), imag(eigenvals))).T
    #remove nans and infs from eigenvalues
    eigarray = eigarray[~isnan(eigarray).any(1), :]
    eigarray = eigarray[~isinf(eigarray).any(1), :]

    savetxt('ev-k'+str(k)+filename[:-7]+'.dat', eigarray)

    #stop the clock
    print 'done in', (time.time()-startTime)

######################TESTS####################################################

###############################################################################
