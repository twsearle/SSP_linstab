#----------------------------------------------------------------------------#
#   Fully Spectral Linear stability analysis Oldroyd B Model
#   Last modified: Tue 25 Feb 15:18:50 2014
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

#set parameters 
#-----------------------------------
N       = 5
M       = 40
Re      = 0.01
beta    = 0.1
Weiss   = 10.0 
Amp     = 0.02
gamma   = pi / 2
k       = float(sys.argv[1])
#-----------------------------------

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

    #Apply Boundary Conditions for u, v, w:
    for i in range(3*(2*N+1)):
        bigM[M*(i+1)-2,:] = hstack((zeros(M*i), BTOP, zeros(10*vLen-M*(i+1))))
        bigM[M*(i+1)-1,:] = hstack((zeros(M*i), BBOT, zeros(10*vLen-M*(i+1))))
    del i

    return bigM

#
# MAIN
#


#Start the clock:
startTime = time.time()

filename = '-N{N}-M{M}-Re{Re}-b{beta}-Wi{Weiss}-amp{Amp}.pickle'.format(\
            N=N,M=M,Re=Re,beta=beta,Weiss=Weiss,Amp=Amp)

# Unpickle the answer from part1 and the V and W vectors

fpickle = open('pf'+filename, 'r')
(U,V,W,conxx,conyy,conzz,conxy,conxz,conyz) = pickle.load(fpickle)
fpickle.close()

# Setup variables:
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

# Use library function to solve for eigenvalues/vectors
print 'In linalg.eig time =', (time.time() - startTime)
eigenvals = linalg.eig(equations_matrix, RHS, left=False, right=False,
                       overwrite_a=True, overwrite_b=True)

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
