#--------------------------------------------------------------------------------------------------
#   Time Iteration Program for the Linear Stability Analysis of the Viscoelastic
#   SSP
#   timeIteration.py
# 
#   Last modified: Tue 10 Sep 16:47:54 2013
#--------------------------------------------------------------------------------------------------

""" Time Iteration Program for the Linear Stability Analysis of the Viscoelastic SSP
Loop over:
    - calculate stress disturbances from conformation disturbances
    - calculate stream-functions disturbances from stress disturbances
    - calculate disturbance velocities from stream-function disturbance
    velocities
    - calculate new conformations from base and disturbance velocities
    - save the size of the disturbance
"""

# MODULES
import sys
import time
from scipy import *
from scipy import linalg
from scipy.sparse import linalg as sparselinalg
import cPickle as pickle
import ConfigParser

# PARAMETERS
startTime = time.time()

config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
N = config.getint('settings', 'N')
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')
kx = config.getfloat('settings', 'k')

fp.close()

filename = '-N{N}-M{M}-Re{Re}-b{beta}-Wi{Weiss}-amp{Amp}.pickle'.format(\
            N=N,M=M,Re=Re,beta=beta,Weiss=Weiss,Amp=Amp)

# FUNCTIONS

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

#
# MAIN
#

#fpickle = open('pf'+filename, 'r')
#(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(fpickle)
#fpickle.close()

# Setup variables:
gamma = pi / 2.
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

#### MAKE OPERATORS

# basic operators
MDY = mk_diff_y()
MDYY = dot(MDY,MDY)
MDZ = mk_diff_z()
MDZZ = dot(MDZ,MDZ)

LAPLACIAN = -(kx**2)*eye((2*N+1)*M, (2*N+1)*M) + MDYY + MDZZ

singleMDY = mk_single_diffy()
singleMDYY = dot(singleMDY, singleMDY)

# Boundary arrays
BTOP = ones(M)
BBOT = ones(M)
BBOT[1:M:2] = -1

DBTOP = ones(M)
DBBOT = ones(M)
for m in range(M):
    DBTOP[m] = dot(BTOP, singleMDY[:,m])
    DBBOT[m] = dot(BBOT, singleMDY[:,m])
del m

# The phi and psi streamfunction operator lists 
phiOps = []
psiOps = []
tempMat = zeros((M,M), dtype='complex')

for i in range(2*N+1):
    n = i-N
    # Solving these equations in one go because they are only 2nd order in fact.
    # remember, need to use the identity matrix for scalars except when
    # multiplying (since * is an elementwise multiplication)

    singleLaplacian = (-kx**2  - (n*gamma)**2)*eye(M,M) + singleMDYY
    tempPhiMat = dot(singleLaplacian, singleLaplacian) * (-kx**2 - (n*gamma)**2)
    tempPsiMat = singleLaplacian * (-kx**2 -(n*gamma)**2) 

    # apply boundary conditions to the operators before inverting (DON'T FORGET TO DO THIS TO LHS)
    # The order is: top, bottom, dtop, dbot
    tempPhiMat[M-4,:] = BTOP
    tempPhiMat[M-3,:] = BBOT
    tempPhiMat[M-2,:] = DBTOP
    tempPhiMat[M-1,:] = DBBOT
    
    tempPsiMat[M-2,:] = BTOP
    tempPsiMat[M-1,:] = BBOT

    tempPhiMat = linalg.inv(tempPhiMat)
    tempPsiMat = linalg.inv(tempPsiMat)

    phiOps.append(tempPhiMat)
    psiOps.append(tempPsiMat)
del i, n


