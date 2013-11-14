#--------------------------------------------------------------------------------------------------
#   Time Iteration Program for the Linear Stability Analysis of the Viscoelastic
#   SSP
#   timeIteration.py
# 
#   Last modified: Wed 13 Nov 22:10:56 2013
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

# PARAMETERS -----------------------------------------------------------------
startTime = time.time()

config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
N = config.getint('settings', 'N')
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Wi = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')
kx = config.getfloat('settings', 'k')

fp.close()

TIMESTEP = 0.0001
numTimeSteps = 10

filename = '-N{N}-M{M}-Re{Re}-b{beta}-Wi{Wi}-amp{Amp}.pickle'.format(\
            N=N,M=M,Re=Re,beta=beta,Wi=Wi,Amp=Amp)

# -----------------------------------------------------------------------------
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

#### SETUP VARIABLES ####

gamma = pi / 2.
zLength = 2.*pi/gamma
vecLen = M*(2*N+1)
vLen = vecLen
print vLen, type(vLen)
oneOverWi = 1./Wi
# Set the oneOverC function: 1/2 for m=0, 1 elsewhere:
oneOverC = ones(M)
oneOverC[0] = 1. / 2.
#set up the CFunc function: 2 for m=0, 1 elsewhere:
CFunc = ones(M)
CFunc[0] = 2.

# make psi and phi vectors
phi = zeros(((2*N+1)*M), dtype='complex')
psi = zeros(((2*N+1)*M), dtype='complex')

# Read in the base profile
(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(open('pf'+filename, 'r'))


#### MAKE OPERATORS ####

# basic operators
II = eye((2*N+1)*M,(2*N+1)*M)
MDY = mk_diff_y()
MDYY = dot(MDY,MDY)
MDZ = mk_diff_z()
MDZZ = dot(MDZ,MDZ)
MDYZ = dot(MDY,MDZ)
MDYZZ = dot(MDY,MDZZ)
MDYYZ = dot(MDYY,MDZ)
MDZZZ = dot(MDZ,MDZZ)


singleMDY = mk_single_diffy()
singleMDYY = dot(singleMDY, singleMDY)

# vectors
MMU = prod_mat(U)
MMV = prod_mat(V)
MMW = prod_mat(W)
MMCXX = prod_mat(Cxx)
MMCYY = prod_mat(Cyy)
MMCZZ = prod_mat(Czz)
MMCXY = prod_mat(Cxy)
MMCXZ = prod_mat(Cxz)
MMCYZ = prod_mat(Cyz)
 
MMDYU = prod_mat(dot(MDY,U))
MMDYV = prod_mat(dot(MDY,V))
MMDYW = prod_mat(dot(MDY,W))
MMDYCXX = prod_mat(dot(MDY,Cxx))
MMDYCYY = prod_mat(dot(MDY,Cyy))
MMDYCZZ = prod_mat(dot(MDY,Czz))
MMDYCXY = prod_mat(dot(MDY,Cxy))
MMDYCXZ = prod_mat(dot(MDY,Cxz))
MMDYCYZ = prod_mat(dot(MDY,Cyz))

MMDZU = prod_mat(dot(MDZ,U))
MMDZV = prod_mat(dot(MDZ,V))
MMDZW = prod_mat(dot(MDZ,W))
MMDZCXX = prod_mat(dot(MDZ,Cxx))
MMDZCYY = prod_mat(dot(MDZ,Cyy))
MMDZCZZ = prod_mat(dot(MDZ,Czz))
MMDZCXY = prod_mat(dot(MDZ,Cxy))
MMDZCXZ = prod_mat(dot(MDZ,Cxz))
MMDZCYZ = prod_mat(dot(MDZ,Cyz))

# Complex operators

LAPLACIAN = -(kx**2)*II + MDYY + MDZZ
VGRAD = 1.j*kx*II*U + dot(MMV,MDY) + dot(MMW,MDZ) 

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
    # Solving these equations in one go because they are only 2nd order.
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

    # use solve or invert?
    #tempPhiMat = linalg.inv(tempPhiMat)
    #tempPsiMat = linalg.inv(tempPsiMat)

    phiOps.append(tempPhiMat)
    psiOps.append(tempPsiMat)
del i, n

# The zeroth mode operators
dw0thOp = -kx**2*eye(M,M) + singleMDYY
dv0thOp = dot((-kx**2*eye(M,M) + singleMDYY), (eye(M,M) +
                                               (1./kx**2)*singleMDYY))

# The zeroth mode boundary conditions
dw0thOp[M-2,:] = BTOP
dw0thOp[M-1,:] = BBOT

dv0thOp[M-2,:] = BTOP
dv0thOp[M-1,:] = BBOT

# use solve or invert?

## SET UP A GUESS FOR THE INITIAL CONFORMATION PROFILE
cvec = zeros((6*vecLen),dtype='complex')
cvec = (random.random(6*vecLen) + 1.j*random.random(6*vecLen))/1000000.0

outFile = open("strength{0}".format(filename), 'w')

#### ITERATE THE INSTABILITY ####
for currTime in range(numTimeSteps):

    dcxx = cvec[:1*vecLen]
    dcyy = cvec[vecLen:2*vecLen]
    dczz = cvec[2*vecLen:3*vecLen]
    dcxy = cvec[3*vecLen:4*vecLen]
    dcxz = cvec[4*vecLen:5*vecLen]
    dcyz = cvec[5*vecLen:6*vecLen]

    # As a test, use only the stresses to recreate a disturbance.
    # need to make a reasonable disturbance using Oldroyd-B2-full.py

    # (du_test,dv_test,dw_test,dp_test,dcxx,dcyy,dczz,dcxy,dcxz,dcyz) = pickle.load(
    #                                     open('lead-evec-k0.7'+filename,'r'))
    # print "magnitude of test disturbance ", linalg.norm(du_test,2)
    # print "magnitude of base profile ", linalg.norm(U, 2)

    ### MAKE STREAM FUNCTIONS

    # Make the vectors from stresses to calculate phi and psi
    # Stress disturbances are equal to confirmation disturbances in Oldroyd-B

    phiStressVec = - kx**2*dot(MDY,dcxx) \
                   + dot( (kx**2*MDY  - MDYZZ), dcyy) + dot(MDYZZ,dczz) \
                   + dot( ((1.j*kx**3)*II - 1.j*kx*MDZZ \
                          + 1.j*kx*MDYY), dcxy) \
                   + dot( 2.j*kx*MDYZ, dcxz) \
                   + dot( (kx**2*MDZ + MDYYZ - MDZZZ), dcyz) 
    phiStressVec = -((1.-beta)/(beta*Wi))*phiStressVec

    psiStressVec = dot(1.j*kx*MDZ, dcxx) + dot(MDYZ,dcxy) \
                   + dot((MDZZ + kx**2*II), dcxz) \
                   - dot(1.j*kx*MDY, dcyz) - dot(1.j*kx*MDZ, dczz)
    psiStressVec = ((1-beta)/(beta*Wi))*psiStressVec

    # insert boundary conditions

    phiStressVec[M-4:(2*N+1)*M:M] = 0
    phiStressVec[M-3:(2*N+1)*M:M] = 0
    phiStressVec[M-2:(2*N+1)*M:M] = 0
    phiStressVec[M-1:(2*N+1)*M:M] = 0

    psiStressVec[M-2:(2*N+1)*M:M] = 0
    psiStressVec[M-1:(2*N+1)*M:M] = 0

    # solve for stream functions

    for i in range(2*N+1):
        # go through each operator and z mode and generate the corresponding z mode
        # of phi and psi
        if i is N: continue #miss out the zeroth mode
        phi[i*M:(i+1)*M] = linalg.solve(phiOps[i], phiStressVec[i*M:(i+1)*M])
        psi[i*M:(i+1)*M] = linalg.solve(psiOps[i], psiStressVec[i*M:(i+1)*M])
    del i

    ### FIND VELOCITY AND PRESSURE DISTURBANCES

    #Find most of the modes 
    du = dot(1.j*kx*MDY, phi) - dot(MDZ, psi)
    dv = kx**2*phi - dot(MDZZ, phi)
    dw = dot(MDYZ, phi) + 1.j*kx*psi

    # calculate dw0 
    dw0vec = - (1.-beta)/(beta*Wi) * (kx*1.j*dcxz[N*M:(N+1)*M] + \
                                      dot(singleMDY, dcyz[N*M:(N+1)*M]))
                                  
    # put in boundary conditions
    dw0vec[M-2] = 0 
    dw0vec[M-1] = 0

    dw[N*M:(N+1)*M] = linalg.solve(dw0thOp, dw0vec)

    # calculate dv0
    dv0vec = (1.-beta)/(beta*Wi)*( dot(singleMDY,dcxx[N*M:(N+1)*M]) \
                                  + dot( ((1./(1.j*kx))*singleMDYY - 1.j*kx), \
                                        dcxy[N*M:(N+1)*M]) \
                                  - dot( singleMDY, dcyy[N*M:(N+1)*M] ) )
    dv0vec[M-2] = 0 
    dv0vec[M-1] = 0

    dv[N*M:(N+1)*M] =  linalg.solve(dv0thOp, dv0vec)

    # calculate du0
    du[N*M:(N+1)*M] = dot( (1./(1.j*kx))*singleMDY, dv[N*M:(N+1)*M] )

    # calculate pressure
    dp = (beta/(1.j*kx))*dot(LAPLACIAN, du) \
        + (1.-beta)/(Wi)*( dcxx + (1./(1.j*kx))*dot(MDY,dcxy) )


    # ---TESTS-----------------------------------------------------------------------
    #print '==================TESTS======================='
    #print 'zeroth modes'
    #print 'du0 has been calculated correctly: ', allclose(du_test[N*M:(N+1)*M], du[N*M:(N+1)*M])
    #print 'dv0 has been calculated correctly: ', allclose(dv_test[N*M:(N+1)*M], dv[N*M:(N+1)*M])
    #print 'dw0 has been calculated correctly: ', allclose(dw_test[N*M:(N+1)*M], dw[N*M:(N+1)*M])
    #print 'dp0 has been calculated correctly: ', allclose(dp_test[N*M:(N+1)*M], dp[N*M:(N+1)*M])
    #print '=============================================='
    #print 'du has been calculated correctly: ', allclose(du_test, du)
    #print 'dv has been calculated correctly: ', allclose(dv_test, dv)
    #print 'dw has been calculated correctly: ', allclose(dw_test, dw)
    #print 'dp has been calculated correctly: ', allclose(dp_test, dp)

    #duData = vstack((du_test[N*M:(N+1)*M], du[N*M:(N+1)*M])).T
    #dvData = vstack((dv_test[N*M:(N+1)*M], dv[N*M:(N+1)*M])).T
    #dwData = vstack((dw_test[N*M:(N+1)*M], dw[N*M:(N+1)*M])).T

    #savetxt("du_comp.dat", (duData))
    #savetxt("dv_comp.dat", (dvData))
    #savetxt("dw_comp.dat", (dwData))

    #print "exit before time step."
    #exit(1)
    # -------------------------------------------------------------------------------

    ### TAKE A TIME STEP IN THE CONFORMATION VECTORS

    # make a vector containing the results of the RHS of the iteration equations.

    iterVec = zeros(6*vecLen, dtype='complex')

    #  cxx
    iterVec[:vecLen] = - oneOverWi*dcxx - dot(VGRAD, dcxx) \
                     - dot(MMDYCXX, dv) - dot(MMDZCXX, dw) + 2.j*kx*dot(MMCXX, du)\
                     + 2*dot(dot(MMCXY, MDY), du) + 2*dot( dot(MMCXZ, MDZ), du) \
                     + 2*dot(MMDYU, dcxy) + 2*dot(MMDZU, dcxz)
    # cyy
    iterVec[1*vecLen:2*vecLen] = - oneOverWi*dcyy - dot(VGRAD, dcyy) \
                     - dot(MMDYCYY, dv) - dot(MMDZCYY, dw) + 2.j*kx*dot(MMCXY, dv)\
                     + 2*dot(dot(MMCYY, MDY), dv)  + 2*dot(dot(MMCYZ, MDZ), dv) \
                     + 2*dot(MMDYV, dcyy) + 2*dot(MMDZV, dcyz)

                            
    # czz
    iterVec[2*vecLen:3*vecLen] = - oneOverWi*dczz - dot(VGRAD, dczz) \
                     - dot(MMDYCZZ, dv) - dot(MMDZCZZ, dw) + 2.j*kx*dot(MMCXZ, dv)\
                     + 2*dot(dot(MMCYZ, MDY), dw) + 2*dot(dot(MMCZZ, MDZ), dw) \
                     + 2*dot(MMDYW, dcyz) + 2*dot(MMDZW, dczz)

    # cxy
    iterVec[3*vecLen:4*vecLen] = - oneOverWi*dcxy - dot(VGRAD, dcxy) \
                     - dot(MMDYCXY, dv) - dot(MMDZCXY, dw) \
                     + dot(dot(MMCYY, MDY) + dot(MMCYZ, MDZ), du) \
                     + dot(MMDYU, dcyy) + dot(MMDZU, dcyz) + 1.j*kx*dot(MMCXX, dv)\
                     + dot(dot(MMCXZ, MDZ), dv) + dot(MMDYV, dcxy) \
                     + dot(MMDZV, dcxz) - dot(dot(MMCXY, MDZ), dw)

    # cxz
    iterVec[4*vecLen:5*vecLen] = - oneOverWi*dcxz - dot(VGRAD, dcxz) \
                     - dot(MMDYCXZ, dv) - dot(MMDZCXZ, dw) \
                     + dot(dot(MMCYZ, MDY), du) + dot(dot(MMCZZ, MDZ), du) \
                     + dot(MMDYU, dcyz) + dot(MMDZU, dczz) + 1.j*dot(MMCXX, dw) \
                     + dot(dot(MMCXY, MDY), dw) + dot(MMDYW, dcxy) \
                     + dot(MMDZW, dcxz) - dot(dot(MMCXZ, MDY), dv)

    # cyz
    iterVec[4*vecLen:5*vecLen] = - oneOverWi*dcyz - dot(VGRAD, dcyz) \
                     - dot(MMDYCYZ, dv) - dot(MMDZCYZ, dw) \
                     + 1.j*dot(MMCXZ, dv) + dot(dot(MMCZZ, MDZ), dv) \
                     + dot(MMDZV, dczz) + 1.j*kx*dot(MMCXY, dw) \
                     + dot(dot(MMCYY, MDY), dw) + dot(MMDYW, dcyy) \
                     - 1.j*kx*dot(MMCYZ, du)

    # do an iteration of the conformations (Euler for now)
    cnew = cvec + TIMESTEP*iterVec
    cvec = cnew

    # save the norm of the disturbance strength and the time to a file
    strength = linalg.norm(concatenate((du,dv,dw,cvec)),2)
    outFile.write('{0:10.5g} {1:10.5g}'.format(currTime, strength))
    outFile.flush()
del curr_time

outFile.close()

# Tidy up and dump the final disturbance vector to a file
dcxx = cvec[:1*vecLen]
dcyy = cvec[vecLen:2*vecLen]
dczz = cvec[2*vecLen:3*vecLen]
dcxy = cvec[3*vecLen:4*vecLen]
dcxz = cvec[4*vecLen:5*vecLen]
dcyz = cvec[5*vecLen:6*vecLen]

output = (du, dv, dw, dp, dcxx, dcyy, dczz, dcxy, dcxz, dcyz)
pickle.dump(output, open('final-dist-vec{0}'.format(filename), 'w'))

print "DONE!"
