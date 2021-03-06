#----------------------------------------------------------------------------#
#   Fully Spectral Newton Raphson Solver                            
#   Oldroyd B Model
#   Last modified: Fri 21 Mar 2014 15:29:22 GMT
#----------------------------------------------------------------------------#
"""Solves system of equations using a fully spectral method. Equations given 
by: V.dU(y,z)/dy + W.dU/dz = 1/Re .del^2."""

# MODULES
import sys
from scipy import *
from scipy import linalg
from scipy import optimize
from scipy import special
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
argparser.add_argument("-piDivide", type=float, default=cfgPiDivide,
                help='Override piDivide from the config file')

args = argparser.parse_args()
N = args.N 
M = args.M
Re = args.Re
beta = args.b
Weiss = args.Wi
Amp = args.amp
piDivide = args.piDivide


filename = 'pf-N{N}-M{M}-Re{Re}-b{beta}-Wi{Wi}-amp{Amp}-gdiv{gdiv}.pickle'.format(\
            N=N,M=M,Re=Re,beta=beta,Wi=Weiss,Amp=Amp, gdiv=piDivide)

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

def mk_no_slip_V():
    V = zeros((M, 2*N+1), dtype = 'complex')
    for m in range(0,M,2):
        V[m,N-1] = 2*oneOverC[m]*( ((-1)**(m/2))*(special.jv(m,p)/cos(p)) - 
                    special.iv(m,gamma)/cosh(gamma) )
        V[m,N+1] = 2*oneOverC[m]*( ((-1)**(m/2))*(special.jv(m,p)/cos(p)) - 
                    special.iv(m,gamma)/cosh(gamma) )
    del m        
    V = 0.5*V       #For the cosine amplitude.
    Normal = ( cos(p)*cosh(gamma) ) / ( cosh(gamma) - cos(p) )
    V = Amp * Normal * V
    return V.T.flatten() #return 1D array

def mk_no_slip_W():
    W = zeros((M, 2*N+1), dtype = 'complex')
    for m in range(0,M,2):
        W[m,N-1] = 2*oneOverC[m]*( ((-1)**(m/2))*(special.jv(m,p)/cos(p)) - 
                    special.iv(m,gamma)/cosh(gamma) )
        W[m,N+1] = 2*oneOverC[m]*( ((-1)**(m/2))*(special.jv(m,p)/cos(p)) - 
                    special.iv(m,gamma)/cosh(gamma) )
    del m
    
    chebdY = mk_single_diffy()
    W[:,N-1] = -dot(chebdY, W[:,N-1])
    W[:,N+1] = -dot(chebdY, W[:,N+1])
    
    W[:,N-1] = W[:,N-1]*0.5j
    W[:,N+1] = W[:,N+1]*-0.5j
    Normal = ( cos(p)*cosh(gamma) ) / ( cosh(gamma) - cos(p) )
    W = Amp * Normal * W / gamma
    return W.T.flatten() #return 1D array

def mk_free_slip_U():
    """
    Create guess of the free slip version of the streamwise base velocity.
    U(y) = sin bb*y
    """

    U = zeros((M*(2*N+1)), dtype='D')

    # Regardless of gamma. cf Waleffe
    bb = pi/2.      

    U = zeros(M*(2*N+1) ,dtype='D')

    for m in range(1,M,2):
        U[M*N+m] = 2*oneOverC[m]*((-1)**((m-1)/2))* special.jv(m,pi/2.)
    del m 

    return U

def mk_free_slip_V():
    """
    Create the free slip version of the wall-normal base velocity.
    V(y) = cos bb*y * cos gamma z
    """

    V = zeros((M*(2*N+1)), dtype='D')

    # Regardless of gamma. cf Waleffe
    bb = pi/2.      

    # Set the y dependence


    for m in range(0,M,2):
        # Only non-zero for even m since it is an even function
        V[M*(N-1)+m] = 2*oneOverC[m]*((-1)**(m/2)) * special.jv(m,bb)
    del m


    # Set the z dependence

    V[M*(N+1): M*(N+2)] = V[M*(N-1):M*N]

    return Amp*V

def mk_free_slip_W():

    V = mk_free_slip_V()/Amp
    W = zeros((M*(2*N+1)), dtype='D')

    chebDY = mk_single_diffy()

    W[M*(N-1): M*N] = 0.5j*(-dot(chebDY, V[M*(N-1): M*N])) 
    W[M*(N+1): M*(N+2)] = -0.5j*(-dot(chebDY, V[M*(N-1): M*N])) 

    return Amp*W / gamma

def mk_pressure_gradient():
    """
    The free slip pressure gradient which sustains the flow
    dPdx = - bb^2 sin(bb * y)
    """

    dPdx = zeros(M*(2*N+1) ,dtype='D')

    for m in range(1,M,2):
        dPdx[M*N+m] = 2*oneOverC[m]*((-1)**((m-1)/2))* special.jv(m,pi/2.)
    del m 

    return -((pi**2) / 4.)*dPdx


def solve_eq1():
    """Oldroyd-B Equation independent of U
    Linearly solves system of equations for vector containing Cyy, Czz, Cyz"""

    #RHS of equation
    RHS = zeros(3*vecLen,dtype='D')
    RHS[N*M] = oneOverWeiss
    RHS[vecLen + N*M] = oneOverWeiss

    # Make the Jacobian:
    jacobian = zeros((3*vecLen, 3*vecLen), dtype='complex')
    
    # (GRAD - 2dVdy + I/Wi)*dCyy
    jacobian[0:vecLen, 0:vecLen] = \
        GRAD - 2*MMDYV + oneOverWeiss*eye(vecLen, vecLen)
    # 0*dCzz
    # (-2dVdz)*dCyz
    jacobian[0:vecLen, 2*vecLen:3*vecLen] = -2*MMDZV
    
    # 2nd block of rows
    # 0*dCyy
    # (GRAD - 2*MMDZW + I/Wi)*dCzz
    jacobian[vecLen:2*vecLen, vecLen:2*vecLen] = \
        GRAD - 2*MMDZW + oneOverWeiss*eye(vecLen, vecLen)
    # (-2*dWdy)*dCyz
    jacobian[vecLen:2*vecLen, 2*vecLen:3*vecLen] = -2*MMDYW

    # 3rd block of rows
    # (-dWdy)*dCyy
    jacobian[2*vecLen:3*vecLen, 0:vecLen] = -MMDYW
    # (-dVdz)*dCzz
    jacobian[2*vecLen:3*vecLen, vecLen:2*vecLen] = -MMDZV
    # (GRAD + I/Wi)*dCyz
    jacobian[2*vecLen:3*vecLen, 2*vecLen:3*vecLen] = \
        GRAD + oneOverWeiss*eye(vecLen, vecLen)

    return linalg.solve(jacobian, RHS)

def solve_eq2(x):
    """use functions to find the residuals vector for the equation containing U"""
    #cut x into U, and the Conformation tensor arrays
    U = (x[0:vecLen])
    conxx = (x[vecLen: 2*vecLen])
    conxy = (x[2*vecLen: 3*vecLen])
    conxz = (x[3*vecLen: 4*vecLen])

    #calculate the stress tensor components, which are used to solve
    #the equations

    tauxx = oneOverWeiss*conxx 
    tauxx[N*M] -= oneOverWeiss
    tauxy = oneOverWeiss*conxy
    tauxz = oneOverWeiss*conxz

    #some useful matrices to have:
    MMDYU = prod_mat(dot(MDY, U))
    MMDZU = prod_mat(dot(MDZ, U))

    #calculate the equations
    #xx 0 = v.grad(Cxx) - 2(Cxy.dy.U + CxzdzU) + tauxx
    resxx = (dot(GRAD, conxx)
        - 2*dot(MMDYU, conxy)
        - 2*dot(MMDZU, conxz)
        + tauxx) 

    # xy 0 = (v.grad)Cxy - CyydyU  - CxydyV - CxzdzV - CyzdzU + tauxy"""
    resxy = ( dot(GRAD, conxy)
        - dot(MMDYU, conyy)
        - dot(MMDZU, conyz)
        - dot(MMDYV, conxy)
        - dot(MMDZV, conxz)
        + tauxy )

    # xz 0 = (v.grad)Cxz - CxydyW - CxzDzW - CyzdyU - CzzdzU + tauxz"""
    resxz = ( dot(GRAD, conxz)
        - dot(MMDYW, conxy)
        - dot(MMDZW, conxz)
        - dot(MMDYU, conyz)
        - dot(MMDZU, conzz)
        + tauxz )

    # Navier-Stokes 0 = -Re*(v.grad)U + beta.laplacian + (1-beta).div.tau"""
    resNS = ( -Re*dot(GRAD, U) + beta*dot(LAPLACIAN, U)
            + (1-beta)*(dot(MDY, tauxy) + dot(MDZ,tauxz))
            - dPdx )

    # Impose Boundary condition equations
    # Free slip so we need dUdy = 0 everywhere, and u_0,0 (the y,z independent
    # part) must be zero?

    for k in range (2*N+1):
        resNS[k*M + M-2] = dot(DERIVTOP, U[k*M:k*M + M])
        resNS[k*M + M-1] = dot(DERIVBOT, U[k*M:k*M + M])
    del k

    resNS[N*M] = 0

    # make residuals vector in 1D again
    residuals = zeros(4*vecLen, dtype='complex')
    residuals[0:vecLen] = resNS
    residuals[vecLen:2*vecLen] = resxx
    residuals[2*vecLen:3*vecLen] = resxy
    residuals[3*vecLen:4*vecLen] = resxz

    #make the Jacobian:
    jacobian = zeros((4*vecLen,4*vecLen), dtype='complex')

    #First row of blocks, for U equation:
    #-Re*(Grad)dU
    jacobian[0:vecLen, 0:vecLen] = -Re*GRAD + beta*LAPLACIAN
    #0*dCxx
    #((1-beta)1/Wi*d/dy)dCxy
    jacobian[0:vecLen, 2*vecLen:3*vecLen] = \
        (1-beta)*oneOverWeiss*MDY
    #((1-beta)1/wi*d/dz)dCxz
    jacobian[0:vecLen, 3*vecLen:4*vecLen] = \
        (1-beta)*oneOverWeiss*MDZ

    #2nd Row of blocks, for conxx equation:
    #(-2conxy*d/dy -2conxz*d/dz)dU
    jacobian[vecLen:2*vecLen, 0:vecLen] = \
        -2*( dot(prod_mat(conxy), MDY) + dot(prod_mat(conxz), MDZ) )
    #(GRAD + I/Wi)dCxx
    jacobian[vecLen:2*vecLen, vecLen:2*vecLen] = \
        GRAD + oneOverWeiss*eye(vecLen,vecLen, dtype='complex')
    #(-2dUdy)dCxy
    jacobian[vecLen:2*vecLen, 2*vecLen:3*vecLen] = \
        -2*( MMDYU )
    #(-2dUdz)dCxz
    jacobian[vecLen:2*vecLen, 3*vecLen:4*vecLen] = \
        -2*( MMDZU )

    #Third row of blocks, for conxy equation:
    #(-Cyy*d/dy -Cyz*d/dz)dU
    jacobian[2*vecLen:3*vecLen, 0:vecLen] = \
        - dot(prod_mat(conyy), MDY) - dot(prod_mat(conyz), MDZ)
    #0*dCxx
    #(GRAD - dVdy)*dCxy
    jacobian[2*vecLen:3*vecLen, 2*vecLen:3*vecLen] = \
        ( GRAD - MMDYV + oneOverWeiss*eye(vecLen,vecLen, dtype='complex') )
    #(-dV/dz)*dCxz
    jacobian[2*vecLen:3*vecLen, 3*vecLen:4*vecLen] = \
        -MMDZV

    #4th row of blocks, conxz equation:
    #(-Czz*d/dz -Cyz*d/dy)*dU
    jacobian[3*vecLen:4*vecLen, 0:vecLen] = \
        -dot(prod_mat(conzz), MDZ) - dot(prod_mat(conyz), MDY)
    #(0)*dCxx
    #(-dW/dy)*dCxy
    jacobian[3*vecLen:4*vecLen, 2*vecLen:3*vecLen] = \
        -MMDYW
    #(Grad - dW/dz)*dCxz
    jacobian[3*vecLen:4*vecLen, 3*vecLen:4*vecLen] = \
        ( GRAD - MMDZW + oneOverWeiss*eye(vecLen,vecLen, dtype='complex') )

    #Apply BC's in the Jacobian
    for n in range(2*N+1):
        jacobian[n*M + M-2, 0 : 4*vecLen ] = \
            concatenate( (zeros(n*M), DERIVTOP, zeros((2*N-n)*M+3*vecLen)) )
        jacobian[n*M + M-1, 0 : 4*vecLen ] = \
            concatenate( (zeros(n*M), DERIVBOT, zeros((2*N-n)*M+3*vecLen)) )
    del n

    # Free slip requires an extra restriction on u_0,0 (the constant flow in
    # y,z) u_0,0 = 0?

    jacobian[N*M, :] = 0
    jacobian[N*M, N*M] = 1


    #my_mc.matrix_checker(jacobian, vecLen, False)

    return (jacobian, residuals)

def NR_solve(the_eq, xguess):
    """use Newton-Raphson to solve """
    while True:  
        (J_x0, f_x0)  = the_eq(xguess)
        dx = linalg.solve(J_x0, -f_x0)
        xguess = xguess + dx
        print linalg.norm(f_x0, 2)
        if (linalg.norm(f_x0,2) < NRdelta): break

    return xguess

def save_pickle(array, name):
    f = open(name, 'w')
    pickle.dump(array, f)
    f.close()
#
# MAIN
#

print """
----------------------------------------
N     = {0}
M     = {1}
Re    = {2}
beta  = {3}
Weiss = {4}
amp   = {5}
gamma = {6}
----------------------------------------
""". format(N, M, Re, beta, Weiss, Amp, pi/piDivide)

gamma = pi / piDivide
p = optimize.fsolve(lambda p: p*tan(p) + gamma*tanh(gamma), 2)
zLength = 2.*pi/gamma
vecLen = M*(2*N+1)
NRdelta = 0.00001
# Set the oneOverC function: 1/2 for m=0, 1 elsewhere:
oneOverC = ones(M)
oneOverC[0] = 1. / 2.
#set up the CFunc function: 2 for m=0, 1 elsewhere:
CFunc = ones(M)
CFunc[0] = 2.

#make the differentiation matrices:
MDY = mk_diff_y()
MDYY = dot(MDY,MDY)
MDZ = mk_diff_z()
MDZZ = dot(MDZ,MDZ)

V = mk_free_slip_V()
W = mk_free_slip_W()
U = zeros(vecLen, dtype='D')
dPdx = mk_pressure_gradient()
pickle.dump(dPdx, open('dPdx.pickle', 'w'))

GRAD = dot(prod_mat(V),MDY) + dot(prod_mat(W),MDZ)
LAPLACIAN = dot(MDY,MDY)  + dot(MDZ,MDZ)
MMDYV = prod_mat(dot(MDY, V))
MMDZV = prod_mat(dot(MDZ, V))
MMDYW = prod_mat(dot(MDY, W))
MMDZW = prod_mat(dot(MDZ, W))


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

oneOverWeiss = 1. / Weiss


#Guess conformation tensor:
conxx = zeros(vecLen, dtype='complex')
conyy = zeros(vecLen, dtype='complex')
conzz = zeros(vecLen, dtype='complex')
conxy = zeros(vecLen, dtype='complex')
conxz = zeros(vecLen, dtype='complex')
conyz = zeros(vecLen, dtype='complex')


#solve first equation:
x1 = solve_eq1()

conyy = x1[0:vecLen]
conzz = x1[vecLen:2*vecLen]
conyz = x1[2*vecLen:3*vecLen]

#solve second equation:
x2 = zeros(4*vecLen, dtype='complex')
#solve equation 2:
while True:  
    (J_x0, f_x0)  = solve_eq2(x2)
    dx = linalg.solve(J_x0, -f_x0)
    x2 = x2 + dx
    print linalg.norm(f_x0, 2)
    if (linalg.norm(f_x0,2) < NRdelta): break

U = x2[0:vecLen]
conxx=x2[1*vecLen:2*vecLen]
conxy=x2[2*vecLen:3*vecLen]
conxz=x2[3*vecLen:4*vecLen]

save_pickle((U,V,W,conxx,conyy,conzz,conxy,conxz,conyz), filename)
