#-------------------------------------------------------------------------------
#   Plot of the streamfunction of the base profile for test
#
#
#-------------------------------------------------------------------------------
"""
Plots only the final state of the time iterated system.
"""


#MODULES
from scipy import *
from scipy import linalg
import cPickle as pickle
import sys
from matplotlib import pyplot as plt
from matplotlib import rc
import ConfigParser
import argparse

#SETTINGS---------------------------------------------------

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

fp.close()

baseFileName = '-k{k}-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    k=k, N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)


inFileName = 'full-evecs-k{k}-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    k=k, N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

# SETUP ARGUMENT PARSER

parser = argparse.ArgumentParser(description='plot Fourier modes from data file')
parser.add_argument('-f','--filename', default=inFileName,
                     help='input filename')
args = parser.parse_args()
inFileName = args.filename
print inFileName

#-----------------------------------------------------------

def Cheb_to_real_transform(vec, y_points) :
    """ calculate the Fourier chebychev transform for the 2D coherent state
    finder"""
    rVec = zeros(numYs, dtype='complex')
    for yIndx in range(numYs):
        y = y_points[yIndx]
        for m in range(M):
            term = vec[m] * cos(m*arccos(y))
            rVec[yIndx] += term
    del y,m

    return rVec

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

#
#   MAIN
#

vecLen = (2*N+1)*M
# Set the oneOverC function: 1/2 for m=0, 1 elsewhere:
oneOverC = ones(M)
oneOverC[0] = 1. / 2.
# Set up the CFunc function: 2 for m=0, 1 elsewhere:
CFunc = ones(M)
CFunc[0] = 2.
Psi = pickle.load(open(inFileName, 'r'))

numYs = 50

y_points = zeros(numYs, dtype='d')
for yIndx in range(numYs):
    y_points[yIndx] = 1.0 - (2.0*yIndx)/(numYs-1.0)
del yIndx

(du,dv,dw,_,dcxx,dcyy,dczz,dcxy,dcxz,dcyz) = pickle.load(open(inFileName, 'r'))


#make plots prettier:
inches_per_Lx = 1.4
inches_per_Ly = 2.2
fig_width = 10 
fig_height = 2*2*inches_per_Ly      
fig_size =  [fig_width,fig_height]
rc('figure', figsize=fig_size)

for n in range(N,2*N+1):
    u = Cheb_to_real_transform(du[n*M: (n+1)*M], y_points)
    v = Cheb_to_real_transform(dv[n*M: (n+1)*M], y_points)
    w = Cheb_to_real_transform(dw[n*M: (n+1)*M], y_points)
    plt.figure()
    ax1 = plt.subplot(311)
    titleString = 'u  n = {mode} mode'.format(mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(u), 'b-')
    plt.plot(y_points, imag(u), 'r-')
    ax1.axhline(linewidth=0.5, linestyle='--', color='k')
    ax1.axvline(linewidth=0.5, linestyle='--', color='k')
    ax2 = plt.subplot(312)
    titleString = 'v  n = {mode} mode'.format(mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(v), 'b-')
    plt.plot(y_points, imag(v), 'r-')
    ax2.axhline(linewidth=0.5, linestyle='--', color='k')
    ax2.axvline(linewidth=0.5, linestyle='--', color='k')
    ax3 = plt.subplot(313)
    titleString = 'w  n = {mode} mode'.format(mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(w), 'b-')
    plt.plot(y_points, imag(w), 'r-')
    ax3.axhline(linewidth=0.5, linestyle='--', color='k')
    ax3.axvline(linewidth=0.5, linestyle='--', color='k')
    plt.savefig('pertb-n{n}{pf}.pdf'.format(n=n-N, pf=baseFileName[:-7]))
    plt.show()

    
#savetxt('test.dat', vstack((real(PSIr1), imag(PSIr1))).T)
