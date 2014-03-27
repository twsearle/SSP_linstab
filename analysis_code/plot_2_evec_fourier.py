#-------------------------------------------------------------------------------
#   Plot of the streamfunction of the base profile for test
#
#
#-------------------------------------------------------------------------------
"""
Plot two ECS solutions side by side 
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
parser.add_argument('-f1','--filename1', default=inFileName,
                     help='input filename 1')
parser.add_argument('-f2','--filename2', default=inFileName,
                     help='input filename 2')
args = parser.parse_args()
inFileName1 = args.filename1
inFileName2 = args.filename2
print "Eigenvectors for:"
print inFileName1
print inFileName2
print "WARNING: MUST HAVE THE SAME RESOLUTION"

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

(du,dv,dw,_,_,_,_,_,_,_) = pickle.load(open(inFileName1, 'r'))
(du2,dv2,dw2,_,_,_,_,_,_,_) = pickle.load(open(inFileName2, 'r'))


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
    u2 = Cheb_to_real_transform(du2[n*M: (n+1)*M], y_points)
    v2 = Cheb_to_real_transform(dv2[n*M: (n+1)*M], y_points)
    w2 = Cheb_to_real_transform(dw2[n*M: (n+1)*M], y_points)

    if n == N:
        if max(abs(u)) > max(abs(u2)):
            print n-N, ": u smaller than u2 "
            uscale = (max(abs(u))/max(abs(u2)))
        else:
            print n-N, ": u bigger than u2 "
            uscale = (max(abs(u2))/max(abs(u)))

        if max(abs(v)) > max(abs(v2)):
            print n-N, ": v smaller than v2 "
            vscale = (max(abs(v))/max(abs(v2)))
        else:
            print n-N, ": v bigger than v2 "
            vscale = (max(abs(v2))/max(abs(v)))

        if max(abs(w)) > max(abs(w2)):
            print n-N, ": w smaller than w2 "
            wscale = (max(abs(w))/max(abs(w2)))
        else:
            print n-N, ": w bigger than w2 "
            wscale = (max(abs(w2))/max(abs(w)))
    
        scales = [uscale, vscale, wscale]
        scale = scales[argmin([log(uscale), log(vscale), log(wscale)])]

    u2 = scale*u2
    v2 = scale*v2
    w2 = scale*w2

    plt.figure()
    ax1 = plt.subplot(311)
    titleString = 'u  n = {mode} mode'.format(mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(u), 'b-')
    plt.plot(y_points, imag(u), 'r-')
    plt.plot(y_points, real(u2), 'c--')
    plt.plot(y_points, imag(u2), 'm--')
    ax1.axhline(linewidth=.5, linestyle='-', color='k')
    ax1.axvline(linewidth=.5, linestyle='-', color='k')
    #ax1.legend(["real file1", "imag file1", "real file1", "imag file2"])
    ax2 = plt.subplot(312)
    titleString = 'v  n = {mode} mode'.format(mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(v), 'b-')
    plt.plot(y_points, imag(v), 'r-')
    plt.plot(y_points, real(v2), 'c--')
    plt.plot(y_points, imag(v2), 'm--')
    ax2.axhline(linewidth=.5, linestyle='-', color='k')
    ax2.axvline(linewidth=.5, linestyle='-', color='k')
    ax3 = plt.subplot(313)
    titleString = 'w  n = {mode} mode'.format(mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(w), 'b-')
    plt.plot(y_points, imag(w), 'r-')
    plt.plot(y_points, real(w2), 'c--')
    plt.plot(y_points, imag(w2), 'm--')
    ax3.axhline(linewidth=.5, linestyle='-', color='k')
    ax3.axvline(linewidth=.5, linestyle='-', color='k')
    plt.savefig('pertb-n{n}{pf}.pdf'.format(n=n-N, pf='comparison'))

    #if n-N > 1: break

    plt.show()

    
#savetxt('test.dat', vstack((real(PSIr1), imag(PSIr1))).T)
