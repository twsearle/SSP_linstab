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
parser.add_argument('-t','--title', default=None,
                     help='description to put in the title')
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
#Psi = pickle.load(open(inFileName, 'r'))

numYs = 50

y_points = zeros(numYs, dtype='d')
for yIndx in range(numYs):
    y_points[yIndx] = 1.0 - (2.0*yIndx)/(numYs-1.0)
del yIndx

(du,dv,dw,_,dcxx,dcyy,dczz,dcxy,dcxz,dcyz) = pickle.load(open(inFileName1, 'r'))
(du2,dv2,dw2,_,dcxx2,dcyy2,dczz2,dcxy2,dcxz2,dcyz2) = pickle.load(open(inFileName2, 'r'))

options1 = {'du': du, 'dv': dv, 'dw':dw, 'dcxx': dcxx, 'dcyy': dcyy, 
            'dczz':dczz, 'dcxy':dcxy, 'dcxz':dcxz, 'dcyz': dcyz}
options2 = {'du': du2, 'dv': dv2, 'dw':dw2, 'dcxx': dcxx2, 'dcyy': dcyy2, 'dczz':
            dczz2, 'dcxy': dcxy2, 'dcxz':dcxz2, 'dcyz': dcyz2}

#make plots prettier:
inches_per_Lx = 1.4
inches_per_Ly = 2.2
fig_width = 10 
fig_height = 2*2*inches_per_Ly      
fig_size =  [fig_width,fig_height]
rc('figure', figsize=fig_size)

# Set scaling based on first v mode 
nscaling = N+1
var = Cheb_to_real_transform(dv[nscaling*M: (nscaling+1)*M], y_points)
var2 = Cheb_to_real_transform(dv2[nscaling*M: (nscaling+1)*M], y_points)

if max(abs(var)) > max(abs(var2)):
    print nscaling-N, ": var bigger than var2 "
    varscale = (max(abs(var))/max(abs(var2)))
else:
    print nscaling-N, ": var smaller than var2 "
    varscale = (max(abs(var2))/max(abs(var)))

for comp in options1: 
    plotNum=311
    plt.figure()
    for n in range(N,2*N+1):

        var = options1[comp]
        var2 = options2[comp]

        var = Cheb_to_real_transform(var[n*M: (n+1)*M], y_points)
        var2 = Cheb_to_real_transform(var2[n*M: (n+1)*M], y_points)

        # apply scaling 
        var2 = varscale*var2
        ax1 = plt.subplot(plotNum)
        titleString = args.title+'{comp}  n = {mode} mode'.format(comp=comp, mode=n-N)
        plt.title(titleString)
        plt.plot(y_points, real(var), 'b-')
        plt.plot(y_points, imag(var), 'r-')
        plt.plot(y_points, real(var2), 'c--')
        plt.plot(y_points, imag(var2), 'm--')
        ax1.axhline(linewidth=.5, linestyle='-', color='k')
        ax1.axvline(linewidth=.5, linestyle='-', color='k')
        #ax1.legend(["real file1", "imag file1", "real file1", "imag file2"])

        plotNum+=1
        if n-N > 1: break

    plt.savefig('pertb_all-{comp}-{pf}.pdf'.format(comp=comp, pf='comparison'))
    plt.show()

    plt.figure()
    plotNum=311

for n in range(N,2*N+1):

    var = dcxx-dcyy
    var2 = dcxx2-dcyy2

    var = Cheb_to_real_transform(var[n*M: (n+1)*M], y_points)
    var2 = Cheb_to_real_transform(var2[n*M: (n+1)*M], y_points)

    ax1 = plt.subplot(plotNum)
    titleString = args.title + '{comp}  n = {mode} mode'.format(comp='N1', mode=n-N)
    plt.title(titleString)
    plt.plot(y_points, real(var), 'b-')
    plt.plot(y_points, imag(var), 'r-')
    plt.plot(y_points, real(var2), 'c--')
    plt.plot(y_points, imag(var2), 'm--')
    ax1.axhline(linewidth=.5, linestyle='-', color='k')
    ax1.axvline(linewidth=.5, linestyle='-', color='k')
    #ax1.legend(["real file1", "imag file1", "real file1", "imag file2"])

    plotNum+=1
    if n-N > 1: break

plt.savefig('pertb-{comp}-{pf}.pdf'.format(comp='N1', pf='comparison'))
plt.show()
    
#savetxt('test.dat', vstack((real(PSIr1), imag(PSIr1))).T)
