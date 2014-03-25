#----------------------------------------------------------------------------#
#   Check for Positive definiteness of                              24-01-13
#   Conformation tensors                              positive_definite_R.py
#   Last modified: Tue 25 Mar 17:16:47 2014
#----------------------------------------------------------------------------#


# MODULES
import sys
import ConfigParser
from scipy import *
from scipy import linalg
from matplotlib import pyplot as plt
import cPickle as pickle

# FUNCTIONS

def y_point(yIndex):
    return -1. + (2./(yDataPts-1.))*yIndex

def z_point(zIndex):
    return zLength*((zIndex/(1.*(zDataPts-1.))) - 0.5)

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

# MAIN
#
#
config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
N = config.getint('settings', 'N')
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')

fp.close()

gamma = pi / 2.
zLength = 2.*pi/gamma

Weissenbergs = r_[2.:20.:2.]

# Loop over all files
for Weiss in Weissenbergs:

    print Weiss
    # Read file data into program
    basefilename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.pickle'.format(
                    N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

    (U,_,_, Conxx,
                Conyy,
                Conzz,
                Conxy,
                Conxz,
                Conyz) = pickle.load(open("real-pf" + basefilename, 'r'))
        
    yDataPts = len(U[0,:])
    zDataPts = len(U[:,0])
    #convert to 1D:
    Conxx = Conxx.flatten()
    Conyy = Conyy.flatten()
    Conzz = Conzz.flatten()
    Conxy = Conxy.flatten()
    Conxz = Conxz.flatten()
    Conyz = Conyz.flatten()

    # Calculate eigenvalues and save to array
    eigenvalues = zeros((zDataPts*yDataPts, 3))

    for i in range(zDataPts*yDataPts):

        pointC = array([[Conxx[i], Conxy[i], Conxz[i]], \
                        [Conxy[i], Conyy[i], Conyz[i]], \
                        [Conxz[i], Conyz[i], Conzz[i]]])
        eigenvalues[i] = linalg.eigvalsh(pointC)

    del i
        
    # check eigenvalues are positive for these settings before looping again.
    booleigs = greater(eigenvalues, zeros((zDataPts*yDataPts, 3)))
    TrueArray = ones((yDataPts*zDataPts, 3))

    #save bool array to file for colour map
    boolmap = zeros(yDataPts*zDataPts)
    for i in range(yDataPts*zDataPts):
        boolmap[i] = array_equiv(booleigs[i], ones(3))
    del i
    
    #make boolmap 2D:
    boolmap = reshape(boolmap, (zDataPts, yDataPts)).T
    
    #Print positive definiteness test as contour map data:
    # save_field(boolmap, 'posdef' + basefilename[:-7])

    plt.figure()
    extent_=[-0.5*zLength, 0.5*zLength,-1,1]
    plt.imshow(boolmap, origin='lower', extent=extent_, cmap='Greys_r', vmin=0,
               vmax=1.0) 
    plt.colorbar(orientation='horizontal')
    keyDict = {'Re':Re, 'Wi':Weiss, 'beta':beta, 'amp':Amp}
    titleString =\
    'positive definite map for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
    plt.title(titleString)
    plt.savefig(r'posdef_map{0}.pdf'.format(basefilename[:-7]))

    # Output result to terminal
    print "For Weiss: ", Weiss
    print "All eigenvalues of Conformation tensor are positive is:"
    print allclose(booleigs,TrueArray)

del Weiss
