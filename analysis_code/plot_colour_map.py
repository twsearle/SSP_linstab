#------------------------------------------------------------------------------
#   colour map plotter for 2D coherent state finder
#
#   Last modified: Fri 21 Feb 16:47:29 2014
#
#------------------------------------------------------------------------------
#TODO check that the axes are the right way up?

#MODULES
import sys
from scipy import *
from scipy import linalg
from scipy import fftpack
import cPickle as pickle
import ConfigParser
from matplotlib import pyplot as plt
from matplotlib import rc

#SETTINGS----------------------------------------

config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
N = config.getint('settings', 'N')
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
Weiss = config.getfloat('settings', 'Weiss')
beta = config.getfloat('settings', 'beta')
Amp = config.getfloat('settings', 'amp')

fp.close()

numYs = 100
numZs = 100
baseFileName = '-N'+str(N)+'-M'+str(M)+'-Re'+str(Re)+'-b'+str(beta)\
          +'-Wi'+str(Weiss)+'-amp'+str(Amp)+'.pickle'

inFileName = 'real-pf-N'+str(N)+'-M'+str(M)+'-Re'+str(Re)+'-b'+str(beta)\
          +'-Wi'+str(Weiss)+'-amp'+str(Amp)+'.pickle'

#------------------------------------------------

# FUNCTIONS

def z_point(zIndex):
    return zLength*((zIndex/(1.*(numZs-1.))) - 0.5)

def y_point(yIndex):
    return -1. + (2./(numYs-1.))*yIndex

# MAIN

# Read in

(U,V,W,Cxx,Cyy,Czz,Cxy,Cxz,Cyz) = pickle.load(open(inFileName, 'r'))
gamma = pi / 2.
zLength = 2.*pi/gamma

z_points = zeros(numZs,dtype='d') 
for zIndx in range(numZs):
    z_points[zIndx] = z_point(zIndx)
del zIndx

y_points = zeros(numYs,dtype='d')
for yIndx in range(numYs):
    y_points[yIndx] = y_point(yIndx)
del yIndx

# make meshes
grid_z, grid_y = meshgrid(z_points, y_points)

#make plots prettier:
#inches_per_Lx = 1.4
#inches_per_Ly = 2.2
#fig_width = 20
#fig_height = 4      
#fig_size =  [fig_width,fig_height]
#rc('figure', figsize=fig_size)

plt.figure()

keyDict = {'Re':Re, 'Wi':Weiss, 'beta':beta, 'amp':Amp}

# Velocity field

extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(U, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
plt.quiver(grid_z, grid_y, W, V, scale=2)
titleString = 'U for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'vel_map{0}.pdf'.format(baseFileName[:-7]))


# Stresses
plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cxx, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'Cxx for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'Cxx_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cyy, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'Cyy for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'Cyy_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cxx-Cyy, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'N1 Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'N1_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cyy-Czz, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'N2 Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'N2_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Czz, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'Cxx for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'Czz_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cxy, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'Cxy for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'Cxy_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cxz, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'Cxz for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'Cxz_map{0}.pdf'.format(baseFileName[:-7]))

plt.figure()
extent_=[-0.5*zLength, 0.5*zLength,-1,1]
plt.imshow(Cyz, origin='lower', extent=extent_) 
plt.colorbar(orientation='horizontal')
titleString = 'Cyz for Re = {Re}, beta = {beta}, Wi = {Wi}, amp = {amp}'.format(**keyDict)
plt.title(titleString)
plt.savefig(r'Cyz_map{0}.pdf'.format(baseFileName[:-7]))

