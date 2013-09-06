#----------------------------------------------------------------------------#
#   Velocity profile convergence analysis                           14-02-13
#   (fully spectral)         
#   Last modified: Mon 04 Mar 2013 10:33:18 GMT
#----------------------------------------------------------------------------#
"""Compare the velocity profiles at different spectral resolutions to find lowest
resolution where convergence occurs. The comparison is performed over the Cxx
component because this is both the largest and the most sensitive to error."""

# MODULES
import sys
from scipy import *
import cPickle as pickle
import ConfigParser

#
# MAIN
#

# Specify profile settings
config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')

Nbegin = config.getint('settings', 'N')
Nend   = config.getint('settings', 'Nend')
Nstep  = config.getint('settings', 'Nstep')
Mbegin = config.getint('settings', 'M')
Mend   = config.getint('settings', 'Mend')
Mstep  = config.getint('settings', 'Mstep')

fp.close()

#Narray = r_[Nbegin:(Nend+Nstep):Nstep]
#Marray = r_[Mbegin:(Mend+Mstep):Mstep]
#settings = vstack((Narray, Marray)).T
settings = array([[5,40], [8,60]])
print 'settings array:'
print settings
 

tolerance = 1.e-5
CompN = 5
CompM = 10
CompArray = zeros( ((2*CompN+1)*CompM, len(settings[:,0])), dtype='complex' )

basefilename = '-Re'+str(Re)+'-b'+str(beta)\
              +'-Wi'+str(Weiss)+'-amp'+str(Amp)+'.pickle'
# Read in profile files

for sitr in range( len(settings[:,0]) ):
    N = settings[sitr,0]
    M = settings[sitr,1]
    filename = 'pf-N{N}-M{M}{basefn}'.format(N=N,M=M,basefn = basefilename)
    (_,_,_,conxx,_,_,_,_,_) = pickle.load(open(filename, 'r'))
    for i in range(CompN):
        CompArray[i*CompM:(i+1)*CompM, sitr] = conxx[(N+i)*M:(N+i)*M + CompM]
    del i
del sitr 

# Perform Comparison over Cxx component
print 'Checking convergence of Cxx with absolute tolerance:', tolerance
for sitr in range(len(settings)-1):
    print settings[sitr,0]
    print 'comparing:\t N = {0} M = {1}'.format(settings[sitr,0],settings[sitr,1])
    print 'with\t\t N = {0} M = {1}'.format(settings[sitr+1,0],settings[sitr+1,1])
    print '\t\t---------------- ',\
            allclose(CompArray[:,sitr], CompArray[:,sitr+1], atol=tolerance)
del sitr

#####################################TEST###################################
# savetxt('compArray.dat', CompArray, fmt='%3.3f')
