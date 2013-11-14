# -----------------------------------------------------------------------------
#
#   Graph the Inter modal separation of an eigenvalue spectrum
#   Last modified: Thu 14 Nov 15:36:35 2013
#
# -----------------------------------------------------------------------------
"""
Chapter 7 of Boyd's book Chebyshev and Fourier spectral methods suggests that a
good way of checking eigenvalue convergence is by comparing the ordinal
difference between adjacent eigenvalues. Should be able to use this to filter
out the physical ones.

"""

# MODULES
from scipy import *
from matplotlib import pylab as plt
from matplotlib import rc
import ConfigParser

# PARAMETERS ------------------------------------------------------------------

config = ConfigParser.RawConfigParser()
fp = open('OB-settings.cfg')
config.readfp(fp)

Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
Amp = config.getfloat('settings', 'Amp')
k = config.getfloat('settings', 'k')

fp.close()

N1 = 7
M1 = 40
N2 = 7
M2 = 60

args = {'k': k, 'N': N1, 'M': M1, 'Re': Re, 'b': beta, 'Wi': Weiss, 'amp': Amp}
filename1 = 'ev-k{k}-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.dat'.format(**args)

args = {'k': k, 'N': N2, 'M': M2, 'Re': Re, 'b': beta, 'Wi': Weiss, 'amp': Amp}
filename2 = 'ev-k{k}-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.dat'.format(**args)

# -----------------------------------------------------------------------------

# FUNCTIONS

def calc_sigma(evs):
    """
    Calculates the scaling factor to calculate the ordinal and nearest
    differences.

    """

    sigma = zeros(vecLen-1, dtype='d')
    sigma[0] = absolute(evs[0] - evs[1])

    for j in range(1, vecLen-1):
        term1 = absolute(evs[j] - evs[j-1])
        term2 = absolute(evs[j+1] - evs[j])
        sigma[j] = 0.5 * (term1 + term2)
        if sigma[j] is 0.0:
            print "Degenerate mode at: ", j
    del j
    return sigma

def ordinal_difference(evs, evs2):
    """
    Finds an array of ordinal differences given two spectra. 

    do_j = |lambda_j(N1) - lambda_j(N2)| / sigma_j

    Only works if eigenvalues are evenly spaced. I think this means that it
    only works if you can guarantee that the eigenvalues in the different
    spectra correspond to each other. We can't guarantee this.

    """
    
    sigma = calc_sigma(evs)
 
    # print sigma
    # ordinal difference
    ordDiff = zeros(vecLen-1, dtype='d')
    for j in range(vecLen-1):
        ordDiff[j] = absolute(evs[j] - evs2[j]) / sigma[j]
    del j
    return ordDiff

def nearest_difference(evs1, evs2):
    """
    Find the distance between the two nearest eigenvalues.

    dn_j = min |lambda_j(N1) - lamda_i(N2)| / sigma_j
          -----
          k in [1, N2]

    """

    sigma = calc_sigma(evs1)
    nearestDiff = zeros((vecLen-1), dtype='d')
    for j in range(vecLen-1):
        minimum = infty
        for i in range(vecLen2):
            diff = absolute(evs1[j] - evs2[i]) / sigma[j]
            if diff < minimum:
                minimum = diff
        del i
        nearestDiff[j] = minimum
    del j

    return nearestDiff


# MAIN 


evs1 = genfromtxt(filename1)
evs1 = evs1[:,0] + 1.j*evs1[:,1]
evs2 = genfromtxt(filename2)
evs2 = evs2[:,0] + 1.j*evs2[:,1]

vecLen = len(evs1)
vecLen2 = len(evs2)

ordDiff = ordinal_difference(evs1, evs2)
nearestDiff = nearest_difference(evs1, evs2)

#make plots prettier:
inches_per_Lx = 1.4
inches_per_Ly = 2.2
fig_width =  8
fig_height = 4*inches_per_Ly      
fig_size =  [fig_width,fig_height]
rc('figure', figsize=fig_size)


plt.semilogy(r_[0:vecLen-1], ordDiff, 'rx')
plt.semilogy(r_[0:vecLen-1], nearestDiff, 'b.')
plt.show()
