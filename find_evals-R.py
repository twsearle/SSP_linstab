#Last modified: Fri 21 Feb 17:15:34 2014
from scipy import *
import sys
import ConfigParser

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

setting_str = sys.argv[1]
ksettings = [float(x) for x in setting_str.split(',')]

leading_eigs = zeros((len(ksettings),3))
leading_eigs[:,0] = ksettings

b_filename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.dat'.format(
               N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

for kindex, k in enumerate(ksettings):
    filename = 'ev-k{k}{base}'.format(k=k, base=b_filename)
    eigs = genfromtxt(filename)
    eig = [100,0]
    while(eig[0] > 50): 
        row_index = argmax(eigs[:,0])
        eig = eigs[row_index,:]
        leading_eigs[kindex, 1:3] = eig
        if (eig[0] > 50):
            eigs[row_index,:] = [0,0]
            eig = [100,0] 
    


savetxt('lead-ev{0}'.format(b_filename), leading_eigs, fmt='%10.4f')
