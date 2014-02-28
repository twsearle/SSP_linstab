# ----------------------------------------------------------------------------
#
#   maximum eigenvalues finder
#   Last modified: Fri 28 Feb 19:36:56 2014
#
# ----------------------------------------------------------------------------
"""
Program to find maximum eigenvalues for various kx at the settings provided in
the settings file. The method is kind of ugly, but ought to produce eigenvalues
ordered by kx and then save them to a data file with the lead- prefix.
upperEigLimit sets the maximum eigenvalue I am prepared to consider as less than
infinity.
"""

from scipy import *
import sys
import ConfigParser
import argparse
import glob
import operator

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
argparser.add_argument("-kl","--kx-list", nargs= '+', type=float, 
                help='run the program in stupid mode, where you have to'+
                       " specify a string containing a csv list of kx's")

args = argparser.parse_args()
N = args.N 
M = args.M
Re = args.Re
beta = args.b
Weiss = args.Wi
Amp = args.amp
ksettings = args.kx_list

baseFilename = '-N{N}-M{M}-Re{Re}-b{b}-Wi{Wi}-amp{amp}.dat'.format(
               N=N, M=M, Re=Re, b=beta, Wi=Weiss, amp=Amp)

upperEigLimit = 100

# -----------------------------------------------------------------------------

# MAIN

if args.kx_list:
    print 'Running in stupid mode'

    leading_eigs = zeros((len(ksettings),3))
    leading_eigs[:,0] = ksettings

    for kindex, k in enumerate(ksettings):
        filename = 'ev-k{k}{base}'.format(k=k, base=baseFilename)
        eigs = genfromtxt(filename)
        eig = [2*upperEigLimit,0]
        while(eig[0] > upperEigLimit): 
            row_index = argmax(eigs[:,0])
            eig = eigs[row_index,:]
            leading_eigs[kindex, 1:3] = eig
            if (eig[0] > upperEigLimit):
                eigs[row_index,:] = [0,0]
                eig = [2*upperEigLimit,0] 
    
    savetxt('lead-ev{0}'.format(baseFilename), leading_eigs, fmt='%10.4f')

else:
    print 'Finding eigenvalues' 
    fileList = glob.glob('*.dat')
    
    fnameKxArr = []
    for filename in fileList:
        splitString = filename.split('-')

        if len(splitString) < 8: continue
        if splitString[0] != 'ev': continue
        sk = splitString[1]
        sN = splitString[2]
        sM = splitString[3]
        sRe = splitString[4]
        sb = splitString[5]
        sWi = splitString[6]
        sAmp = splitString[7]
       
        if sk[0] != 'k': continue
        if sN[0] != 'N': continue
        if sM[0] != 'M': continue
        if sRe[:2] != 'Re': continue
        if sb[0] != 'b': continue
        if sWi[:2] != 'Wi': continue
        if sAmp[:3] != 'amp': continue

        fk = float(sk[1:]) 
        fN = int(sN[1:])
        fM = int(sM[1:])
        fRe = float(sRe[2:])
        fb = float(sb[1:])
        fWi = float(sWi[2:])
        fAmp = float(sAmp[3:-4])

        kwargs = {'fk':fk, 'fN':fN, 'fM':fM, 'fRe':fRe, 'fb':fb, 'fWi':fWi, 'fAmp':fAmp}
        #print 'kx = {fk}, N={fN}, M={fM}, Re={fRe}, b={fb}, Wi={fWi}, amp={fAmp}'.format(**kwargs)

        if fN != N: continue
        if fM != M: continue
        if fRe != Re: continue
        if fRe != Re: continue
        if fWi != Weiss: continue
        if fAmp != Amp: continue


        fnameKxArr.append([filename, fk])
    del filename 

    fnameKxArr.sort(key=operator.itemgetter(1))
    fnameKxArr = array(fnameKxArr) 

    fileList = fnameKxArr[:,0]
    kxList = fnameKxArr[:,1]

    leading_eigs = 2*upperEigLimit*ones((len(fileList), 3), dtype='d')
    for i, filename in enumerate(fileList):
        dispersionData = genfromtxt(filename)

        while leading_eigs[i,1] > upperEigLimit:
            fk = float(kxList[i])
            maxIndx = argmax(dispersionData[:,0])
            leading_eigs[i,:] = [fk, dispersionData[maxIndx,0],
                                 dispersionData[maxIndx,1]]
            if leading_eigs[i,1] > upperEigLimit:
                leading_eigs[i,1] = 2*upperEigLimit
                dispersionData[maxIndx, :] = 0 

        print 'kx = ', fk

    del i, filename

    savetxt('lead-ev{0}'.format(baseFilename), leading_eigs, fmt='%10.4f')
