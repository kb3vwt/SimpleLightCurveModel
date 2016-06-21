#!/usr/bin/python

#SLCTM: Simple Light Curve Transit Model
#Outputs a file with transit times with noise
#**This script takes commandline arguments, input/output in years.**
# >$python SLCTTM.py [outfile] [N of transits] [t0] [c1] [c2] [porb] [pttv] [noisefactor_transit] [b] [vtan] [noisefactor_flux]
import math
import sys
import numpy as np


#############################################
#Classes and Functions

#Transit Times with Noise. Units of Days

class SLCTM:
    def __init__(self,t0,c1,c2,porb,pttv,noisett_e,b,vtan,noisef_e):
        self.t0 = t0 #Time zero
        self.c1 = c1 #Fit parameter
        self.c2 = c2 #Fit parameter
        self.porb= porb #Period of Orbit
        self.pttv = pttv #TTV Period
        self.noisett_e = noisett_e #Noise factor for transit times
        self.b = b #'impact parameter' - closest to center of star, projected
        self.vtan = vtan #tangential velocity of planet
        self.noisef_e = noisef_e #flux noise factor

        #Transit Times:
        self.transtimes = []

        #Fluxes with times
        self.time = []
        self.flux = []
        self.normfluxes = []
    def populate_transtimes(self,n):
        #Populates transit times array with n transits
        #according to t_i = t_0 + n*Porb + c1*sin((t_0+n*Porb)/(2pi*P_ttv)) + c2*sin((t_0+n*Porb)/(2pi*P_ttv))
        for i in range(0,n):
            self.transtimes.append(np.random.normal(0,1)*self.noisett_e+(self.t0 + i*self.porb + self.c1*math.sin((self.t0 + i*self.porb)/(2*math.pi*self.pttv))+ self.c2*math.cos((self.t0 + i*self.porb)/(2*math.pi*self.pttv))))
    def populate_fluxes(self):
        ##To Do##
        #Populates fluxes around each transit.

    def normalize_fluxes(self):
        ##To Do##
    def writetranstimes(self,ofname):
        ofile = open(ofname+'_ttimes.txt','wb')
        for i in range(0,len(self.transtimes)):
            ofile.write(str(self.transtimes[i]) + '\n')
        ofile.close()

    def writetranstimes(self,ofname):
        ofile = open(ofname+'_normfluxes.txt','wb')
        for i in range(0,len(self.normfluxes)):
            ofile.write(str(self.normfluxes[i]) + '\n')
        ofile.close()





##########################################################################################
#####TESTING##############################################################################

#Take Arguments:
try:
    arg_ofilename = str(sys.argv[1])
except:
    arg_ofilename = ''
try:
    arg_numtrans = int(sys.argv[2])
except:
    arg_numtrans = 0
try:
    arg_t0 = float(sys.argv[3])
except:
    arg_t0 = 0
try:
    arg_c1 = float(sys.argv[4])
except:
    arg_c1 = 0
try:
    arg_c2 = float(sys.argv[5])
except:
    arg_c2 = 0
try:
    arg_porb = float(sys.argv[6])
except:
    arg_porb = 0
try:
    arg_pttv = float(sys.argv[7])
except:
    arg_pttv = 0
try:
    arg_noisett_e = float(sys.argv[8])
except:
    arg_noisett_e = 0.0
try:
    arg_b = float(sys.argv[9])
except:
    arg_b = 0.0
try:
    arg_vtan = float(sys.argv[10])
except:
    arg_vtan = 0.0
try:
    arg_noisef_e = float(sys.argv[11])
except:
    arg_noisef_e = 0.0



#Test Output
if (arg_ofilename == 'h' or arg_ofilename == ''):
    print('SLCTM: Simple Light Curve Transit Model:')
    print('Outputs files with transit times and fluxes.')
    print('-->This script takes these parameters: ')
    print('   $SLCTTM.py [outfile (omit .txt)] [N of transits] [t0] [c1] [c2] [porb] [pttv] [noisefactor_t] [b] [vtan] [noisefactor_flux]')
    print('-->Input arguments in units of YEARS, METERS, METERS/SEC. Outputs file with units of Normalized Flux.')
    print('\n')
else:
    print('oops')


    OutputSLCTM = SLCTTM(arg_t0,arg_c1,arg_c2,arg_porb,arg_pttv,arg_noisett_e)
    OutputSLCTM.populate_transtimes(int(arg_numtrans))
#    OutputSLCTM.writetranstimes(arg_ofilename)
