#!/usr/bin/python

#SLCTM: Simple Light Curve Transit Model with Quadratic Limb Darkening
#Outputs a file with transit times with noise
#**This script takes commandline arguments, input/output in years.**
# >$python SLCTTM.py [outfile] [N of transits] [t0] [c1] [c2] [porb] [pttv] [noisefactor_transit] [b] [vtan] [noisefactor_flux]
import math
import sys
import numpy as np
import batman #Limb Darkening Model from
              #http://astro.uchicago.edu/~kreidberg/batman/quickstart.html


#############################################

#Transit Times with Noise. Units of Days
#To Use SLCTM:
#First Initialize / Create Instance of Object
#Then set Model parameters with SLCTM.SetModelParams and SLCTM.setBatmanParams
#Then Populate Transit times array with SLCTM.populate_transtimes
#Finally, populate flux arrays with SLCTM.populate_fluxes.
#Access Transit times array as SLCTM.transtimes[], fluxes as SLCTM.flux[],
#normalized fluxes as SLCTM.normflux[], and flux times as SLCTM.time[].

class SLCTM:
    def __init__(self):

        #Initialize Simple Model Parameters:
        self.t0 = 0.0        #Time zero
        self.c1 = 0.0        #Fit parameter
        self.c2 = 0.0        #Fit parameter
        self.porb = 0.0      #Period of Orbit
        self.pttv = 0.0      #TTV Period
        self.noisett_e = 0.0 #Noise factor for transit times
        self.b = 0.0         #'impact parameter' - closest to center of star, projected
        self.vtan = 0.0      #tangential velocity of planet
        self.noisef_e = 0.0  #flux noise factor

        #Initialize Batman Parameters:
        self.bmparams = batman.TransitParams()
        self.bmparms.t0 = 0.0                #time of inferior conjunction
        self.bmparms.per = 0.0               #orbital period
        self.bmparms.rp = 0.0                #planet radius (stellar radii)
        self.bmparms.a = 0.0                 #semi major axis (stellar radii)
        self.bmparms.inc = 0.0               #orbital inclination (degrees)
        self.bmparms.ecc = 0.0               #eccentricity
        self.bmparms.w = 0.0                 #longitude of periastron (degrees)
        self.bmparms.u = [0.0,0.0]           #limb darkening coeffs
        self.bmparms.limb_dark = "quadratic" #Use Quadratic Limb Darkening

        #Data Arrays:
        self.transtimes = []        #Transit Times
        self.transitStartIndex = [] #Start array index of ith transit
        self.transitEndIndex = []   #End array index of ith transit
        self.time = []              #Timestamp of Flux Datapoint
        self.flux = []              #Flux Datapoint
        self.normflux = []          #Normalized Flux Datapoint

    def SetModelParams(self,t0,c1,c2,porb,pttv,noisett_e,b,vtan,noisef_e):
        self.t0 = t0
        self.c1 = c1
        self.c2 = c2
        self.porb= porb
        self.pttv = pttv
        self.noisett_e = noisett_e
        self.b = b
        self.vtan = vtan
        self.noisef_e = noisef_e

    def setBatmanParams(self,t0,per,rp,a,inc,ecc,w,u):
        self.bmparms.t0 = t0
        self.bmparms.per = per
        self.bmparms.rp = rp
        self.bmparms.a = a
        self.bmparms.inc = inc
        self.bmparms.ecc = ecc
        self.bmparms.w = w
        self.bmparms.u = u

    def populate_transtimes(self,n):
        #Populates transit times array with n transits
        #according to t_i = t_0 + n*Porb + c1*sin((t_0+n*Porb)/(2pi*P_ttv)) + c2*sin((t_0+n*Porb)/(2pi*P_ttv))
        for i in range(0,n):
            self.transtimes.append(np.random.normal(0,1)*self.noisett_e+(self.t0 + i*self.porb + self.c1*math.sin((self.t0 + i*self.porb)/(2*math.pi*self.pttv))+ self.c2*math.cos((self.t0 + i*self.porb)/(2*math.pi*self.pttv))))
    #def populate_fluxes_naive(self,bmparams):
        ##To Do##
        #Populates fluxes around each transit.
        #Generate times to calculate model around: (self.time[]):

        #Run Batman around each elem of transtimes with model params bmparams



    #def normalize_fluxes(self):
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

#Removed CLI argument parsing for simpler testing. Trying bogus values:

OutputSLCTM = SLCTM()
OutputSLCTM.SetModelParams(0.0,2.2,2.3,20.0,0.01,0.2)
OutputSLCTM.setBatmanParams(0.0,20.0,0.3,120,87,0.0,90.0,[0.1,0.3])
OutputSLCTM.populate_transtimes(int(20))
OutputSLCTM.writetranstimes("testfile")
