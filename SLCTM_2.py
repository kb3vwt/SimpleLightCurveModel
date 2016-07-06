#!/usr/bin/python

#SLCTM: Simple Light Curve Transit Model with Quadratic Limb Darkening
#Outputs a file with transit times with noise and a file with times and fluxes
import math
import sys
import numpy as np
import batman #Limb Darkening Model from
              #http://astro.uchicago.edu/~kreidberg/batman/quickstart.html

#For Testing:
import matplotlib.pyplot as plt

#######################
#Things to Do:
#o Implement Error on fluxes
#o Multiplanet (vectorize Porb,...)
#o Check Validity of Noise on Light Curve
#o Named Parameters (take vector/dictionary)
#o Possibly adaptively assign # of datapoints to a transit (save RAM)
#  s.t. time between points << transit duration
#o Remove member (test) functions, create 'regular' functions
#  that take SLCTM as input
#######################

#############################################
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
        #self.noisef_e = 0.0  #flux noise factor

        #Data Arrays:
        self.transtimes = []          #Transit Times
        self.fluxtimes = []           #Timestamp of Flux Datapoint
        self.noiselessfluxes = []     #Retainer for Fluxes before noise addition
        self.fluxes = []              #Flux Datapoint
        self.normfluxes = []          #Normalized Flux Datapoint
        self.transitStartIndex = []   #Start array index of ith transit
        self.transitEndIndex = []     #End array index of ith transit
        self.sigma = []               #Uncertainty (for flux) at index

        #ChiSquared Test
        self.ChiSqComp = []           #Components of ChiSq, to be summed
        self.ChiSq = 0.0

        #Eventual Multiplanet Support Requires:
        #self.transitWhichPlanet = [] #Which Planet is transiting
        #self.porb = []
        #self.pttv = []
        #self.c1 = []
        #self.c2 = []
        #self.b = []
        #self. vtan = []

    def SetModelParams(self,t0,c1,c2,porb,pttv,noisett_e,b,vtan):
        self.t0 = t0
        self.c1 = c1
        self.c2 = c2
        self.porb= porb
        self.pttv = pttv
        self.noisett_e = noisett_e
        self.b = b
        self.vtan = vtan
        #self.noisef_e = noisef_e

    def setBatmanParams(self,BatmanParams,t0,per,rp,a,inc,ecc,w,u1,u2):
        BatmanParams.t0 = t0
        BatmanParams.per = per
        BatmanParams.rp = rp
        BatmanParams.a = a
        BatmanParams.inc = inc
        BatmanParams.ecc = ecc
        BatmanParams.w = w
        BatmanParams.u = [u1, u2]
        #Set to Quadratic Limb Darkening - Possibly open to change?
        BatmanParams.limb_dark = "quadratic"

    def test_PopTransTimes(self,n):
        #Populates transit times array with n transits
        #according to t_i = t_0 + n*Porb + c1*sin((t_0+n*Porb)/(2pi*P_ttv)) + c2*sin((t_0+n*Porb)/(2pi*P_ttv))
        self.__NTransits = n
        for i in range(0,n):
            self.transtimes.append(np.random.normal(0,1)*self.noisett_e+(self.t0 + i*self.porb + self.c1*math.sin((self.t0 + i*self.porb)/(2*math.pi*self.pttv))+ self.c2*math.cos((self.t0 + i*self.porb)/(2*math.pi*self.pttv))))
    def test_PopFluxesNaive(self,BatmanParams,NFluxPoints):
        #Populates fluxes for each transit.
        #Naive - uses length of transit in plane ( inc = 0 deg )

        #Flux Times for This transit:
        #Calculate Interval duration:
        self.__HalfTransitAngle = 2.0 * math.asin((BatmanParams.rp + 1) / (2.0 * BatmanParams.a))
        self.__HalfTransitTime = 1.5*(self.__HalfTransitAngle / (2.0 * math.pi)) * self.porb


        #For ith transit, calculate fluxes using Batman
        for i in range(0,self.__NTransits):
            #Write Start Index:
            if len(self.fluxtimes) == 0:
                self.transitStartIndex.append(0)
            else:
                self.transitStartIndex.append(len(self.fluxtimes) - 1)

            #Set t0 for this transit:
            BatmanParams.t0 = self.transtimes[i]
            self.__BMFluxTimes = np.linspace(self.transtimes[i] - self.__HalfTransitTime, self.transtimes[i] + self.__HalfTransitTime, NFluxPoints) #Temporary array of times for computation
            #append times to master flux times:
            self.fluxtimes.extend(self.__BMFluxTimes)

            #Initialize Batman model, compute fluxes for above interval:
            bmmodel = batman.TransitModel(BatmanParams,self.__BMFluxTimes)
            self.__BMFluxes = bmmodel.light_curve(BatmanParams) #Temporary flux array
            #Append new fluxes to master flux array:
            self.fluxes.extend(self.__BMFluxes)

            #Write End Index:
            if len(self.fluxtimes) == 0:
                self.transitEndIndex.append(0)
            else:
                self.transitEndIndex.append(len(self.fluxtimes) - 1)

    def add_norm_lcnoise(self,lcnoise):
        #Retain noiseless fluxes:
        for i in range(0,len(self.fluxes)):
            self.noiselessfluxes.append(self.fluxes[i])

        #Apply Noise to Fluxes
        for i in range(0,len(self.fluxes)):
            self.fluxes[i] = self.fluxes[i] + np.random.normal(0,1)*lcnoise
        #Renormalize after noise:
        self.normfluxes = np.divide(self.fluxes,np.amax(self.fluxes))
        self.fluxes = self.normfluxes

        #Compute ChiSq:
        #for i in range(0,len(self.fluxes)):
        #    self.ChiSqComp.append((self.fluxes[i] - self.noiselessfluxes[i])**2 / (self.noiselessfluxes[i]))

        self.ChiSqComp = (self.fluxes - self.noiselessfluxes)**2 / (self.noiselessfluxes)

        self.ChiSq = np.sum(self.ChiSqComp)

    def test_WriteTransTimes(self,ofname):
        ofile = open(ofname+'_ttimes.txt','wb')
        for i in range(0,len(self.transtimes)):
            ofile.write(str(self.transtimes[i]) + '\n')
        ofile.close()
    def test_WriteColsTimesFluxes(self,ofname):
        ofile = open(ofname+'_fluxes.txt','wb')
        for i in range(0,len(self.fluxtimes)):
            ofile.write(str(self.fluxtimes[i]) + ',' + str(self.fluxes[i]) + '\n')
        ofile.close()





##########################################################################################
#####TESTING##############################################################################

#Testing Single Light Curve:
#Removed CLI argument parsing for simpler testing. Trying bogus values:

#OutputSLCTM = SLCTM()
#bmparams = batman.TransitParams()
#OutputSLCTM.SetModelParams(0.0,0.05,0.08,10.0,100,0.0001,100.0,200.0)
#OutputSLCTM.setBatmanParams(bmparams,0.0,10.0,0.2,12,90,0.0,90.0,0.1,0.3)
#OutputSLCTM.test_PopTransTimes(int(10))
#OutputSLCTM.test_PopFluxesNaive(bmparams,1000)
#OutputSLCTM.add_norm_lcnoise(0.0005)
#OutputSLCTM.test_WriteTransTimes("testfile")
#OutputSLCTM.test_WriteColsTimesFluxes("testfile")

#plt.plot(OutputSLCTM.fluxtimes,OutputSLCTM.fluxes)
#plt.show()


#Testing Array of Light Curves with various TTV periods. Plots ChiSq(PTTV) vs PTTVs
#TTV Periods:
PTTV_arr = np.linspace(80,120,50)
#ChiSq Array:
ChiSqs = []
#Batman Params:
bmparams = batman.TransitParams()
#List of LCs:
LightCurves = [SLCTM() for i in range(len(PTTV_arr))]
for i in range(len(PTTV_arr)):
    LightCurves[i].SetModelParams(0.0,0.05,0.08,10.0,PTTV_arr[i],0.0001,100.0,200.0)
    LightCurves[i].setBatmanParams(bmparams,0.0,10.0,0.2,12,90,0.0,90.0,0.1,0.3)
    LightCurves[i].test_PopTransTimes(int(1))
    LightCurves[i].test_PopFluxesNaive(bmparams,1000)
    LightCurves[i].add_norm_lcnoise(0.0005)
    ChiSqs.append(LightCurves[i].ChiSq)

#Plots ChiSq(PTTV) vs PTTVs
plt.plot(PTTV_arr,ChiSqs)
plt.title("$\chi^2$ vs $P_{ttv}$")
plt.xlabel("$P_{ttv}$ [Days]")
plt.ylabel("$\chi^2$")
plt.show()

#Plots ChiSq Components of a single model against time
#plt.plot(LightCurves[0].fluxtimes,LightCurves[0].ChiSqComp)
#plt.title("$\chi^2$ Components vs $P_{ttv}$")
#plt.xlabel("Time [Days]")
#plt.ylabel("$\chi^2$ Component")
#plt.show()
