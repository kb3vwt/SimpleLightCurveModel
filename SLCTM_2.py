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
#o Possibly adaptively assign # of datapoints to a transit (save RAM)
#  s.t. time between points << transit duration
#######################

#######CLASSES##########################################################################
########################################################################################
class SLCTM:
    def __init__(self):
        #Initialize Simple Model Parameters:
        self.t0 = 0.0                 #Time zero
        self.c1 = 0.0                 #Fit parameter
        self.c2 = 0.0                 #Fit parameter
        self.porb = 0.0               #Period of Orbit
        self.pttv = 0.0               #TTV Period
        self.noisett_e = 0.0          #Noise factor for transit times
        self.b = 0.0                  #'impact parameter' - closest to center of star, projected
        self.vtan = 0.0               #tangential velocity of planet
        #Data Arrays:
        self.transtimes = []          #Transit Times
        self.fluxtimes = []           #Timestamp of Flux Datapoint
        self.fluxes = []              #Fluxes - either a model or Simulated Data
        self.interpfluxtimes = []     #Interpolated Flux Timestamps
        self.interpfluxes = []        #Interpolated Fluxes
        self.sigma = []               #Uncertainty (for flux) at index
        self.isNoisy = False          #Assigned if the model has noise added (simulated data)
        self.transitStartIndex = []   #Start array index of ith transit
        self.transitEndIndex = []     #End array index of ith transit
    def SetModelParams(self,INDIC):
        self.t0 = INDIC['t0']
        self.c1 = INDIC['c1']
        self.c2 = INDIC['c2']
        self.porb= INDIC['porb']
        self.pttv = INDIC['pttv']
        self.noisett_e = INDIC['noisett_e']
        self.b = INDIC['b']
        self.vtan = INDIC['vtan']
    def setBatmanParams(self,BatmanParams,BMIn):
        BatmanParams.t0 = BMIn['t0']
        BatmanParams.per = BMIn['per']
        BatmanParams.rp = BMIn['rp']
        BatmanParams.a = BMIn['a']
        BatmanParams.inc = BMIn['inc']
        BatmanParams.ecc = BMIn['ecc']
        BatmanParams.w = BMIn['w']
        BatmanParams.u = [BMIn['u1'], BMIn['u2']]
        BatmanParams.limb_dark = "quadratic"         #Set to Quadratic Limb Darkening - Possibly open to change?

#####FUNCTIONS##########################################################################
########################################################################################
def PopTransTimes(SLCTM,n):
    #Populates transit times array with n transits
    #according to t_i = t_0 + n*Porb + c1*sin((t_0+n*Porb)/(2pi*P_ttv)) + c2*sin((t_0+n*Porb)/(2pi*P_ttv))
    SLCTM.NTransits = n
    for i in range(0,n):
        SLCTM.transtimes.append(np.random.normal(0,1)*SLCTM.noisett_e+(SLCTM.t0 + i*SLCTM.porb + SLCTM.c1*math.sin((SLCTM.t0 + i*SLCTM.porb)/(2*math.pi*SLCTM.pttv))+ SLCTM.c2*math.cos((SLCTM.t0 + i*SLCTM.porb)/(2*math.pi*SLCTM.pttv))))

def PopFluxesNaive(SLCTM,BatmanParams,NFluxPoints):
    #Populates fluxes for each transit.
    #Naive - uses length of transit in plane ( inc = 0 deg )

    #Flux Times for This transit:
    #Calculate Interval duration:
    timefactor = 2.0 #Use 1.0 if not calculating time outside occlusion
    HalfTransitAngle = 2.0 * math.asin((BatmanParams.rp + 1) / (2.0 * BatmanParams.a))
    HalfTransitTime = timefactor*(HalfTransitAngle / (2.0 * math.pi)) * SLCTM.porb


    #For ith transit, calculate fluxes using Batman
    for i in range(0,SLCTM.NTransits):
        #Write Start Index:
        if len(SLCTM.fluxtimes) == 0:
            SLCTM.transitStartIndex.append(0)
        else:
            SLCTM.transitStartIndex.append(len(SLCTM.fluxtimes) - 1)

        #Set t0 for this transit:
        BatmanParams.t0 = SLCTM.transtimes[i]
        BMFluxTimes = np.linspace(SLCTM.transtimes[i] - HalfTransitTime, SLCTM.transtimes[i] + HalfTransitTime, NFluxPoints) #Temporary array of times for computation
        #append times to master flux times:
        SLCTM.fluxtimes.extend(BMFluxTimes)

        #Initialize Batman model, compute fluxes for above interval:
        bmmodel = batman.TransitModel(BatmanParams,BMFluxTimes)
        BMFluxes = bmmodel.light_curve(BatmanParams) #Temporary flux array
        #Append new fluxes to master flux array:
        SLCTM.fluxes.extend(BMFluxes)

        #Write End Index:
        if len(SLCTM.fluxtimes) == 0:
            SLCTM.transitEndIndex.append(0)
        else:
            SLCTM.transitEndIndex.append(len(SLCTM.fluxtimes) - 1)

def add_norm_lcnoise(SLCTM,lcnoise):
    #Note that this is simulated data:
    SLCTM.isNoisy = True
    #Apply Noise to Fluxes
    for i in range(0,len(SLCTM.fluxes)):
        SLCTM.fluxes[i] = SLCTM.fluxes[i] + np.random.normal(0,1)*lcnoise
    #Renormalize after noise:
    SLCTM.fluxes = np.divide(SLCTM.fluxes,np.amax(SLCTM.fluxes))

#Returns ChiSquared for a given set of data and a model, both interpolated.
#ASSUMPTION: Same number of points, same data range
def ComputeChiSqInter(DataSLCTM,ModelSLCTM):
    ChiSqComp = []
    #for ii in range(len(DataSLCTM.fluxes)):
    #    ChiSqComp.append((DataSLCTM.fluxes[ii] - ModelSLCTM.fluxes[ii])**2 / (ModelSLCTM.fluxes[ii]))
    ChiSqComp = (DataSLCTM.interpfluxes - ModelSLCTM.interpfluxes)**2 / (ModelSLCTM.interpfluxes)
    if DataSLCTM.isNoisy == False:
        print "WARNING: Simulated Data does not have Flux Noise!"
    return np.sum(ChiSqComp)

def InterpolateSLCTM(SLCTM,InterpolatedTimes):
    #Takes array of interpolated times, interpolates, saves to SLCTM object
    SLCTM.interpfluxes = np.interp(InterpolatedTimes,SLCTM.fluxtimes,SLCTM.fluxes)
    SLCTM.interpfluxtimes = InterpolatedTimes



##########################################################################################
#####TESTING##############################################################################
###################TESTING################################################################
######################################TESTING#############################################
##########################################################################################

#Testing Array of Light Curves with various TTV periods. Plots ChiSq(PTTV) vs PTTVs
#TTV Periods:
PTTVLowerBound = 50
PTTVUpperBound = 150
PTTVSegments = 200
PTTV_arr = np.linspace(PTTVLowerBound,PTTVUpperBound,PTTVSegments)
#Orbital Periods:
PorbLowerBound = 0.5
PorbUpperBound = 20
PorbSegments = 200
Porb_arr = np.linspace(PorbLowerBound,PorbUpperBound,PorbSegments)
#Interpolated Flux TimeStamps:
FTimeLow = 0
FTimeHigh= 100
FTimeSegments=1000000
InterpFluxTimes = np.linspace(FTimeLow,FTimeHigh,FTimeSegments)


#SLCTM / Batman Initial Params:
SLCTMInputParams = {
    't0': 0.0,
    'c1': 0.05,
    'c2':0.08,
    'porb':10.0,
    'pttv': (PTTVUpperBound - PTTVLowerBound) / 2.0,
    'noisett_e': 0.0001,
    'b':100.0,
    'vtan':200.0
}
BatmanInputParams = {
    't0': 0.0,
    'per':10.0,
    'rp':0.2,
    'a':12.0,
    'inc':90.0,
    'ecc':0.0,
    'w':90.0,
    'u1':0.1,
    'u2':0.3
}


#Batman Params:
bmparams = batman.TransitParams()
#SIMULATED DATA, with PTTV somewhere in middle of PTTV bounds:
DataLightCurve = SLCTM()
DataLightCurve.SetModelParams(SLCTMInputParams)
DataLightCurve.setBatmanParams(bmparams,BatmanInputParams)
PopTransTimes(DataLightCurve,int(10))
PopFluxesNaive(DataLightCurve,bmparams,500)
add_norm_lcnoise(DataLightCurve,0.005)
InterpolateSLCTM(DataLightCurve,InterpFluxTimes)

ChiSqs = []



#Compute Models:
#List of LC Models (Vary PTTV):
LightCurves = [SLCTM() for i in range(len(PTTV_arr))]
for i in range(len(PTTV_arr)):
    SLCTMInputParams['pttv'] = PTTV_arr[i]
    LightCurves[i].SetModelParams(SLCTMInputParams)
    LightCurves[i].setBatmanParams(bmparams,BatmanInputParams)
    PopTransTimes(LightCurves[i],int(10))
    PopFluxesNaive(LightCurves[i],bmparams,500)
    InterpolateSLCTM(LightCurves[i],InterpFluxTimes)
    #Compute ChiSqs for each model: #PROBLEM: Lines up with flux point index, NOT time index
    ChiSqs.append(ComputeChiSqInter(DataLightCurve,LightCurves[i]))

#List of LC Models (Vary Porb):
#LightCurves = [SLCTM() for i in range(len(Porb_arr))]
#for i in range(len(Porb_arr)):
#    SLCTMInputParams['porb'] = Porb_arr[i]
#    BatmanInputParams['per'] = Porb_arr[i]
#    LightCurves[i].SetModelParams(SLCTMInputParams)
#    LightCurves[i].setBatmanParams(bmparams,BatmanInputParams)
#    PopTransTimes(LightCurves[i],int(10))
#    PopFluxesNaive(LightCurves[i],bmparams,500)
#    InterpolateSLCTM(LightCurves[i],InterpFluxTimes)
#    #Compute ChiSqs for each model: #PROBLEM: Lines up with flux point index, NOT time index
#    ChiSqs.append(ComputeChiSqInter(DataLightCurve,LightCurves[i]))


#Plots ChiSq(PTTV) vs PTTVs
plt.plot(PTTV_arr,ChiSqs)
plt.title("$\chi^2$ vs $P_{ttv}$")
plt.xlabel("$P_{ttv}$ [Days]")
plt.ylabel("$\chi^2$")
plt.show()

#plt.plot(LightCurves[2].interpfluxtimes,LightCurves[2].interpfluxes)
#plt.title("Model Light Curve")
#plt.xlabel("Time [Days]")
#plt.ylabel("Normalized Flux")
#plt.show()

#plt.plot(LightCurves[100].interpfluxtimes,LightCurves[100].interpfluxes)
#plt.title("Model Light Curve")
#plt.xlabel("Time [Days]")
#plt.ylabel("Normalized Flux")
#plt.show()

#print Porb_arr
