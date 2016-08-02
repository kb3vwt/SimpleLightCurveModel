#!/usr/bin/python

#SLCTM: Simple Light Curve Transit Model with Quadratic Limb Darkening
#Outputs a file with transit times with noise and a file with times and fluxes
import math
import sys
import numpy as np
import batman #Limb Darkening Model from
              #http://astro.uchicago.edu/~kreidberg/batman/quickstart.html

#For Testing, can remove later:
import matplotlib.pyplot as plt

#######################
#Things to Do:
#o Implement Error on fluxes
#o Multiplanet (vectorize Porb,...)
#o Possibly adaptively assign # of datapoints to a transit (save RAM)
#  s.t. time between points << transit duration
#o Make Model Compute on Data's Timestamps
#o Implement emcee model fit
#######################









#######CLASSES/DATA STRUCTURES##########################################################
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

        #SelfIdentification
        self.modelnumber = 0          #Model Number for self identification (optional)

        #Data Arrays:
        self.transTimes = []          #Transit Times
        self.fluxTimes = []           #Timestamp of Flux Datapoint
        self.fluxes = []              #Fluxes 1D Array - either a model or Simulated Data
        #self.interpfluxTimes = []     #Interpolated Flux Timestamps
        #self.interpfluxes = []        #Interpolated Fluxes
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
        SLCTM.transTimes.append(np.random.normal(0,1)*SLCTM.noisett_e+(SLCTM.t0 + i*SLCTM.porb + SLCTM.c1*math.sin((SLCTM.t0 + i*SLCTM.porb)/(2*math.pi*SLCTM.pttv))+ SLCTM.c2*math.cos((SLCTM.t0 + i*SLCTM.porb)/(2*math.pi*SLCTM.pttv))))

def PopFluxesNaive_Data(SLCTM,BatmanParams,NFluxPoints):
    #Populates fluxes for each transit. Generates Flux Timestamps
    #Naive - uses length of transit in plane ( inc = 0 deg )

    #Flux Times for This transit:
    #Calculate Interval duration:
    TIMEFACTOR = 2.0 #Use 1.0 if not calculating time outside occlusion
    HalfTransitAngle = 2.0 * math.asin((BatmanParams.rp + 1) / (2.0 * BatmanParams.a))
    HalfTransitTime = TIMEFACTOR*(HalfTransitAngle / (2.0 * math.pi)) * SLCTM.porb

    #Do for each transit...
    for i in range(0,SLCTM.NTransits):

        #Write Start Index:
        if len(SLCTM.fluxTimes) == 0:
            SLCTM.transitStartIndex.append(0)
        else:
            SLCTM.transitStartIndex.append(len(SLCTM.fluxTimes) - 1)

        #For ith transit, calculate fluxes using Batman
        BatmanParams.t0 = SLCTM.transTimes[i]                   #Set t0 for this transit:
        BMfluxTimes = np.linspace(SLCTM.transTimes[i] - HalfTransitTime, SLCTM.transTimes[i] + HalfTransitTime, NFluxPoints) #Temporary array of times for computation
        SLCTM.fluxTimes.extend(BMfluxTimes)                     #append times to master flux times:
        bmmodel = batman.TransitModel(BatmanParams,BMfluxTimes) #Initialize Batman model, compute fluxes for above interval
        BMFluxes = bmmodel.light_curve(BatmanParams)            #Temporary flux array
        SLCTM.fluxes.extend(BMFluxes)                           #Append new fluxes to master flux array:

        #print "Data, T: "+str(i)+" Flux Times... LEN:" + str(len(BMfluxTimes))
        #print BMfluxTimes

        #Write End Index:
        if len(SLCTM.fluxTimes) == 0:
            SLCTM.transitEndIndex.append(0)
        else:
            SLCTM.transitEndIndex.append(len(SLCTM.fluxTimes) - 1)

def PopFluxesNaive_Model(MSLCTM,DSLCTM,BatmanParams):
    #Populates fluxes for each transit. Uses given array of times, number of flux points per transit
    #MSLCTM - Model SLCTM
    #DSLCTM - Data SLCTM (Read only) - copies times so to skip interpolation
    #Naive - Assumes times and fluxes are in two 1D arrays, no overlap

    #Copy Start/End Times for Transits:
    MSLCTM.transitStartIndex = DSLCTM.transitStartIndex
    MSLCTM.transitEndIndex = DSLCTM.transitEndIndex #End Indices minus 1

    #Flux Times for This transit:
    #For ith transit, calculate fluxes using Batman
    for i in range(MSLCTM.NTransits):
        BatmanParams.t0 = MSLCTM.transTimes[i] #Set t0 for this transit:
        ThisTransitStart = DSLCTM.transitStartIndex[i]
        ThisTransitEnd   = DSLCTM.transitEndIndex[i]
        if i != 0:
            ThisTransitStart = ThisTransitStart + 1
        ThisTransitEnd = ThisTransitEnd + 1
        BMfluxTimes = DSLCTM.fluxTimes[ThisTransitStart:ThisTransitEnd] #Temporary array of times for computation
        MSLCTM.fluxTimes.extend(BMfluxTimes)   #append times to master flux times:
        BMfluxTimes = np.asarray(BMfluxTimes)  #Convert Flux Times to ndarray from list.
        bmmodel = batman.TransitModel(BatmanParams,BMfluxTimes)  #Initialize Batman model, compute fluxes for above interval
        BMFluxes = bmmodel.light_curve(BatmanParams)             #Temporary flux array
        MSLCTM.fluxes.extend(BMFluxes)                           #Append new fluxes to master flux array:

        #print "M: "+str(MSLCTM.modelnumber)+" T: "+str(i)+" Flux Times... LEN:" + str(len(BMfluxTimes))
        #print BMfluxTimes



def Add_Norm_LCnoise(SLCTM,lcnoise):
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
    if DataSLCTM.isNoisy == False:
        print "WARNING: Simulated Data does not have Flux Noise!"
        print "         Expect Unexpected Results!!"
    else:
        if len(DataSLCTM.fluxes) == len(ModelSLCTM.fluxes):
            ChiSqComp = (DataSLCTM.fluxes - ModelSLCTM.fluxes)**2 / (ModelSLCTM.fluxes)
        else:
            print "ERROR: DataSLCTM / ModelSLCTM Length Mismatch!"
            print "     --> data length: " + str(len(DataSLCTM.fluxes)) + ";   Model " + str(ModelSLCTM.modelnumber) + " Length: " + str(len(ModelSLCTM.fluxes))
    return np.sum(ChiSqComp)

#def InterpolateSLCTM(SLCTM,InterpolatedTimes):
#    #Takes array of interpolated times, interpolates, saves to SLCTM object
#    SLCTM.interpfluxes = np.interp(InterpolatedTimes,SLCTM.fluxTimes,SLCTM.fluxes)
#    SLCTM.interpfluxTimes = InterpolatedTimes










##TESTS##############################################################
#####################################################################


#Testing Array of Light Curves with various TTV periods. Plots ChiSq(PTTV) vs PTTVs
#TTV Periods:
PTTVLowerBound = 50
PTTVUpperBound = 150
PTTVSegments = 200
PTTV_arr = np.linspace(PTTVLowerBound,PTTVUpperBound,PTTVSegments)
#Orbital Periods:
PorbLowerBound = 8
PorbUpperBound = 12
PorbSegments = 200
Porb_arr = np.linspace(PorbLowerBound,PorbUpperBound,PorbSegments)


#SLCTM / Batman Initial Params:
SLCTMInputParams = {
    't0': 0.0,
    'c1': 0.3,
    'c2':0.01,
    'porb':10.0,
    'pttv': (PTTVUpperBound + PTTVLowerBound) / 2.0,
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

#Number of Datapoints per transit:
DATAPOINTSPERTRANSIT = 500
NUMBEROFTRANSITS = 4

#Batman Params:
bmparams = batman.TransitParams()
#SIMULATED DATA, with PTTV somewhere in middle of PTTV bounds:
DataLightCurve = SLCTM()
DataLightCurve.SetModelParams(SLCTMInputParams)
DataLightCurve.setBatmanParams(bmparams,BatmanInputParams)
PopTransTimes(DataLightCurve,NUMBEROFTRANSITS)
print "SLCTM Test in Progress..."
print "Varying P_ttv. N = " + str(len(PTTV_arr)) + " samples from " + str(PTTVLowerBound) + " to " + str(PTTVUpperBound) + " days."
print "------------------------"
print "Data Progress:"
print "    o Transit Times Populated."
PopFluxesNaive_Data(DataLightCurve,bmparams,DATAPOINTSPERTRANSIT)
print "    o Fluxes Populated."
Add_Norm_LCnoise(DataLightCurve,0.005)
print "    o Noise Added."
print "    o PTTV of SimuData:" + str((PTTVUpperBound + PTTVLowerBound) / 2.0) + " days"
print "------------------------"

ChiSqs = []



#Compute Models:
#List of LC Models (Vary PTTV):
print "MODELS: "
LightCurves = [SLCTM() for i in range(len(PTTV_arr))]
modelscalculated = 0
for i in range(len(PTTV_arr)):
    #Set model number:
    LightCurves[i].modelnumber = i
    SLCTMInputParams['pttv'] = PTTV_arr[i]
    LightCurves[i].SetModelParams(SLCTMInputParams)
    LightCurves[i].setBatmanParams(bmparams,BatmanInputParams)
    PopTransTimes(LightCurves[i], NUMBEROFTRANSITS)
    PopFluxesNaive_Model(LightCurves[i],DataLightCurve,bmparams)
    #Compute ChiSqs for each model:
    ComputeChiSqInter(DataLightCurve,LightCurves[i])
    ChiSqs.append(ComputeChiSqInter(DataLightCurve,LightCurves[i]))
    modelscalculated = modelscalculated + 1
print "    o Ran " + str(modelscalculated) + " out of " + str(PTTVSegments) + " Models."
#print "Model Flux Times:"
#print LightCurves[0].fluxTimes

#print "Data Flux Times:"
#print DataLightCurve.fluxTimes

#print type(DataLightCurve.transitEndIndex)

#print LightCurves[1].transitStartIndex[3]
#print LightCurves[1].transitEndIndex[3]
#print LightCurves[1].fluxTimes[LightCurves[1].transitStartIndex[3]+1:LightCurves[1].transitEndIndex[3]]


#Plots ChiSq(PTTV) vs PTTVs
plt.plot(PTTV_arr,ChiSqs)
plt.title("$\chi^2$ vs $P_{TTV}$")
plt.xlabel("$P_{TTV}$ [Days]")
plt.ylabel("$\chi^2$")
plt.show()

#plt.plot(LightCurves[2].fluxTimes,LightCurves[2].fluxes)
#plt.title("Model Light Curve")
#plt.xlabel("Time [Days]")
#plt.ylabel("Normalized Flux")
#plt.show()

#print Porb_arr
