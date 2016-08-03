#!/usr/bin/python
import math
import sys
import numpy as np
import batman #Limb Darkening Model from
              #http://astro.uchicago.edu/~kreidberg/batman/quickstart.html
import emcee
#SLCTM: Simple Light Curve Transit Model with Quadratic Limb Darkening
#Outputs a file with transit times with noise and a file with times and fluxes


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
class LightCurve:
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

def Add_Norm_LCnoise(SLCTM,lcnoise):
    #Note that this is simulated data:
    SLCTM.isNoisy = True
    #Apply Noise to Fluxes
    for i in range(0,len(SLCTM.fluxes)):
        SLCTM.fluxes[i] = SLCTM.fluxes[i] + np.random.normal(0,1)*lcnoise
    #Renormalize after noise:
    SLCTM.fluxes = np.divide(SLCTM.fluxes,np.amax(SLCTM.fluxes))

def ComputeChiSqInter(DataSLCTM,ModelSLCTM):
    #Returns ChiSquared for a given set of data and a model, both interpolated.
    #ASSUMPTION: Same number of points, same data range
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

def LogP():
    return 0
def RunSampler():
    return 0
