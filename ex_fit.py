#!/usr/bin/python
#This Example shows how to use emcee with the SLCTM class.

import math
import sys
import numpy as np
import batman
import matplotlib.pyplot as plt
import emcee

import SLCTM


#Start Modifiable Constants for Test:
#WARNING: Too many light curves + too many data points can crash computer.
PTTV_Actual = 100       #Set Simulated Data's Actual PTTV
PTTVLowerBound = 1      #Lower Bound of Models
PTTVUpperBound = 250    #Upper Bound of Models
PTTVSegments = 1000     #Number of Model Light Curves to Calculate
PTTV_arr = np.linspace(PTTVLowerBound,PTTVUpperBound,PTTVSegments)
DATAPOINTSPERTRANSIT = 200 #Number of Datapoints per transit
NUMBEROFTRANSITS = 80      #Number of Transits in Dataset / Each Model
#End Modifiable Constants for Test





#To Start, create simulated data and models. Same as in ex_fit.py example, see for more detail.
########################################################################################
########################################################################################
SLCTMInputParams = {
    't0': 0.0,
    'c1': 0.2,
    'c2':0.05,
    'porb':10.0,
    'pttv': PTTV_Actual,
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
bmparams = batman.TransitParams()
DataLightCurve = SLCTM.LightCurve()
DataLightCurve.SetModelParams(SLCTMInputParams)
DataLightCurve.setBatmanParams(bmparams,BatmanInputParams)
SLCTM.PopTransTimes(DataLightCurve,NUMBEROFTRANSITS)
SLCTM.PopFluxesNaive_Data(DataLightCurve,bmparams,DATAPOINTSPERTRANSIT)
SLCTM.Add_Norm_LCnoise(DataLightCurve,0.005)
print "----SLCTM emcee Test (ex_fit.py)----"
print "  o Simulated Data Generated"
LightCurves = [SLCTM.LightCurve() for i in range(len(PTTV_arr))]
ChiSqs = []
for i in range(len(PTTV_arr)):
    LightCurves[i].modelnumber = i
    SLCTMInputParams['pttv'] = PTTV_arr[i]
    LightCurves[i].SetModelParams(SLCTMInputParams)
    LightCurves[i].setBatmanParams(bmparams,BatmanInputParams)
    SLCTM.PopTransTimes(LightCurves[i], NUMBEROFTRANSITS)
    SLCTM.PopFluxesNaive_Model(LightCurves[i],DataLightCurve,bmparams)
    SLCTM.ComputeChiSqInter(DataLightCurve,LightCurves[i])
    ChiSqs.append(SLCTM.ComputeChiSqInter(DataLightCurve,LightCurves[i]))
    prog = (100 * i / len(PTTV_arr)) + 1.0
    print "  o Model Progress: /" + int(prog/5)*"/" + (20-int(prog/5))*" " + "/ " + str(prog)+ "% \r",
print ""
########################################################################################
########################################################################################
#Okay, Data and Models are generated. These are accessed by objects:
#  o DataLightCurve
#  o LightCurves[i] from i = 0 to i = PTTVSegments - 1
#The Chisquareds are accessed in ChiSqs[i] with the same indices as LightCurves[i].

#emcee / MCMC setup / Use:
