#!/usr/bin/python
#This example showcases the production of a simulated light curve and a large
#number of model light curves varying one parameter, which then have their
#ChiSquareds evaluated and plotted.

import math
import sys
import numpy as np
import batman
import matplotlib.pyplot as plt
import emcee

import SLCTM


#Start Modifiable Constants for Test:
#WARNING: Too many light curves + too many data points can crash computer.
PTTV_Actual = 80     #Set Simulated Data's Actual PTTV
PTTVLowerBound = 1      #Lower Bound of Models
PTTVUpperBound = 200    #Upper Bound of Models
PTTVSegments = 500     #Number of Model Light Curves to Calculate
PTTV_arr = np.linspace(PTTVLowerBound,PTTVUpperBound,PTTVSegments)

DATAPOINTSPERTRANSIT = 100 #Number of Datapoints per transit
NUMBEROFTRANSITS = 50      #Number of Transits in Dataset / Each Model
#End Modifiable Constants for Test


#SLCTM / Batman Initial Param Dictionaries:
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

#Optionally Print Settings to Commandline
print "----SLCTM ChiSq Test (ex_chisqs.py)----"
print "  o Varying P_ttv. N = " + str(len(PTTV_arr)) + " samples from " + str(PTTVLowerBound) + " to " + str(PTTVUpperBound) + " days. Data's PTTV is: " + str(PTTV_Actual) + " days."

bmparams = batman.TransitParams() #Initialize Batman's native parameter object.
                                  #This is used for both the simulated data and
                                  #each model. Held constant under this test.

#Create Simulated Data / Light Curve
DataLightCurve = SLCTM.LightCurve() #Initialize a light curve for the simulated data.
DataLightCurve.SetModelParams(SLCTMInputParams) #Set the Parameters by passing our parameter dictionary.
DataLightCurve.setBatmanParams(bmparams,BatmanInputParams) #Set the batman parameters by passing our batman parameter dictionary.
SLCTM.PopTransTimes(DataLightCurve,NUMBEROFTRANSITS) #Populate NUMBEROFTRANSITS transit times. This fills an array full of transit times.
SLCTM.PopFluxesNaive_Data(DataLightCurve,bmparams,DATAPOINTSPERTRANSIT) #This populates each transit's light curve with a given number of points and generates the list of flux time points.
SLCTM.Add_Norm_LCnoise(DataLightCurve,0.005) #This adds a simulated amount of noise to each data point. Sets LightCurve.isNoisy to true.
#End of Processing Simulated Data
print "  o Simulated Data Generated"


#Compute Models:
LightCurves = [SLCTM.LightCurve() for i in range(len(PTTV_arr))] #Create a Model Light Curve for each PTTV Time in PTTV_arr (same as PTTVSegments)
ChiSqs = []
#For each model...
for i in range(len(PTTV_arr)):
    LightCurves[i].modelnumber = i            #Set model number
    SLCTMInputParams['pttv'] = PTTV_arr[i]    #Set this model to an element of PTTV_arr
    LightCurves[i].SetModelParams(SLCTMInputParams) #Set the (revised) Parameters by passing our parameter dictionary.
    LightCurves[i].setBatmanParams(bmparams,BatmanInputParams) #Set the batman parameters by passing our batman parameter dictionary.
    SLCTM.PopTransTimes(LightCurves[i], NUMBEROFTRANSITS) #Populate NUMBEROFTRANSITS transit times. This fills an array full of transit times.
    SLCTM.PopFluxesNaive_Model(LightCurves[i],DataLightCurve,bmparams) #This populates each transit's light curve using same time points as the data's light curve.
    SLCTM.ComputeChiSqInter(DataLightCurve,LightCurves[i])    #Compute ChiSqs for each model:
    ChiSqs.append(SLCTM.ComputeChiSqInter(DataLightCurve,LightCurves[i])) #Append calculated ChiSq to master array

    #Optionally print out the progress. This can take a while so it's nice to see something working.
    prog = (100 * i / len(PTTV_arr)) + 1.0
    print "  o Model Progress: /" + int(prog/5)*"/" + (20-int(prog/5))*" " + "/ " + str(prog)+ "% \r",
print ""

#Plots ChiSq(PTTV) vs PTTVs
plt.plot(PTTV_arr,ChiSqs, color = 'k')
plt.axvline(x=PTTV_Actual, linewidth=1, color='r')
plt.title("$\chi^2$ vs $P_{TTV}$ for a Simulated $P_{TTV}$ = " + str(PTTV_Actual) + " days")
plt.xlabel("$P_{TTV}$ [Days]")
plt.ylabel("$\chi^2$")
plt.show()

print ".....SLCTM Terminated....."
