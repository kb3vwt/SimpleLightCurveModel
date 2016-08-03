#!/usr/bin/python
#This example showcases the production of a light curve inside the SLCTM.LightCurve class.

import math
import sys
import numpy as np
import batman
import matplotlib.pyplot as plt
import emcee

import SLCTM


#Start Modifiable Constants for Test:
#WARNING: Too many data points can crash computer.
DATAPOINTSPERTRANSIT = 200 #Number of Datapoints per transit
NUMBEROFTRANSITS = 2      #Number of Transits in Dataset / Each Model
#End Modifiable Constants for Test


#SLCTM / Batman Initial Param Dictionaries:
SLCTMInputParams = {
    't0': 0.0,
    'c1': 0.2,
    'c2':0.05,
    'porb':10.0,
    'pttv': 100,
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


#Plots ChiSq(PTTV) vs PTTVs
plt.plot(DataLightCurve.fluxTimes,DataLightCurve.fluxes, color = 'k')
plt.title("Light Curve Example")
plt.xlabel("Time Stamp [Days]")
plt.ylabel("Normalized Flux")
plt.show()
