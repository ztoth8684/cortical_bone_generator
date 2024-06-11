# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 23:04:06 2023

@author: ztoth
"""
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import os
import cProfile
import pstats
from skimage import io
from winsound import Beep

from ccl_gh_pages.ccl import connected_component_labelling

class Struct:
    pass

#%% Get Array from TIFF file

def importBone(fpath, fname):
    # reads TIFF file
    rawBone = io.imread(fpath + fname)
    # if TIFF file is RGB, flattens to B/W
    if rawBone.ndim == 4:
        Bone_bool = np.squeeze(np.delete(rawBone, (1,2), 3))
    elif rawBone.ndim == 3:
        Bone_bool = rawBone
    # if image has [1] for matrix and [0] for pore, switches it
    if np.mean(Bone_bool) > 0.5:
        Bone = np.float32(~Bone_bool)
    else:
        Bone = np.float32(Bone_bool)
    
    return Bone

# %%  Get Data for ,tif files in folder

def GetPoreData(file_directory = None):

    if file_directory is None: 
        file_directory = "./pore_files/bones_7-25/"
    if file_directory == 'Real': 
        file_directory = "./pore_files/real_bone/" 
    
    fname = [f for f in os.listdir(file_directory) if f.endswith('.tif')]
    # fname = ['pytest.tif']
    
    counts = []
    values = []
    
    for f in fname:
        Bone = importBone(file_directory, f)
        # for each slice, label and get pixel counts for each 2Dpore
        for n in range(len(Bone)):    
            layer = Bone[n,:,:]
            result = connected_component_labelling(layer, 8)
            counts.append(np.unique(result, return_counts=True)[1])
    
    # for each 2Dpore, get diameter and add to values list
    for n in range(len(counts)):
        for m in range(1,len(counts[n])):
            A = counts[n][m]
            R = np.sqrt(A/np.pi)
            values.append(round(2*R))
                
    values.sort(reverse=True)
    
    Beep(500, 500)  # when done
    return counts, values

#%% Save Pore Data

def SavePoreData(filename):
    with open("./stats_datasets/" + filename, 'wb') as f:
        try:
            pickle.dump([counts, values], f)
        except NameError:
            print("Variable 'counts' or 'values' is not defined")
    
#%% Read Pore Data

def ReadPoreData(targetfile):
    targetfile = "./stats_datasets/" + targetfile
    if (os.path.isfile(targetfile) & (os.path.getsize(targetfile) > 0)):
        with open(targetfile, 'rb') as f:
            *_, porevalues = pickle.load(f)
    return porevalues

# %% Create Histogram

def PoreHist(values, title):
    pd.Series(values).plot.hist(grid=True, bins=20, range=(0,60), rwidth=0.9, density=True, figure=plt.figure())
    plt.title(title)
    plt.xlabel('Diameter')
    plt.ylabel('Proportion of Total Pores')
    plt.grid(axis='y', alpha=0.75)
    plt.xlim([0,60])
    plt.ylim([0,0.175])

#%% Plot Diameter Graphs

def PlotDiameterGraphs(titles=None, varnames=None):
    if titles is None: 
        titles = ['CT Scan Pore Diameter','4.8-15.7 Pore Diameter','7-18 Pore Diameter','7-15.7 Pore Diameter','Matching Pore Diameter','7-25 Pore Diameter']
    if varnames is None: 
        varnames = ['values_Real','values_48_157','values_7_18','values_7_157','values_Matching','values_7_25']
    if type(titles) is str: titles = [titles]
    if type(varnames) is str: varnames = [varnames]
    if len(titles) != len(varnames):
        raise Exception

    for c,b in zip(titles, varnames):
        if b in locals():
            exec("PoreHist(%s, c)" % (b))

#%% Unbalanced T-Test

def UnbalancedTTest(values_test, values_real):
    
    variance_test = pd.DataFrame(values_test).var()[0]
    variance_real = pd.DataFrame(values_real).var()[0]
    
    mu_real = np.mean(values_real)
    mu_test = np.mean(values_test)
    sigma_real = np.std(values_real)
    sigma_test = np.std(values_test)
    n_real = len(values_real)
    n_test = len(values_test)
    ssqn_real = np.square(sigma_real)/n_real
    ssqn_test = np.square(sigma_test)/n_test
    
    t = (mu_real - mu_test)/np.sqrt(ssqn_real + ssqn_test)
    
    df = np.square(ssqn_real + ssqn_test)/((np.square(ssqn_real)/(n_real - 1)) + (np.square(ssqn_test)/(n_test - 1)))
    
    return t, df, variance_real, variance_test

#%% Percent Bin Height Change

def PercentBinHeightChange(Best_values = None, Literature_values = None, Real_values = None):
    # Best_values : diameter list from parameters optimized for closeness to CT
    if Best_values is None: 
        Best_values = ReadPoreData('7-18PoreInfo.pkl')
    # Literature_values: diameter list from parameters from the literature
    if Literature_values is None: 
        Literature_values = ReadPoreData('GenPoreInfo.pkl')
    # Real_values : diameter list from CT scans
    if Real_values is None: 
        Real_values = ReadPoreData('RealPoreinfo.pkl')
    
    Real_incidence = []
    Best_incidence = []
    Literature_incidence = []
    
    range_ = 60
    bins = 20
    
    for b in range(bins):
        Real_incidence.append(len(list(x for x in Real_values if range_*b/bins <= x < range_*(b+1)/bins))/len(Real_values))
        Best_incidence.append(len(list(x for x in Best_values if range_*b/bins <= x < range_*(b+1)/bins))/len(Best_values))
        Literature_incidence.append(len(list(x for x in Literature_values if range_*b/bins <= x < range_*(b+1)/bins))/len(Literature_values))
    
    
    Delta_literature = 100*abs(np.subtract(Real_incidence, Literature_incidence))
    Delta_best = 100*abs(np.subtract(Real_incidence, Best_incidence))
    
    Best_dev = np.mean(Delta_best)
    Best_sigma = np.std(Delta_best)
    Lit_dev = np.mean(Delta_literature)
    Lit_sigma = np.std(Delta_literature)
    
    print("Change Real Sample -> Values from Literature: "+str(round(Lit_dev,ndigits=2))+"±"+str(round(Lit_sigma,ndigits=2))+"%")
    print("Change Real Sample -> Optimized Values:       "+str(round(Best_dev,ndigits=2))+"±"+str(round(Best_sigma,ndigits=2))+"%")

#%% Runtime Profiler

def RuntimeProfiler(filename = 'PoreGenerator_MK6.py'):
    cProfile.run(filename, 'file')
    p = pstats.Stats('file')
    p.sort_stats('cumulative').print_stats(10)

# %% Import all values lists from .pkl files
# Note: This is hacky and for testing purposes

if __name__ == "__main__":
    targetfile = ['RealPoreInfo.pkl','GenPoreInfo.pkl','7-18PoreInfo.pkl','7-15.7PoreInfo.pkl','MatchingPoreInfo.pkl','7-25PoreInfo.pkl']
    varnames = ['values_Real','values_48_157','values_7_18','values_7_157','values_Matching','values_7_25']
    
    targetfile = list(map(lambda x: "./stats_datasets/" + x, targetfile))
    
    for a,b in zip(targetfile, varnames):
        if (os.path.isfile(a) & (os.path.getsize(a) > 0)):
            a = a.removeprefix("./stats_datasets/")
            exec("%s = ReadPoreData(a)" % (b))
    del a, b, targetfile, varnames
