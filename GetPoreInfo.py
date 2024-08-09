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
from skimage import io
import tifffile as tf

from ccl.ccl import connected_component_labelling


#%% Get Array from TIFF file

def importBone(fpath, fname):
    '''
    Imports 'Bone' array back into python for processing
    '''
    # reads TIFF file
    rawBone = io.imread(fpath + fname)
    # if TIFF file is RGB, flattens to B/W
    if rawBone.ndim == 4:
        Bone_bool = np.squeeze(np.delete(rawBone, (1,2), 3))
        Bone_bool = Bone_bool.astype(bool)
    elif rawBone.ndim == 3:
        Bone_bool = rawBone.astype(bool)
    # if image has [1] for matrix and [0] for pore, switches it
    if np.mean(Bone_bool) > 0.5:
        Bone = np.float32(~Bone_bool)
    else:
        Bone = np.float32(Bone_bool)
    
    return Bone

# %% Functions for modifying images for 2D viewing in ImageJ

def rotateBone(fpath, fname):
    '''Rotates bone file so the view plane is along the long axis '''
    Bone = importBone(fpath, fname)
    Bone = np.rot90(Bone, axes=(0,2))
    tf.imsave(fpath+fname, Bone)
    
def invertBone(fpath, fname):
    '''Inverts bone file to swap black / white locations'''
    Bone = importBone(fpath, fname)
    Bone = ~(Bone.astype(bool))
    Bone = np.float32(Bone)
    tf.imsave(fpath+fname, Bone)

# %%  Get Data for .tif files in folder

def GetPoreData(file_directory = './pore_files/'):
    '''
    Runs through file directory, cataloguing all .tif files
    
    'counts' is list of lists of area of each pore
    'values' is list of all observed pore diameters
    '''
    # gets list of all .tif files in directory
    fname = [f for f in os.listdir(file_directory) if f.endswith('.tif')]
    
    counts = []
    values = []
    
    # for each file
    for f in fname:
        Bone = importBone(file_directory, f)
        # for each slice of file, label and get pixel counts for each 2Dpore
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
    
    if os.name == 'nt':
        from winsound import Beep
        Beep(500, 500)  # when done
        
    return counts, values

#%% Save Pore Data

def SavePoreData(counts, values, filename):
    '''Saves data from GetPoreData as pkl file'''

    fpath = "./stats_datasets/"
    if os.path.isdir(fpath) is False:
        os.makedirs(fpath)
    with open(fpath + filename, 'wb') as f:
        try:
            pickle.dump([counts, values], f)
        except NameError:
            print("Variable 'counts' or 'values' is not defined")
    
#%% Read Pore Data

def ReadPoreData(targetfile):
    '''Imports values dataset from GetPoreData back in'''
    
    targetfile = "./stats_datasets/" + targetfile
    if (os.path.isfile(targetfile) & (os.path.getsize(targetfile) > 0)):
        with open(targetfile, 'rb') as f:
            *_, porevalues = pickle.load(f)
    return porevalues

# %% Create Histogram

def PoreHist(values, title):
    '''Creates histogram of pore diameters (0–60 µm, 20 bins)'''
    
    pd.Series(values).plot.hist(grid=True, bins=20, range=(0,60), rwidth=0.9, density=True, figure=plt.figure())
    plt.title(title)
    plt.xlabel('Diameter')
    plt.ylabel('Proportion of Total Pores')
    plt.grid(axis='y', alpha=0.75)
    plt.xlim([0,60])
    plt.ylim([0,0.175])

#%% Plot Diameter Graphs

def PlotDiameterGraphs(targetfiles, titles):
    '''Creates histogram of pore diameters directly from PoreData file'''
    
    if type(targetfiles) is str: targetfiles = [targetfiles]
    if type(titles) is str: titles = [titles]
    if len(titles) != len(targetfiles):
        raise Exception('Number of titles and data files to plot do not match')
            
    for a,b in zip(targetfile, titles):
        PoreHist(ReadPoreData(a), b)

#%% Percent Bin Height Change

def PercentBinHeightChange(values_1, values_2):
    '''Gets the average percent change in bin heights between samples'''
    
    values_1_incidence = []
    values_2_incidence = []
    
    range_ = 60
    bins = 20
    
    for b in range(bins):
        values_1_incidence.append(len(list(x for x in values_1 if range_*b/bins <= x < range_*(b+1)/bins))/len(values_1))
        values_2_incidence.append(len(list(x for x in values_2 if range_*b/bins <= x < range_*(b+1)/bins))/len(values_2))
    
    
    delta = 100*abs(np.subtract(values_2_incidence, values_1_incidence))
    
    average = np.mean(delta)
    sigma = np.std(delta)
    
    print("Change Sample 1 ->  Sample 2:       "+str(round(average,ndigits=2))+"±"+str(round(sigma,ndigits=2))+"%")

# %% Gets graphs and bin height changes for all datasets
# This serves as an example for how to use these functions
# Note: This uses my (author's) datasets

if __name__ == "__main__":
    
    targetfile = ['RealPoreInfo.pkl',
                  'GenPoreInfo.pkl',
                  '7-18PoreInfo.pkl',
                  '7-15.7PoreInfo.pkl',
                  'MatchingPoreInfo.pkl',
                  '7-25PoreInfo.pkl']
            
    titles = ['Pore Diameter from CT Scans',
              'Pore Diameter from Literature Values',
              'Pore Diameter from Optimized Values',
              '7-15.7 Pore Diameter',
              'Matching Pore Diameter',
              '7-25 Pore Diameter']

    PlotDiameterGraphs(targetfile, titles)
    
    # Get bin height change between Real and Literature / Optimized datasets
    print("Real Sample -> Values from Literature:")
    PercentBinHeightChange(ReadPoreData(targetfile[0]), ReadPoreData(targetfile[1]))
    print("Real Sample -> Optimized Values:")
    PercentBinHeightChange(ReadPoreData(targetfile[0]), ReadPoreData(targetfile[2]))