# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 10:35:10 2023

@author: ztoth
"""
# set PYTHONPATH = %PYTHONPATH%; C:\Users\ztoth\Documents\Python\

import os.path
from winsound import Beep
import random
import pickle
import numpy as np
import pandas as pd

import tifffile as tf

from numpy import cos, sin, tan, pi

from PoreGenerator_funcs import nameFig
from PoreGenerator_funcs import getPD
from PoreGenerator_funcs import getRC
from PoreGenerator_funcs import networkPore
from PoreGenerator_funcs import getXY
from PoreGenerator_funcs import erodePores
from PoreGenerator_funcs import mergePores
from PoreGenerator_funcs import getTextOutput
from PoreGenerator_funcs import make3DModel

class Struct:
    pass


# %% Settings

option = Struct()
# set to True to keep the same rng used last generation. Any other value will be used as a seed
option.debug = 'GivenPorosity'
# 'Timestamp' or name to be used
option.namestyle = '20%_Porosity' # 'Timestamp'
# True if variation in diameter and circularity should be linked
# eg, large pored would tend to also be oblong
option.varLink = True
# mindiameter only in effect for linked diameter and circularity
option.mindiameter = 0

# 'Circle' (or 1) for circular grid, 'Square' (or 2) for square grid, else random
option.LocationType = 0
# If LocationType is not random, distance between grid lines
option.Spacing = 16  # 16
# If LocationType is not random, variation of pores from grid lines
# values 0-5 ignore TargetPorosity
option.location_err = 6  # 6
# If LocationType is square grid, whether to generate pores along the edge.
option.ignoreborder = False

# List of diameters to randomly choose from
# ignored if left empty
# Linear chooses from list randomly
# Weighted chooses from list according to distribution
option.LinearDiscreteDiameters = []
option.WeightedDiscreteDiameters = []
# List of circularities (values 0-1) to randomly choose from
# ignored if left empty
# Linear chooses from list randomly
# Weighted chooses from list according to distribution
option.LinearDiscreteCircularities = []
option.WeightedDiscreteCircularities = []

# True if nearby pores should be smoothed and merged
option.mergePores = True
# True if pores should have varying shape (cylinder,cone,ellipsoid,hyperboloid)
option.variedPoreShape = True

# 10 micrometers/voxel
option.ArraySize = 200  # 200

# 'Exp' to choose value from experimental distribution
TargetPorosity = 'Exp'  # 'Exp'

option.maxosteonlength = 220/3  # 220/3

export = Struct()
export.xcls = True
export.txt = False
export.tiff = True
export.stl = False

# %% Parameters

mu = Struct()
sigma = Struct()
weighting = Struct()

# SED related parameters below from DOI: 10.1002/jbmr.3561
# Parameters for SED per hole
mu.SED = 8  # 8
sigma.SED = 4.5  # 4.5
# SED weightings for combined distribution
weighting.SED = [0.9774, 0.0226]  # [0.9774 0.0226]

# Parameters for Normal SED
mu.Ndiameter = 7  # 4.8
sigma.Ndiameter = 0.4  # 0.4
mu.Ncircularity = 0.66  # 0.66
sigma.Ncircularity = 0.03  # 0.03

# Parameters for High SED
mu.Hdiameter = 18  # 15.7
sigma.Hdiameter = 3.75  # 3.75
mu.Hcircularity = 0.42  # 0.42
sigma.Hcircularity = 0.045  # 0.045

# Parameters for Canal Length
mu.osteonlength = 100  # 100
sigma.osteonlength = 75/2  # 75/2

# Parameters for Porosity
mu.porosity = 0.075046512  # 0.075046512
sigma.porosity = 0.036744908  # 0.036744908

# Parameters for phi value selection DOI: 10.1016/8756-3282(94)90288-7
weighting.phi_values = [0, pi/12, pi/2]  # [0,pi/12,pi/2]
weighting.phi_probs = [0.5, 0.5]  # [0.5, 0.5]

# Parameters for theta value selection
weighting.theta_values = [0,2*pi]  # [0,2*pi]
weighting.theta_probs = [1]  # [1]

# Proportions of each pore shape: [Cylinder, Proximal Cone, Distal Cone,
#                                                   Ellipsoid, Hyperboloid]
shape_proportions = [0.392, 0.094, 0.351, 0.122, 0.041]  # [0.392, 0.094, 0.351, 0.122, 0.041] DOI: 10.1111/j.1439-0264.2009.00973.x

pores_before_networking = 75  # 75
top_branches = [0,2] # [0,2]
bottom_branches = [0,2] # [0,2]
sealed_osteon_chance = 0.068  # 0.068 # DOI: 10.1002/ar.21309
transverse_flag_onset = pi/4  # pi/4

# %% Initialization

# set file name
[fpath, fname] = nameFig(option)

# sets rng based on option.debug
if option.debug == 1:
    if os.path.isfile(fpath+'saved_rng.pkl'):
        with open('saved_rng.pkl') as f:
            RNGkey = pickle.load(f)
        random.setstate(RNGkey)
elif option.debug == 0:
    random.seed()
    RNGkey = random.getstate()
    with open('saved_rng.pkl', 'wb') as f:
        pickle.dump(RNGkey, f)
else:
    rngseed = []
    for character in option.debug:
        rngseed.append(ord(character))
    random.seed(sum(rngseed))
    RNGkey = random.getstate()
    with open('saved_rng.pkl', 'wb') as f:
        pickle.dump(RNGkey, f)

# Initializes array of proper size
Bone = np.ones((option.ArraySize, option.ArraySize, option.ArraySize), dtype=np.float32)
# Sets up indexing of Bone array for use in calculations
[ii,ij,ik] = np.unravel_index(np.arange(option.ArraySize**3), [option.ArraySize, option.ArraySize, option.ArraySize], 'F')

PD = getPD(mu, sigma, weighting, option)

# Chooses TargetPorosity Value from experimental distribution
if TargetPorosity == 'Exp':
    TargetPorosity = PD.porosity.rvs(1)[0]

# readjusts TargetPorosity to account for loss when mergePores is used
if option.mergePores is True:
    TargetPorosity = TargetPorosity/0.739

# creates log for use in pore networking
valueslog = np.zeros([12, int(round(7000000*TargetPorosity/mu.osteonlength, ndigits=-3))])
    # [R; C; theta; phi; x; y; minz; z; maxz; isfilled; A; B]
iteration = 0

# Pre-sets some variables for choosing pore location
# these are set for pore(n+1) when pore(n) is generated
XYprimer = Struct()
if option.LocationType == 1:
    XYprimer.it = 1
    XYprimer.iu = 1
    XYprimer.AngleList = np.linspace(0,2*pi,XYprimer.it*int(np.sqrt(option.ArraySize/option.Spacing)))
else:
    XYprimer.it = 1+option.ignoreborder
    XYprimer.iu = 1+option.ignoreborder

XYprimer.SpaceList = np.linspace(0,option.ArraySize,int(option.ArraySize/option.Spacing))
XYprimer.ignoreTP = 0
XYprimer.griddone = 0

# %% Main Body

while ((1-np.mean(Bone) < TargetPorosity) and (XYprimer.ignoreTP == 0)) or (XYprimer.griddone == 0):

    [R, C] = getRC(option, PD)
    z = random.random() * option.ArraySize
    # theta is angle of trajectory, phi is angle of depression
    theta = PD.theta.rvs(1)
    phi = PD.phi.rvs(1)

    # Correction factor (phi is also correcting things in eqn)
    C = C/cos(phi)
    # cuts off cylinder
    minz = z - cos(phi)*0.5*PD.osteonlength.rvs(1)[0]
    maxz = z + cos(phi)*0.5*PD.osteonlength.rvs(1)[0]

    if sum(valueslog[10,:]) > pores_before_networking and random.random() > sealed_osteon_chance:
        [x,y,minz,z,maxz,valueslog,iteration] = networkPore(valueslog,minz,z,maxz,iteration)
    else:
        [x, y, XYprimer] = getXY(option, XYprimer)

    if phi > transverse_flag_onset:
        miny = y - sin(phi) * PD.osteonlength.rvs(1)[0] * 0.25
        maxy = y + sin(phi) * PD.osteonlength.rvs(1)[0] * 0.25
        minx = x - sin(phi) * PD.osteonlength.rvs(1)[0] * 0.25
        maxx = x + sin(phi) * PD.osteonlength.rvs(1)[0] * 0.25
    else:
        miny = 0
        minx = 0
        maxy = option.ArraySize
        maxx = option.ArraySize

    if option.variedPoreShape == 1:
        shapechoice = random.random()
        if shapechoice < shape_proportions[0]:
            shape = 1                                          # Cylinder
        elif shapechoice < sum(shape_proportions[0:2]):
            shape = 3*(5+abs(minz-ik))/245                    # Proximal Cone
        elif shapechoice < sum(shape_proportions[0:3]):
            shape = 3*(5+abs(maxz-ik))/245                    # Distal Cone
        elif shapechoice < sum(shape_proportions[0:4]):
            shape = 3*(5+np.minimum(abs(maxz-ik), abs(minz-ik)))/130  # Ellipsoid
        elif shapechoice <= 1:
            shape = 3*(5+abs(z-ik))/130                       # Hyperboloid
    else:
        shapechoice = 0
        shape = 1


    Bone = np.multiply(Bone, np.reshape(np.transpose(\
        ((cos(theta)**2 + sin(theta)**2 / C**2)*(ii-x+ (ik-z)*tan(phi)*sin(theta))**2 \
        + (2*sin(theta)*cos(theta)*(1- (1/ C**2))*(ii-x+ (ik-z)*tan(phi)*sin(theta)))*(ij-y- (ik-z)*tan(phi)*cos(theta)) \
        + (sin(theta)**2 + cos(theta)**2 / C**2)*(ij-y- (ik-z)*tan(phi)*cos(theta))**2 \
        > R**2 * shape * ((ii>minx) & (ii<maxx)) * ((ij>miny) & (ij<maxy)) * ((ik>minz) & (ik<maxz)) ) \
        ), Bone.shape, order = 'F'))

    A = random.randint(bottom_branches[0], bottom_branches[1])
    B = random.randint(top_branches[0], top_branches[1])

    if A+B != 0:
        valueslog[:,iteration] = [np.squeeze(i) for i in [R,C,theta,phi,x,y,minz,z,maxz,1,A,B]]
        iteration += 1

    # For Debugging Purposes Only
    # Uncomment following line to only generate one pore
    # XYprimer.griddone = 1; XYprimer.ignoreTP = 1;

if option.mergePores == 1:
    for n in range(10):
        if random.randint(0, 1) == 0:
            Bone = erodePores(Bone)
        else:
            Bone = mergePores(Bone)
    # reverts TargetPorosity for bookkeeping
    TargetPorosity = TargetPorosity*0.739

Bone = ~(Bone.astype(bool))

# %% Save Results

fullcell, sheetprep, sheetcell = getTextOutput(option, mu, sigma, weighting, TargetPorosity, pores_before_networking, \
                   sealed_osteon_chance, transverse_flag_onset, shape_proportions, RNGkey, fname)
if export.xcls is True:
    spreadsheet_name = 'Metadata.xlsx'
    if os.path.isfile(fpath+spreadsheet_name) is False:
        headerslist = {sheetprep[n]: '' for n in range(len(sheetprep))}
        headers = pd.DataFrame(headerslist, index=[0])
        headers.to_excel(fpath+spreadsheet_name, sheet_name='sheet1', index=False)

    current_params = {sheetprep[n]:sheetcell[0,n] for n in range(len(sheetprep))}
    new_row = pd.DataFrame(current_params, index=[0])
    # append to excel
    past_rows = pd.read_excel(fpath+spreadsheet_name)
    whole_sheet = pd.concat([past_rows, new_row], ignore_index=True)
    whole_sheet.to_excel(fpath+spreadsheet_name, sheet_name='sheet1', index=False)

if export.txt is True:
    txtfile = open(fpath+fname.removesuffix('.tif')+'.txt', "w")
    txtfile.writelines(fullcell)
    txtfile.close()

if export.tiff is True:
    tf.imsave(fpath+fname, Bone)

if export.stl is True:
    make3DModel(fpath, fname, Bone)


Beep(500,500)




