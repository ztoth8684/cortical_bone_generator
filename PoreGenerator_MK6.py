# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 10:35:10 2023

@author: ztoth
"""
# set PYTHONPATH = %PYTHONPATH%; C:\Users\ztoth\Documents\Python\

import os.path
import sys
from winsound import Beep
import random
import pickle
import numpy as np
import pandas as pd
import tifffile as tf

from numpy import cos, sin, tan, pi

from PoreGenerator_funcs import nameFig, getPD, getRC, networkPore, getXY, \
    poreBlast, poreClast, getTextOutput, make3DModel
import PoreGenerator_classes as PGc

from LoadParameters import LoadParameters

class Struct:
    pass


# %% Initialization

option, target_porosity, export, mu, sigma, weighting, params = LoadParameters()


# set file name
[fpath, fname] = nameFig(option)

# sets rng based on option.debug
if option.debug == 1:
    if os.path.isfile('./saved_rng.pkl'):
        with open('saved_rng.pkl') as f:
            RNGkey = pickle.load(f)
        random.seed(RNGkey)
elif option.debug == 0:
    RNGkey = random.randrange(sys.maxsize)
    random.seed(RNGkey)
    with open('saved_rng.pkl', 'wb') as f:
        pickle.dump(RNGkey, f)
else:
    rngseed = []
    for character in option.debug:
        rngseed.append(ord(character))
    RNGkey = sum(rngseed)
    random.seed(RNGkey)
    with open('saved_rng.pkl', 'wb') as f:
        pickle.dump(RNGkey, f)

# Initializes array of proper size
Bone = np.ones((option.ArraySize, option.ArraySize, option.ArraySize), dtype=np.float32)
# Sets up indexing of Bone array for use in calculations
[ii,ij,ik] = np.unravel_index(np.arange(option.ArraySize**3), [option.ArraySize, option.ArraySize, option.ArraySize], 'F')

PD = getPD(mu, sigma, weighting, option)

# Chooses target_porosity Value from experimental distribution
if target_porosity == 'Exp':
    target_porosity = PD.porosity.rvs(1)[0]

# readjusts target_porosity to account for loss when mergePores is used
if option.smoothPores is True:
    target_porosity = target_porosity/0.739

# creates log for use in pore networking
valueslog = np.zeros([12, int(round(7_000_000*target_porosity/mu.osteonlength, ndigits=-3))])
    # [R; C; theta; phi; x; y; minz; z; maxz; isfilled; A; B]
iteration = 0

# Pre-sets some variables for choosing pore location
XYprimer = PGc.XYprimer(option)


# %% Main Body

while ((1-np.mean(Bone) < target_porosity) and (XYprimer.ignore_target_porosity == 0)) or (XYprimer.grid_complete == 0):

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

    if sum(valueslog[10,:]) > params.pores_before_networking and random.random() > params.sealed_osteon_chance:
        [x,y,minz,z,maxz,valueslog,iteration] = networkPore(valueslog,minz,z,maxz,iteration)
    else:
        [x, y, XYprimer] = getXY(option, XYprimer)

    if phi > params.transverse_flag_onset:
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
        if shapechoice < params.shape_proportions[0]:
            shape = 1                                          # Cylinder
        elif shapechoice < sum(params.shape_proportions[0:2]):
            shape = 3*(5+abs(minz-ik))/245                    # Proximal Cone
        elif shapechoice < sum(params.shape_proportions[0:3]):
            shape = 3*(5+abs(maxz-ik))/245                    # Distal Cone
        elif shapechoice < sum(params.shape_proportions[0:4]):
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

    A = random.randint(params.bottom_branches[0], params.bottom_branches[1])
    B = random.randint(params.top_branches[0], params.top_branches[1])

    if A+B != 0:
        valueslog[:,iteration] = [np.squeeze(i) for i in [R,C,theta,phi,x,y,minz,z,maxz,1,A,B]]
        iteration += 1

    # For Debugging Purposes Only
    # Uncomment following line to only generate one pore
    # XYprimer.grid_complete = 1; XYprimer.ignore_target_porosity = 1;

if option.smoothPores == 1:
    for n in range(10):
        if random.randint(0, 1) == 0:
            Bone = poreBlast(Bone)
        else:
            Bone = poreClast(Bone)
    # reverts target_porosity for bookkeeping
    target_porosity = target_porosity*0.739

Bone = ~(Bone.astype(bool))

# %% Save Results

fullcell, sheetprep, sheetcell = getTextOutput(option, mu, sigma, weighting, params, target_porosity, RNGkey, fname)
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




