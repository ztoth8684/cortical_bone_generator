# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:09:28 2023

@author: ztoth
"""

import os.path
import sys
import datetime
import pickle
import numpy as np
import scipy
import scipy.stats as stats
import random
import meshlib.mrmeshpy as mr
import meshlib.mrmeshnumpy as mrn

from PoreGenerator_classes import MixtureModel
class Struct:
    pass

#%%

def normalizeParameters(option, mu, sigma):
    '''Converts inputs from µm to pixel units'''
    
    option.mindiameter /= 10

    option.Spacing /= 10

    option.location_err /= 10

    option.LinearDiscreteDiameters = list(np.array(option.LinearDiscreteDiameters)/10)
    option.WeightedDiscreteDiameters = list(np.array(option.WeightedDiscreteDiameters)/10)

    option.ArraySize = int(option.ArraySize/10)

    mu.Ndiameter /= 10
    sigma.Ndiameter /= 10

    mu.Hdiameter /= 10
    sigma.Hdiameter /= 10

    mu.osteonlength /= 10
    sigma.osteonlength /= 10
    option.maxosteonlength /= 10
   
    return option, mu, sigma

#%%

def revertParameters(option, mu, sigma):
    '''Converts inputs back from pixel to µm units'''
    
    option.mindiameter *= 10

    option.Spacing *= 10

    option.location_err *= 10

    option.LinearDiscreteDiameters = list(np.array(option.LinearDiscreteDiameters)*10)
    option.WeightedDiscreteDiameters = list(np.array(option.WeightedDiscreteDiameters)*10)

    option.ArraySize = int(option.ArraySize*10)

    mu.Ndiameter *= 10
    sigma.Ndiameter *= 10

    mu.Hdiameter *= 10
    sigma.Hdiameter *= 10

    mu.osteonlength *= 10
    sigma.osteonlength *= 10
    option.maxosteonlength *= 10
   
    return option, mu, sigma
#%%

def nameFig(option):
    '''sets name for the file generated'''
    
    # File name:
    fpath = './pore_files/'
    if option.namestyle == 'Timestamp':
        clock_ = datetime.datetime.now()
        fname = clock_.strftime("%Y_%m_%d_%H_%M_%S") + '.tif'
    else:
        banlist = ['<','>',':','\"','/','\\','|','?','*']
        option.namestyle = option.namestyle.translate({ord(i): None for i in banlist})
        if '.tif' in option.namestyle:
            fname = option.namestyle
        else:
            fname = option.namestyle + '.tif'
            
    return fpath,fname
#%%

def setRNG(option):
    '''sets RNG seed value for generation'''
    
    if option.rng_method == 1:
        if os.path.isfile('./saved_rng.pkl'):
            with open('saved_rng.pkl', 'rb') as f:
                RNGkey = pickle.load(f)
            random.seed(RNGkey)
    elif option.rng_method == 0:
        RNGkey = random.randrange(sys.maxsize)
        random.seed(RNGkey)
        with open('saved_rng.pkl', 'wb') as f:
            pickle.dump(RNGkey, f)
    else:
        rngseed = []
        for character in option.rng_method:
            rngseed.append(ord(character))
        RNGkey = sum(rngseed)
        random.seed(RNGkey)
        with open('saved_rng.pkl', 'wb') as f:
            pickle.dump(RNGkey, f)
            
    return RNGkey 

#%%

def getPD(mu,sigma,weighting,option):
    PD = Struct()
    #getPD(mu,sigma,weighting,mindiameter) Outputs Probablility Distributions
    #
    #   Uses PD and option structs to choose method of generating radius and
    #   circularity
    #
    #   mu, sigma, and PD structs contain:
    #       Ncircularity, Ndiameter, Hcircularity, Hdiameter, SED, TOTdiameter, 
    #       TOTcircularity
    
    def span_dist(values, probs, rand_range):
        '''
        Creates a distribution broken into weighted self-uniform spans
        
        values defines the bounds between each span
        probs defines the weightings of each span
        rand_range is a two-length vector that defines the 'random' generation limits
        '''
        if values == 'rand':
            dist = stats.uniform(loc=rand_range[0], scale=rand_range[1])
        elif len(values) == 1:
            dist = stats.uniform(loc=values, scale = 0)
        else:
            dist_vector = len(probs)*[0]
        
            for n in range(len(dist_vector)):
                dist_vector[n] = stats.uniform(loc=values[n], scale=(values[n+1]-values[n]))
            dist = MixtureModel(dist_vector, weights=probs)
    
        return dist
    
    # Probability distributions for linked diameters and circularities
    # Normal SED (small, regular pores)
    PD.Ncircularity = stats.norm(loc=mu.Ncircularity, scale=sigma.Ncircularity)
    PD.Ndiameter = stats.truncnorm(a=((option.mindiameter-mu.Ndiameter)/sigma.Ndiameter) ,b=np.Inf,loc=mu.Ndiameter, scale=sigma.Ndiameter)
    # High SED (larger, irregular pores)
    PD.Hcircularity = stats.norm(loc=mu.Hcircularity, scale=sigma.Hcircularity)
    PD.Hdiameter = stats.truncnorm(a=((option.mindiameter-mu.Hdiameter)/sigma.Hdiameter), b=np.Inf,loc=mu.Hdiameter, scale=sigma.Hdiameter)
    # SED distribution
    PD.SED = stats.norm(loc=mu.SED, scale=sigma.SED)

    # Calculates weighting split for unlinked distributions
    split = lambda num: [num, 1-num]
    SED_weighting = split(stats.norm(mu.SED, sigma.SED).cdf(option.SED_limit))

    # Probability distributions for unlinked diameters and circularities    
    PD.TOTdiameter = MixtureModel([PD.Ndiameter, PD.Hdiameter], SED_weighting)
    PD.TOTcircularity = MixtureModel([PD.Ncircularity, PD.Hcircularity],SED_weighting)

    # Probability distributions for number/ length of pores
    # Minimum porosity clipped at 0.01
    PD.porosity = stats.truncnorm(a=((0.01-mu.porosity)/sigma.porosity) ,b=np.Inf,loc=mu.porosity, scale=sigma.porosity)
    # Maximum clipped at maxosteonlength
    PD.osteonlength = stats.truncnorm(a=(-mu.osteonlength/sigma.osteonlength) ,b=((option.maxosteonlength-mu.osteonlength)/sigma.osteonlength),loc=mu.osteonlength, scale=sigma.osteonlength)

    # Probability distributions for azimuthal angle
    PD.phi = span_dist(weighting.phi_values, weighting.phi_probs, [0, 0.5*np.pi])
    # Probability distributions for radial angle
    PD.theta = span_dist(weighting.theta_values, weighting.theta_probs, [0, 2*np.pi])

    return PD


#%%

def getRC(option, PD):
    '''
    Outputs radius and circularity values
    Uses PD and option structs to choose method of generating radius and circularity
    '''
    
    if option.varLink == 1:
        # linked radius and circularity
        if PD.SED.rvs() <= option.SED_limit:
            C = PD.Ncircularity.rvs()
            R = 0.5* PD.Ndiameter.rvs()
        else:
            C = PD.Hcircularity.rvs()
            R = 0.5* PD.Hdiameter.rvs()
    else:
        C = PD.TOTcircularity.rvs()
        R = 0.5* PD.TOTdiameter.rvs()
    
    if option.LinearDiscreteDiameters != []:
        R = 0.5* random.choice(option.LinearDiscreteDiameters)
    elif option.WeightedDiscreteDiameters != []:
        R = min(option.WeightedDiscreteDiameters, key=lambda x: abs(np.subtract(x,2*R)))
    
    if option.LinearDiscreteCircularities != []:
        C = random.choice(option.LinearDiscreteCircularities)
    elif option.WeightedDiscreteCircularities != []:
        C = min(option.WeightedDiscreteCircularities, key=lambda x: abs(np.subtract(x,C)))
    
    return R,C
#%%

def networkPore(valueslog, minz, z, maxz, iteration):
    '''
    Uses the log of previously generated pores to create new pores that branch off of them.
       DOI: 10.1007%2Fs11999-009-0806-x
    [R; C; theta; phi; x; y; minz; z; maxz; isfilled; A; B]
    '''
    
    # minz and maxz expressed as offset from z
    maxz = maxz - z
    minz = minz - z
    # picks a random pore
    select = random.randint(0,sum(valueslog[9,:])-1)
    # pool of 'spaces' on top and bottom of pore
    Mpool = [0]*int(valueslog[10,select]) + [1]*int(valueslog[11,select])
    # picks from the pool
    M = random.choice(Mpool)
    # indexing shift for valueslog
    N = 6 + 2*M
    # randomly choose branching direction
    # 0.75 chance for branching and side directions to align
    D = random.randint(0,3)
    
    # set pore coordinates to end of chosen pore
    x = valueslog[4,select] + (-1)**(M) *np.sin(valueslog[2,select])*np.tan(valueslog[3,select]) \
    *abs(valueslog[7,select]-valueslog[N,select])
    
    y = valueslog[5,select] + (-1)**(M+1) *np.cos(valueslog[2,select])*np.tan(valueslog[3,select]) \
    *abs(valueslog[7,select]-valueslog[N,select])
    
    z = valueslog[N,select]
    
    # reduces 'pool' according to choice
    valueslog[10+M,select] -= 1
    
    # actually sets branching direction
    if M == 0:
        if D == 0:
            minz = valueslog[7,select]
            maxz = maxz + z
        else:
            maxz = valueslog[7,select]
            minz = minz + z
    else:
        if D == 0:
            maxz = valueslog[7,select]
            minz = minz + z
        else:
            minz = valueslog[7,select]
            maxz = maxz + z
    
    # removes a pore from branching consideration after Mpool is empty
    if valueslog[10,select]+valueslog[11,select] == 0:
        valueslog = np.delete(valueslog,select,1)
        iteration -= 1
    
    return x,y,minz,z,maxz,valueslog,iteration
#%%

def getXY(option, XYprimer):
    '''
    Outputs location values and data to be used on loop
    Uses XYprimer and option structs to choose method of generating location
    '''
    
    if option.LocationType not in XYprimer.Locations:
        XYprimer.grid_complete = 1;
        # random pore location
        x = option.ArraySize*random.random()
        y = option.ArraySize*random.random()

    elif XYprimer.Locations[option.LocationType] in [1, 2]:
        
        Circle = XYprimer.Locations[option.LocationType] == 1
        Square = XYprimer.Locations[option.LocationType] == 2
        
        XYprimer.iu += 1
        if (XYprimer.iu == (len(XYprimer.AngleList))*Circle 
            + (len(XYprimer.SpaceList) + 1 - option.ignoreborder)*Square ):
            if XYprimer.it == len(XYprimer.SpaceList) - option.ignoreborder*Square:
                XYprimer.it = option.ignoreborder*Square
                XYprimer.grid_complete = 1
                if option.ignore_target_porosity:
                    XYprimer.ignore_target_porosity = 1
            XYprimer.it += 1
            XYprimer.iu = 1 + option.ignoreborder*Square
            if Circle:
                XYprimer.AngleList = np.linspace(0, 2*np.pi, XYprimer.it*int(np.sqrt(option.ArraySize/option.Spacing)))
        if Circle:
            # circular grid for pore location
            x = (0.5*option.ArraySize + XYprimer.SpaceList[XYprimer.it - 1]*np.cos(XYprimer.AngleList[XYprimer.iu - 1])) + option.location_err*(2*random.random() - 1)
            y = (0.5*option.ArraySize + XYprimer.SpaceList[XYprimer.it - 1]*np.sin(XYprimer.AngleList[XYprimer.iu - 1])) + option.location_err*(2*random.random() - 1);
        elif Square:
            # square grid for pore location
            x = XYprimer.SpaceList[XYprimer.it - 1] + option.location_err*(2*random.random() - 1)
            y = XYprimer.SpaceList[XYprimer.iu - 1] + option.location_err*(2*random.random() - 1)


    return x, y, XYprimer
#%%              
              
def poreBlast(Bone):
    '''deposits bone matrix -> decrease porosity (more 1's)'''
        
    threshold = 15
    
    b = 1/threshold
    kernel = np.ones((3,3,3))*b
    kernel[1,1,1] = 1
    
    conv = scipy.ndimage.convolve(Bone, kernel)
    rounded = np.floor(conv)
    
    rand = np.random.randint(0, 2, Bone.shape)
    salted = np.multiply(rand,rounded).astype(dtype=np.float32)
    mixed = Bone + salted    
    
    flattened = mixed.astype(bool).astype(dtype=np.float32)
    
    return flattened

#%%

def poreClast(Bone):
    '''removes bone matrix -> increase porosity (more 0's)'''
    
    threshold = 12
    
    b = 1/(3**3)  # borders
    m = 1 - b*(threshold+1)  # middle
    kernel = np.ones((3,3,3))*b
    kernel[1,1,1] = m
    
    conv = scipy.ndimage.convolve(Bone, kernel, mode='constant', cval=0)
    rounded = np.floor(conv)
    
    return rounded

#%%

def getTextOutput(option, mu, sigma, weighting, params, target_porosity, porosity, RNGkey, fname):
    
    varlist = [
    'option.varLink' ,
    'option.mindiameter' ,
    'option.LocationType' ,
    'option.Spacing' ,
    'option.location_err' ,
    'option.ignore_target_porosity',
    'option.ignoreborder' ,
    'option.LinearDiscreteDiameters' ,
    'option.WeightedDiscreteDiameters' ,
    'option.LinearDiscreteCircularities' ,
    'option.WeightedDiscreteCircularities' ,
    'option.smoothPores' ,
    'option.variedPoreShape' ,
    'option.ArraySize' ,
    'option.maxosteonlength' ,
    'option.SED_limit',
    'option.TP_CORRECTION_FACTOR',
    'option.experimental_porosity',
    'mu.SED' ,
    'mu.Ndiameter' ,
    'mu.Ncircularity' ,
    'mu.Hdiameter' ,
    'mu.Hcircularity' ,
    'mu.osteonlength' ,
    'mu.porosity',
    'sigma.SED' ,
    'sigma.Ndiameter' ,
    'sigma.Ncircularity' ,
    'sigma.Hdiameter' ,
    'sigma.Hcircularity' ,
    'sigma.osteonlength' ,
    'sigma.porosity',
    'weighting.phi_values' ,
    'weighting.phi_probs' ,
    'weighting.theta_values' ,
    'weighting.theta_probs' ,
    'target_porosity' ,
    'params.pores_before_networking' ,
    'params.top_branches' ,
    'params.bottom_branches' ,
    'params.sealed_osteon_chance' ,
    'params.transverse_flag_onset' ,
    'params.shape_proportions' ,
    'RNGkey' ,
    'porosity'
    ]
    
    fullcell = list()
    for n in range(0, len(varlist)):
        fullcell.append(varlist[n]+' is <<'+str(eval(varlist[n]))+'>>'+'\n')
    
    sheetprep = list()
    sheetprep.append('Filename')
    sheetprep.append('')
    for n in range(2,len(varlist)+2):
        sheetprep.append(varlist[n-2])
        
    sheetcell = np.zeros([1,len(varlist)+2], dtype=object)
    sheetcell[0,0] = fname
    for n in range(2,len(varlist)+2):
        sheetcell[0,n] = str(eval(varlist[n-2]))
         
    
    return fullcell, sheetprep, sheetcell
#%%

def make3DModel(fpath, fname, Bone):
    '''
    Converts array to STL
    Credit to Aleksandr Burakov, https://stackoverflow.com/questions/69524209/how-to-generate-3d-mesh-from-a-numpy-binary-mask
    '''
    # Mesh elements per voxel side length
    res = 1
    
    Bone32 = np.float32(~Bone)
    simpleVolume = mrn.simpleVolumeFrom3Darray(Bone32)
    floatGrid = mr.simpleVolumeToDenseGrid(simpleVolume)
    mesh = mr.gridToMesh(floatGrid, mr.Vector3f(1/res, 1/res, 1/res), 0.5)
    mr.saveMesh(mesh, fpath+fname.removesuffix('.tif')+'.stl')
    