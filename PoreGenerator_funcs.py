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
import fnmatch

from PoreGenerator_classes import MixtureModel
class Struct:
    pass

#%%

def chooseExports(exports):
    '''packages "exports" input into Struct'''
    class Struct:
        pass
    export = Struct()

    dictionary = {
        'xlsx' : 1,
        'XLSX' : 1,
        'excel' : 1,
        'spreadsheet' : 1,
        'txt' : 2,
        'TXT' : 2,
        'text' : 2,
        'tiff' : 3,
        'tif' : 3,
        'TIFF' : 3,
        'TIF' : 3,
        'stl' : 4,
        'STL' : 4
        }
    
    export.xlsx = False
    export.txt = False
    export.tiff = False
    export.stl = False
    
    lst = []
    
    for exp in exports:
        if exp in dictionary:
            lst.append(dictionary[exp])
        else:
            raise(ValueError('"exports" input includes invalid type.'))

    if 1 in lst:
        export.xlsx = True
    if 2 in lst:
        export.txt = True
    if 3 in lst:
        export.tiff = True
    if 4 in lst:
        export.stl = True
            
    if (export.xlsx + export.txt + export.tiff + export.stl) == 0:
        raise(Exception('No file outputs selected.'))
        
    return export

# %%

def scaleParameters(option, mu, sigma, scale):
    '''
    if scale is 1/10: Converts inputs from µm to voxel units
    
    if scale is 10: Converts inputs back from voxel to µm units
    '''
    
    option.mindiameter *= scale

    option.Spacing *= scale

    option.location_err *= scale

    option.LinearDiscreteDiameters = list(np.array(option.LinearDiscreteDiameters)*scale)
    option.WeightedDiscreteDiameters = list(np.array(option.WeightedDiscreteDiameters)*scale)

    option.ArraySize = int(option.ArraySize*scale)

    mu.Ndiameter *= scale
    sigma.Ndiameter *= scale

    mu.Hdiameter *= scale
    sigma.Hdiameter *= scale

    mu.osteonlength *= scale
    sigma.osteonlength *= scale
    option.maxosteonlength *= scale
   
    return option, mu, sigma

#%%

def nameFig(namestyle):
    '''sets name for the file generated'''
    
    # File name:
    fpath = './pore_files/'
    clock_ = datetime.datetime.now()
    if namestyle == 'Timestamp':
        fname = clock_.strftime("%Y_%m_%d_%H_%M_%S") + '.tif'
    else:
        banlist = ['<','>',':','\"','/','\\','|','?','*']
        namestyle = namestyle.translate({ord(i): None for i in banlist})
        # You can prefix the Timestamp name style
        if len(fnmatch.filter([namestyle],"prefix(*)")) > 0:
            fname = namestyle.removeprefix("prefix(").removesuffix(")") + clock_.strftime("%Y_%m_%d_%H_%M_%S") + '.tif'
        elif namestyle.endswith('.tif') or namestyle.endswith(".tiff"):
            fname = namestyle
        else:
            fname = namestyle + '.tif'
            
    return fpath,fname
#%%

def setRNG(rng_method):
    '''sets RNG seed value for generation'''
    
    if (rng_method == 1) and os.path.isfile('./saved_rng.pkl'):
        with open('saved_rng.pkl', 'rb') as f:
            RNGkey = pickle.load(f)
        random.seed(RNGkey)
    elif (rng_method == 0) or (not os.path.isfile('./saved_rng.pkl')):
        RNGkey = random.randrange(sys.maxsize)
        random.seed(RNGkey)
        with open('saved_rng.pkl', 'wb') as f:
            pickle.dump(RNGkey, f)
    else:
        rngseed = []
        for character in rng_method:
            rngseed.append(ord(character))
        RNGkey = sum(rngseed)
        random.seed(RNGkey)
        with open('saved_rng.pkl', 'wb') as f:
            pickle.dump(RNGkey, f)
            
    return RNGkey 

#%%

def getPD(mu, sigma, weighting, option):
    PD = Struct()
    '''
    Creates probability distributions from passed parameters
    '''
    
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
            if sum(probs) != 1:
                raise(ValueError('weighting.phi_probs and weighting.theta_probs must each sum to 1.'))
            dist_vector = len(probs)*[0]
        
            for n in range(len(dist_vector)):
                dist_vector[n] = stats.uniform(loc=values[n], scale=(values[n+1]-values[n]))
            dist = MixtureModel(dist_vector, weights=probs)
    
        return dist
    
    # Probability distributions for linked diameters and circularities
    # Normal SED (small, regular pores)
    PD.Ncircularity = stats.norm(loc=mu.Ncircularity, scale=sigma.Ncircularity)
    PD.Ndiameter = stats.truncnorm(a=((option.mindiameter-mu.Ndiameter)/sigma.Ndiameter) ,b=np.inf,loc=mu.Ndiameter, scale=sigma.Ndiameter)
    # High SED (larger, irregular pores)
    PD.Hcircularity = stats.norm(loc=mu.Hcircularity, scale=sigma.Hcircularity)
    PD.Hdiameter = stats.truncnorm(a=((option.mindiameter-mu.Hdiameter)/sigma.Hdiameter), b=np.inf,loc=mu.Hdiameter, scale=sigma.Hdiameter)
    # SED distribution
    PD.SED = stats.norm()

    # Calculates weighting split for unlinked distributions
    split = lambda num: [num, 1-num]
    SED_weighting = split(stats.norm().cdf(option.SED_limit))

    # Probability distributions for unlinked diameters and circularities    
    PD.TOTdiameter = MixtureModel([PD.Ndiameter, PD.Hdiameter], SED_weighting)
    PD.TOTcircularity = MixtureModel([PD.Ncircularity, PD.Hcircularity],SED_weighting)

    # Probability distributions for number/ length of pores
    # Minimum porosity clipped at 0.01
    PD.porosity = stats.truncnorm(a=((0.01-mu.porosity)/sigma.porosity) ,b=np.inf,loc=mu.porosity, scale=sigma.porosity)
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
    numpores = sum(valueslog[9,:])
    numpores = numpores.astype(int)
    select = random.randint(0,numpores-1)
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
                else:
                    XYprimer.iter += 1
                    if XYprimer.iter == XYprimer.max_iters:
                        XYprimer.ignore_target_porosity = 1
                        print(f'Grid populated {XYprimer.max_iters} times without meeting target porosity. Aborting further population attempts.')
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
    else:
        raise(ValueError('XYprimer.Locations dictionary contains references to non-0 or -1 values.'))


    return x, y, XYprimer
#%%              
              
def poreBlast(Bone):
    '''deposits bone matrix -> decrease porosity (more 1's)'''
        
    # Create convolution kernel
    threshold = 15
    
    b = 1/threshold
    kernel = np.ones((3,3,3))*b
    kernel[1,1,1] = 1
    
    # Convolve Bone with kernel
    conv = scipy.ndimage.convolve(Bone, kernel)
    # Floor function to revert to binary
    rounded = np.floor(conv)
    
    # Only change some voxels this way
    # Because poreBlast is greedier than poreClast
    rand = np.random.randint(0, 2, Bone.shape)
    salted = np.multiply(rand,rounded).astype(dtype=np.float32)
    mixed = Bone + salted    
    
    # Revert back to binary values
    flattened = mixed.astype(bool).astype(dtype=np.float32)
    
    return flattened

#%%

def poreClast(Bone):
    '''removes bone matrix -> increase porosity (more 0's)'''
    
    # Create convolution kernel
    threshold = 12
    
    b = 1/(3**3)  # borders
    m = 1 - b*(threshold+1)  # middle
    kernel = np.ones((3,3,3))*b
    kernel[1,1,1] = m
    
    # Convolve Bone with kernel
    conv = scipy.ndimage.convolve(Bone, kernel, mode='constant', cval=0)
    # Floor function to revert to binary
    rounded = np.floor(conv)
    
    return rounded

#%%

def getTextOutput(option, mu, sigma, weighting, params, target_porosity, porosity, RNGkey, fname):
    '''
    Generates list of all parameter values for export
    '''
    
    # List of all parameters to be exported
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
    'mu.Ndiameter' ,
    'mu.Ncircularity' ,
    'mu.Hdiameter' ,
    'mu.Hcircularity' ,
    'mu.osteonlength' ,
    'mu.porosity',
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
    
    # creates list for text file export
    fullcell = list()
    for n in range(0, len(varlist)):
        fullcell.append(varlist[n]+' is <<'+str(eval(varlist[n]))+'>>'+'\n')
    
    # creates header for excel file export
    sheetprep = list()
    sheetprep.append('Filename')
    sheetprep.append('')
    for n in range(0,len(varlist)):
        sheetprep.append(varlist[n])
    
    # creates new row for excel file
    sheetcell = np.zeros([1,len(varlist)+2], dtype=object)
    sheetcell[0,0] = fname
    for n in range(0,len(varlist)):
        sheetcell[0,n+2] = str(eval(varlist[n]))
         
    
    return fullcell, sheetprep, sheetcell
#%%

def make3DModel(fpath, fname, Bone):
    '''
    Converts array to STL
    Credit to Aleksandr Burakov, https://stackoverflow.com/questions/69524209/how-to-generate-3d-mesh-from-a-numpy-binary-mask
    '''
    # Mesh elements per voxel side length
    res = 1
    
    Bone32 = np.float32(~Bone.astype(bool))
    simpleVolume = mrn.simpleVolumeFrom3Darray(Bone32)
    floatGrid = mr.simpleVolumeToDenseGrid(simpleVolume)
    mesh = mr.gridToMesh(floatGrid, mr.Vector3f(1/res, 1/res, 1/res), 0.5)
    mr.saveMesh(mesh, fpath+fname.removesuffix('.tif')+'.stl')
    