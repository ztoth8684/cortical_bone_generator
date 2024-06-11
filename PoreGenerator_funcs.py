# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:09:28 2023

@author: ztoth
"""

import datetime
import numpy as np
import scipy.stats as stats
import random
import meshlib.mrmeshpy as mr
import meshlib.mrmeshnumpy as mrn

from MixtureModel import MixtureModel
class Struct:
    pass


#%%

def nameFig(option):
    # sets name for the file generated
    # File name:
    # fpath = 'C:\\Users\\ztoth\\Documents\\GitHub\\cortical_bone_generator\\pore_files\\'
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
    
    # Probability distributions for linked diameters and circularities
    
    PD.Ncircularity = stats.norm(loc=mu.Ncircularity, scale=sigma.Ncircularity)
    PD.Ndiameter = stats.truncnorm(a=((option.mindiameter-mu.Ndiameter)/sigma.Ndiameter) ,b=np.Inf,loc=mu.Ndiameter, scale=sigma.Ndiameter)
    
    PD.Hcircularity = stats.norm(loc=mu.Hcircularity, scale=sigma.Hcircularity)
    PD.Hdiameter = stats.truncnorm(a=((option.mindiameter-mu.Hdiameter)/sigma.Hdiameter) ,b=np.Inf,loc=mu.Hdiameter, scale=sigma.Hdiameter)
    
    PD.SED = stats.norm(loc=mu.SED, scale=sigma.SED)
    
    # Probability distributions for unlinked diameters and circularities    
    PD.TOTdiameter = MixtureModel([PD.Ndiameter, PD.Hdiameter],weighting.SED)
    PD.TOTcircularity = MixtureModel([PD.Ncircularity, PD.Hcircularity],weights= weighting.SED)
    
    PD.porosity = stats.truncnorm(a=((0.01-mu.porosity)/sigma.porosity) ,b=np.Inf,loc=mu.porosity, scale=sigma.porosity)
    PD.osteonlength = stats.truncnorm(a=(-mu.porosity/sigma.porosity) ,b=((option.maxosteonlength-mu.osteonlength)/sigma.osteonlength),loc=mu.osteonlength, scale=sigma.osteonlength)
    
    PD.phi1 = stats.uniform(loc=weighting.phi_values[0],scale=weighting.phi_values[1]-weighting.phi_values[0])
    PD.phi2 = stats.uniform(loc=weighting.phi_values[1],scale=weighting.phi_values[2]-weighting.phi_values[1])
    PD.phi = MixtureModel([PD.phi1, PD.phi2],weights= weighting.phi_probs)
    
    return PD
#%%

def getRC(option,PD):
    #getRC(PD,option) Outputs radius and circularity values
    #
    #   Uses PD and option structs to choose method of generating radius and
    #   circularity
    #
    #   PD struct contains:
    #       Ncircularity, Ndiameter, Hcircularity, Hdiameter, SED, TOTdiameter, 
    #       TOTcircularity
    #   option struct contains:
    #       varLink, LinearDiscreteDiameters, WeightedDiscreteDiameters, 
    #       LinearDiscreteCircularities, WeightedDiscreteCircularities
    
    if option.varLink == 1:
        # linked radius and circularity
        E = PD.SED.rvs()
        if E <= 17:
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

def networkPore(valueslog,minz,z,maxz,iteration):
    # Uses the log of previously generated pores to create new pores that branch off of them.
    #   DOI: 10.1007%2Fs11999-009-0806-x
    #[R; C; theta; phi; x; y; minz; z; maxz; isfilled; A; B]
    
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

def getXY(option,XYprimer):
    #getXY(option,XYprimer) Outputs location values and data to be used on loop
    #
    #   Uses XYprimer and option structs to choose method of generating
    #   location
    #
    #   XYprimer struct contains:
    #       iu, it, SpaceList, AngleList, flag, break
    #   option struct contains:
    #       LocationType, location_err, ArraySize, Spacing, ignoreborder
    if (option.LocationType in [1,'Circle','Radial']):
        XYprimer.iu += 1
        if XYprimer.iu == len(XYprimer.AngleList):
            if XYprimer.it == len(XYprimer.SpaceList):
                XYprimer.it = 0
                XYprimer.griddone = 1
                if option.location_err < 6:
                    XYprimer.ignoreTP = 1
            XYprimer.it += 1
            XYprimer.iu = 1
            XYprimer.AngleList = np.linspace(0,2*np.pi,XYprimer.it*int(np.sqrt(option.ArraySize/option.Spacing)))
        
        # circular grid for pore location
        x = (0.5*option.ArraySize + XYprimer.SpaceList[XYprimer.it -1]*np.cos(XYprimer.AngleList[XYprimer.iu -1])) + option.location_err*(2*random.random() -1)
        y = (0.5*option.ArraySize + XYprimer.SpaceList[XYprimer.it -1]*np.sin(XYprimer.AngleList[XYprimer.iu -1])) + option.location_err*(2*random.random() -1);

    elif (option.LocationType in [2,'Square']):
        XYprimer.iu += 1
        if XYprimer.iu == len(XYprimer.SpaceList) + 1 - option.ignoreborder:
            if XYprimer.it == len(XYprimer.SpaceList) - option.ignoreborder:
                XYprimer.it = option.ignoreborder
                XYprimer.griddone = 1
                if option.location_err < 6:
                    XYprimer.ignoreTP = 1
            XYprimer.it += 1
            XYprimer.iu = 1 + option.ignoreborder
            
        # square grid for pore location
        x = XYprimer.SpaceList[XYprimer.it -1] + option.location_err*(2*random.random() -1)
        y = XYprimer.SpaceList[XYprimer.iu -1] + option.location_err*(2*random.random() -1)
        
    else:
        XYprimer.griddone = 1;
        # random pore location
        x = option.ArraySize*random.random()
        y = option.ArraySize*random.random()



    return x,y,XYprimer
#%%              
              
def erodePores(Bone):
    #makes area adjacent to high amount of pores more porous
        
    adjacency = [(i,j,k) for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1) if not (i == j == k == 0)] #the adjacency matrix

    BoneMerger = np.zeros(Bone.shape, dtype=object)
    BoneOrig = Bone

    for ia in range(2,Bone.shape[0]-1):
        for ib in range(2,Bone.shape[1]-1):
            for ic in range(2,Bone.shape[2]-1):
                if BoneOrig[ia,ib,ic] == 0:
                    BoneMerger[ia,ib,ic] = [BoneOrig[ia+dx, ib+dy, ic+dz] for dx, dy, dz in adjacency]
                    BoneMerger[ia,ib,ic] = sum(list(map(int,BoneMerger[ia,ib,ic])))
                if BoneMerger[ia,ib,ic] >= 15 and random.randint(0,1) == 0:
                    Bone[ia,ib,ic] = 1

    return Bone
#%%

def mergePores(Bone):
    #makes area adjacent to high amount of pores more porous
    
    adjacency = [(i,j,k) for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1) if not (i == j == k == 0)] #the adjacency matrix

    BoneMerger = np.zeros(Bone.shape, dtype=object)
    BoneOrig = Bone
    
    for ia in range(2,Bone.shape[0]-1):
        for ib in range(2,Bone.shape[1]-1):
            for ic in range(2,Bone.shape[2]-1):
                if BoneOrig[ia,ib,ic] == 1:
                    BoneMerger[ia,ib,ic] = [BoneOrig[ia+dx, ib+dy, ic+dz] for dx, dy, dz in adjacency]
                    BoneMerger[ia,ib,ic] = sum(list(map(int,BoneMerger[ia,ib,ic])))
                if BoneMerger[ia,ib,ic] <= 12:
                    Bone[ia,ib,ic] = 0

    return Bone              
#%%

def getTextOutput(option, mu, sigma, weighting, TargetPorosity, pores_before_networking, \
                  sealed_osteon_chance, transverse_flag_onset, shape_proportions, RNGkey, fname):
    
    varlist = [
    'option.varLink' ,
    'option.mindiameter' ,
    'option.LocationType' ,
    'option.Spacing' ,
    'option.location_err' ,
    'option.ignoreborder' ,
    'option.LinearDiscreteDiameters' ,
    'option.WeightedDiscreteDiameters' ,
    'option.LinearDiscreteCircularities' ,
    'option.WeightedDiscreteCircularities' ,
    'option.mergePores' ,
    'option.variedPoreShape' ,
    'option.ArraySize' ,
    'option.maxosteonlength' ,
    'mu.SED' ,
    'mu.Ndiameter' ,
    'mu.Ncircularity' ,
    'mu.Hdiameter' ,
    'mu.Hcircularity' ,
    'mu.osteonlength' ,
    'sigma.SED' ,
    'sigma.Ndiameter' ,
    'sigma.Ncircularity' ,
    'sigma.Hdiameter' ,
    'sigma.Hcircularity' ,
    'sigma.osteonlength' ,
    'weighting.SED' ,
    'weighting.phi_values' ,
    'weighting.phi_probs' ,
    'TargetPorosity' ,
    'pores_before_networking' ,
    'sealed_osteon_chance' ,
    'transverse_flag_onset' ,
    'shape_proportions' ,
    'RNGkey' ,
    ]
    
    fullcell = list()
    for n in range(0, len(varlist)):
        fullcell.append('The value of '+varlist[n]+' is <<'+str(eval(varlist[n]))+'>>'+'\n')
    
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
    
    Bone32 = np.float32(~Bone)
    simpleVolume = mrn.simpleVolumeFrom3Darray(Bone32)
    floatGrid = mr.simpleVolumeToDenseGrid(simpleVolume)
    mesh = mr.gridToMesh(floatGrid, mr.Vector3f(0.1,0.1,0.1),0.5)
    mr.saveMesh(mesh, fpath+fname.removesuffix('.tif')+'.stl')
    
    