# -*- coding: utf-8 -*-
# LoadParameters.py
"""
Created on Tue Jun 18 12:16:52 2024

@author: ztoth
"""
import pandas as pd
from numpy import pi, inf
class Struct:
    pass

def LoadParameters(param_file, exports):
    
    option = Struct()
    mu = Struct()
    sigma = Struct()
    weighting = Struct()
    params = Struct()
    export = Struct()
    
    param_file = './param_files/' + param_file
    
    export = choose_exports(exports)
    
    if param_file is None:
        
        '''
        Options
        '''
        # True if variation in diameter and circularity should be linked
        # eg, large pored would tend to also be oblong
        option.varLink = True
        # Minimum pore diameter permissible (µm)
        option.mindiameter = 0

        # 'Circle' (or 1) for circular grid, 'Square' (or 2) for square grid, else random
        option.LocationType = 0
        # If LocationType is not random, distance between grid lines (µm)
        option.Spacing = 160  # 160
        # If LocationType is not random, variation of pores from grid lines (µm)
        option.location_err = 0  # 0
        # If LocationType is not random, stops generation after grid is complete
        # Can cause problems if False while location_err is low
        option.ignore_target_porosity = False
        # If LocationType is square grid, whether to generate pores along the edge.
        option.ignoreborder = False
        
        # List of diameters to randomly choose from (µm)
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
        option.smoothPores = True
        
        # True if pores should have varying shape (cylinder,cone,ellipsoid,hyperboloid)
        option.variedPoreShape = True
        
        # Array size (µm): 10 micrometers/voxel
        option.ArraySize = 2000  # 2000
        target_porosity = 0
        
        # Choose value from experimental distribution instead of using target_porosity
        option.experimental_porosity = True
        
        # readjusts target_porosity to account for loss when smoothPores is used
        option.TP_CORRECTION_FACTOR = 1
        
        '''
        Parameters
        '''
        # SED related parameters below from DOI: 10.1002/jbmr.3561
        # Parameters for SED per hole
        # Controls the distribution of pore size/irregularity
        mu.SED = 8  # 8
        sigma.SED = 4.5  # 4.5
        option.SED_limit = 17  # 17
        
        # Parameters for Normal SED
        # diameter (µm)
        mu.Ndiameter = 48  # 48
        sigma.Ndiameter = 4  # 4
        # circularity
        mu.Ncircularity = 0.66  # 0.66
        sigma.Ncircularity = 0.03  # 0.03
        
        # Parameters for High SED
        # diameter (µm)
        mu.Hdiameter = 157  # 157
        sigma.Hdiameter = 37.5  # 37.5
        # circularity
        mu.Hcircularity = 0.42  # 0.42
        sigma.Hcircularity = 0.045  # 0.045
        
        # Parameters for Canal Length (µm)
        mu.osteonlength = 1000  # 1000
        sigma.osteonlength = 375  # 375
        option.maxosteonlength = 2200/3 # 2200/3
        
        # Parameters for Porosity
        mu.porosity = 0.075046512  # 0.075046512
        sigma.porosity = 0.036744908  # 0.036744908
        
        # Parameters for phi value selection DOI: 10.1016/8756-3282(94)90288-7
        weighting.phi_values = [0,pi/12,pi/2]  # [0,pi/12,pi/2]
        weighting.phi_probs = [0.5, 0.5]  # [0.5, 0.5]
        
        # Parameters for theta value selection
        weighting.theta_values = [0, 2*pi]  # [0,2*pi]
        weighting.theta_probs = [1]  # [1]
        
        # Proportions of each pore shape: [Cylinder, Proximal Cone, Distal Cone,
        #                                                   Ellipsoid, Hyperboloid]
        params.shape_proportions = [0.392, 0.094, 0.351, 0.122, 0.041]  # [0.392, 0.094, 0.351, 0.122, 0.041] DOI: 10.1111/j.1439-0264.2009.00973.x
        
        params.pores_before_networking = 75  # 75
        params.top_branches = [0,2] # [0,2]
        params.bottom_branches = [0,2] # [0,2]
        params.sealed_osteon_chance = 0.068  # 0.068 # DOI: 10.1002/ar.21309
        params.transverse_flag_onset = pi/4  # pi/4

    else:
        
        def clean(string):
            out = string.replace(' is ', '').replace('>>', '').split('<<')
            return out

        df = pd.read_table(param_file, header=None, dtype='string')
        lst = list()
        for n in range(0, len(df[0])): lst.append(df[0].loc[n])
        df = pd.DataFrame(list(map(clean, lst)), columns=['Variable', 'Value'])

        option.varLink = eval(df.iloc[0].Value)
        option.mindiameter = eval(df.iloc[1].Value)
        option.LocationType = eval(df.iloc[2].Value)
        option.Spacing = eval(df.iloc[3].Value)
        option.location_err = eval(df.iloc[4].Value)
        option.ignore_target_porosity = eval(df.iloc[5].Value)
        option.ignoreborder = eval(df.iloc[6].Value)
        option.LinearDiscreteDiameters = eval(df.iloc[7].Value)
        option.WeightedDiscreteDiameters = eval(df.iloc[8].Value)
        option.LinearDiscreteCircularities = eval(df.iloc[9].Value)
        option.WeightedDiscreteCircularities = eval(df.iloc[10].Value)
        option.smoothPores = eval(df.iloc[11].Value)
        option.variedPoreShape = eval(df.iloc[12].Value)
        option.ArraySize = eval(df.iloc[13].Value)
        option.maxosteonlength = eval(df.iloc[14].Value)
        option.SED_limit = eval(df.iloc[15].Value)
        option.TP_CORRECTION_FACTOR = eval(df.iloc[16].Value)
        option.experimental_porosity = eval(df.iloc[17].Value)
        mu.SED = eval(df.iloc[18].Value)
        mu.Ndiameter = eval(df.iloc[19].Value)
        mu.Ncircularity = eval(df.iloc[20].Value)
        mu.Hdiameter = eval(df.iloc[21].Value)
        mu.Hcircularity = eval(df.iloc[22].Value)
        mu.osteonlength = eval(df.iloc[23].Value)
        mu.porosity = eval(df.iloc[24].Value)
        sigma.SED = eval(df.iloc[25].Value)
        sigma.Ndiameter = eval(df.iloc[26].Value)
        sigma.Ncircularity = eval(df.iloc[27].Value)
        sigma.Hdiameter = eval(df.iloc[28].Value)
        sigma.Hcircularity = eval(df.iloc[29].Value)
        sigma.osteonlength = eval(df.iloc[30].Value)
        sigma.porosity = eval(df.iloc[31].Value)
        weighting.phi_values = eval(df.iloc[32].Value)
        weighting.phi_probs = eval(df.iloc[33].Value)
        weighting.theta_values = eval(df.iloc[34].Value)
        weighting.theta_probs = eval(df.iloc[35].Value)
        target_porosity = eval(df.iloc[36].Value)
        params.pores_before_networking = eval(df.iloc[37].Value)
        params.top_branches = eval(df.iloc[38].Value)
        params.bottom_branches = eval(df.iloc[39].Value)
        params.sealed_osteon_chance = eval(df.iloc[40].Value)
        params.transverse_flag_onset = eval(df.iloc[41].Value)
        params.shape_proportions = eval(df.iloc[42].Value)

    
    return option, target_porosity, export, mu, sigma, weighting, params

def choose_exports(exports):
    
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
    
    if type(exports) is str:
        exports = [exports]
    
    lst = []
    for exp in exports: lst.append(dictionary[exp])
    
    if 1 in lst:
        export.xlsx = True
    if 2 in lst:
        export.txt = True
    if 3 in lst:
        export.tiff = True
    if 4 in lst:
        export.stl = True
        
    return export