# -*- coding: utf-8 -*-
# LoadParameters.py
"""
Created on Tue Jun 18 12:16:52 2024

@author: ztoth
"""
from numpy import pi
class Struct:
    pass

def LoadParameters(param_file = None):
    
    option = Struct()
    mu = Struct()
    sigma = Struct()
    weighting = Struct()
    params = Struct()
    export = Struct()
    
    if param_file is None:
        
        '''
        Options
        '''
        
        # Set to False to use random seed.
        # Set to True to keep the same seed used last generation. 
        # Any other value will be used as a seed.
        option.rng_method = False
        # 'Timestamp' or name to be used
        option.namestyle = 'Timestamp' # 'Timestamp'

        # True if variation in diameter and circularity should be linked
        # eg, large pored would tend to also be oblong
        option.varLink = True
        # Minimum pore diameter permissible
        option.mindiameter = 0

        # 'Circle' (or 1) for circular grid, 'Square' (or 2) for square grid, else random
        option.LocationType = 0
        # If LocationType is not random, distance between grid lines
        option.Spacing = 16  # 16
        # If LocationType is not random, variation of pores from grid lines
        # values 0-5 ignore target_porosity
        option.location_err = 6  # 6
        # If LocationType is not random, stops generation after grid is complete
        # Can cause problems if False while location_err is low
        option.ignore_target_porosity = False
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
        option.smoothPores = True
        
        # True if pores should have varying shape (cylinder,cone,ellipsoid,hyperboloid)
        option.variedPoreShape = True
        
        # 10 micrometers/voxel
        option.ArraySize = 200  # 200
        # 'Exp' to choose value from experimental distribution
        target_porosity = 'Exp'  # 'Exp'

        # Files to Export: Excel, Text, TIFF Stack, STL
        export.xcls = True
        export.txt = False
        export.tiff = True
        export.stl = False
        
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
        sigma.osteonlength = 37.5  # 37.5
        option.maxosteonlength = 220/3  # 220/3
        
        # Parameters for Porosity
        mu.porosity = 0.075046512  # 0.075046512
        sigma.porosity = 0.036744908  # 0.036744908
        
        # Parameters for phi value selection DOI: 10.1016/8756-3282(94)90288-7
        weighting.phi_values = [0, pi/12, pi/2]  # [0,pi/12,pi/2]
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
        raise Exception
    
    return option, target_porosity, export, mu, sigma, weighting, params