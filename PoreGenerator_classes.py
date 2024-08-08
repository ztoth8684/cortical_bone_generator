# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:45:18 2023

"""
import numpy as np
import scipy.stats as stats
#%%

class MixtureModel(stats.rv_continuous):
    """
    @author: Jenny Shoars, ramzeek: https://stackoverflow.com/questions/
    47759577/creating-a-mixture-of-probability-distributions-for-sampling
    """
    
    def __init__(self, submodels, *args, weights = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.submodels = submodels
        if weights is None:
            weights = [1 for _ in submodels]
        if len(weights) != len(submodels):
            raise(ValueError(f'There are {len(submodels)} submodels and {len(weights)} weights, but they must be equal.'))
        self.weights = [w / sum(weights) for w in weights]
        
    def _pdf(self, x):
        pdf = self.submodels[0].pdf(x) * self.weights[0]
        for submodel, weight in zip(self.submodels[1:], self.weights[1:]):
            pdf += submodel.pdf(x)  * weight
        return pdf
            
    def _sf(self, x):
        sf = self.submodels[0].sf(x) * self.weights[0]
        for submodel, weight in zip(self.submodels[1:], self.weights[1:]):
            sf += submodel.sf(x)  * weight
        return sf

    def _cdf(self, x):
        cdf = self.submodels[0].cdf(x) * self.weights[0]
        for submodel, weight in zip(self.submodels[1:], self.weights[1:]):
            cdf += submodel.cdf(x)  * weight
        return cdf        

    def rvs(self, size):
        submodel_choices = np.random.choice(len(self.submodels), size=size, p = self.weights)
        submodel_samples = [submodel.rvs(size=size) for submodel in self.submodels]
        rvs = np.choose(submodel_choices, submodel_samples)
        return rvs
    
# %%

class XYprimer:
    '''
    Pre-sets some variables for choosing pore location
    these are set for pore(n+1) when pore(n) is generated
    '''
    def __init__(self, option):
        self.Locations = {
            'Circle' : 1,
            'circle' : 1,
            'Radial' : 1,
            'radial' : 1,
            1 : 1,
            'Square' : 2,
            'square' : 2,
            2 : 2}

        if option.LocationType not in self.Locations:
            pass
        
        elif self.Locations[option.LocationType] in {1, 2}:
            self.Square = self.Locations[option.LocationType] == 2
            
            self.it = 1 + option.ignoreborder*self.Square
            self.iu = 1 + option.ignoreborder*self.Square
            self.SpaceList = np.linspace(0, option.ArraySize, int(option.ArraySize/option.Spacing))
            self.AngleList = np.linspace(0, 2*np.pi, self.it*int(np.sqrt(option.ArraySize/option.Spacing)))
        else:
            raise(ValueError('XYprimer.Locations dictionary contains references to non-0 or -1 values.'))
        
        self.ignore_target_porosity = 0
        self.grid_complete = 0
        
        # controls failsafe for low option.location_err and False option.ignore_target_porosity
        self.iter = 0
        # max times grid can be filled before aborting further generation attempts
        self.max_iters = 5
