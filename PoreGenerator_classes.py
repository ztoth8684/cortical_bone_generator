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
    #%%

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
        
        elif self.Locations[option.LocationType] == 1:
            self.it = 1
            self.iu = 1
            self.AngleList = np.linspace(0, 2*np.pi, XYprimer.it*int(np.sqrt(option.ArraySize/option.Spacing)))
            self.SpaceList = np.linspace(0, option.ArraySize, int(option.ArraySize/option.Spacing))
            
        elif self.Locations[option.LocationType] == 2:
            self.it = 1+option.ignoreborder
            self.iu = 1+option.ignoreborder
            self.SpaceList = np.linspace(0, option.ArraySize, int(option.ArraySize/option.Spacing))
        
        self.ignore_target_porosity = 0
        self.grid_complete = 0
#%%

class probability_dist:
    def __init__(self, mu, sigma, weighting, option):

        # Probability distributions for linked diameters and circularities 
        # Normal SED (small, regular pores)
        self.Ncircularity = stats.norm(loc=mu.Ncircularity, scale=sigma.Ncircularity)
        self.Ndiameter = stats.truncnorm(a=option.mindiameter, b=np.Inf, loc=mu.Ndiameter, scale=sigma.Ndiameter)
        # High SED (larger, irregular pores)
        self.Hcircularity = stats.norm(loc=mu.Hcircularity, scale=sigma.Hcircularity)
        self.Hdiameter = stats.truncnorm(a=option.mindiameter, b=np.Inf, loc=mu.Hdiameter, scale=sigma.Hdiameter)
        # SED distribution
        self.SED = stats.norm(loc=mu.SED, scale=sigma.SED)
        
        # Probability distributions for unlinked diameters and circularities    
        self.TOTdiameter = MixtureModel([self.Ndiameter, self.Hdiameter], weighting.SED)
        self.TOTcircularity = MixtureModel([self.Ncircularity, self.Hcircularity], weighting.SED)
        
        # Probability distributions for number/ length of pores
        self.porosity = stats.truncnorm(a=0.01, b=np.Inf, loc=mu.porosity, scale=sigma.porosity)
        self.osteonlength = stats.truncnorm(a=1, b=option.maxosteonlength, loc=mu.osteonlength, scale=sigma.osteonlength)
        
        # Probability distributions for azimuthal angle
        if weighting.phi_values == 'rand':
            self.phi = stats.uniform(loc=0, scale = 0.5*np.pi)
        elif len(weighting.phi_values) == 1:
            self.phi = stats.uniform(loc=weighting.phi_values, scale = 0)
        else:
            self.phi = len(weighting.phi_probs)*[0]
        
            for n in range(len(self.phi)):
                self.phi[n] = stats.uniform(loc=weighting.phi_values[n], scale=weighting.phi_values[n+1]-weighting.phi_values[n])
            self.phi = MixtureModel(self.phi,weights= weighting.phi_probs)
        
        # Probability distributions for radial angle
        if weighting.theta_values == 'rand':
            self.theta = stats.uniform(loc=0, scale = 2*np.pi)
        elif len(weighting.theta_values) == 1:
            self.theta = stats.uniform(loc=weighting.theta_values, scale = 0)
        else:
            self.theta = len(weighting.theta_probs)*[0]
            
            for n in range(len(self.theta)):
                self.theta[n] = stats.uniform(loc=weighting.theta_values[n], scale=weighting.theta_values[n+1]-weighting.theta_values[n])
            self.theta = MixtureModel(self.theta,weights= weighting.theta_probs)