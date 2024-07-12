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
        
        elif self.Locations[option.LocationType] in {1, 2}:
            self.Square = self.Locations[option.LocationType] == 2
            
            self.it = 1 + option.ignoreborder*self.Square
            self.iu = 1 + option.ignoreborder*self.Square
            self.SpaceList = np.linspace(0, option.ArraySize, int(option.ArraySize/option.Spacing))
            self.AngleList = np.linspace(0, 2*np.pi, self.it*int(np.sqrt(option.ArraySize/option.Spacing)))
        
        self.ignore_target_porosity = 0
        self.grid_complete = 0
#%%
'''
class probability_dist:
    def trunc_dist(self, mu, sigma, lower, upper):
        # reindexes limits according to truncnorm documentation
        a = (lower - mu)/sigma
        b = (upper - mu)/sigma
        dist = stats.truncnorm(a=a, b=b, loc=mu, scale=sigma)
        
        return dist
    
    def span_dist(self, values, probs, rand_range):
        ''''''
        Creates a distribution broken into weighted self-uniform spans
        
        values defines the bounds between each span
        probs defines the weightings of each span
        rand_range is a two-length vector that defines the 'random' generation limits
        ''''''
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
    
    def __init__(self, mu, sigma, weighting, option):

        # Probability distributions for linked diameters and circularities 
        # Normal SED (small, regular pores)
        self.Ncircularity = stats.norm(loc=mu.Ncircularity, scale=sigma.Ncircularity)
        self.Ndiameter = self.trunc_dist(option.mindiameter, np.Inf, mu.Ndiameter, sigma.Ndiameter)
        # High SED (larger, irregular pores)
        self.Hcircularity = stats.norm(loc=mu.Hcircularity, scale=sigma.Hcircularity)
        self.Hdiameter = self.trunc_dist(option.mindiameter, np.Inf, mu.Hdiameter, sigma.Hdiameter)
        # SED distribution
        self.SED = stats.norm(loc=mu.SED, scale=sigma.SED)
        
        # Calculates weighting split for unlinked distributions
        split = lambda num: [num, 1-num]
        SED_weighting = split(stats.norm(mu.SED, sigma.SED).cdf(option.SED_limit))
        
        # Probability distributions for unlinked diameters and circularities    
        self.TOTdiameter = MixtureModel([self.Ndiameter, self.Hdiameter], SED_weighting)
        self.TOTcircularity = MixtureModel([self.Ncircularity, self.Hcircularity], SED_weighting)
        
        # Probability distributions for number/ length of pores
        # Minimum porosity clipped at 0.01
        self.porosity = self.trunc_dist(mu.porosity, sigma.porosity, 0.01, np.Inf)
        # Maximum clipped at maxosteonlength
        self.osteonlength = self.trunc_dist(mu.osteonlength, sigma.osteonlength, 0, option.maxosteonlength)
        
        # Probability distributions for azimuthal angle
        self.phi = self.span_dist(weighting.phi_values, weighting.phi_probs, [0, 0.5*np.pi])
        
        # Probability distributions for radial angle
        self.theta = self.span_dist(weighting.theta_values, weighting.theta_probs, [0, 2*np.pi])
'''