"""
	Program Name Goes Here
	Zachary Toth, August 2024

    Description goes here.
    
        namestyle
            # 'Timestamp' or name to be used
                'Timestamp' saves file as 'YYYY_MM_DD_hh_mm_ss.tif'
                namestyle = 'file_name.tif'
        
        param_file
            # File passes set of parameters to use
            # Choose one from ./param_files
            # or None to generate from LoadParameters.py file
        
        exports
            # list of filetypes to export
            # subset of ['tiff', 'txt', 'stl', 'xlsx']
            # txt files exported can be used as parameter files
                    
        rng_method
            # Set to False to use random seed.
            # Set to True to keep the same seed used last generation. 
            # Any other value will be used as a seed.

	Usage:
		Python:
            >>> from PoreGenerator_main import cortical_bone_generator
            >>> exports = ['tiff', 'txt', 'stl', 'xlsx']
            >>> rng_method = 'some seed value'
			>>> Bone = cortical_bone_generator(None, 'Timestamp', exports, rng_method)

		Terminal:
			$  python PoreGenerator_main.py Defaults.txt  'Timestamp'  ['tiff'] 'seed value'
"""

import os.path
import random
import numpy as np
import pandas as pd
import tifffile as tf

from numpy import cos, sin, tan

import PoreGenerator_funcs as PGf
import PoreGenerator_classes as PGc
from LoadParameters import LoadParameters

# %% Initialization

BEEP = False

def cortical_bone_generator(param_file = None, namestyle = 'Timestamp', exports = ['tiff'], rng_method = False):

    option, target_porosity, export, mu, sigma, weighting, params = LoadParameters(param_file, exports)
    option, mu, sigma = PGf.normalizeParameters(option, mu, sigma)
    # set file name
    [fpath, fname] = PGf.nameFig(namestyle)
    
    # sets rng based on rng_method
    RNGkey = PGf.setRNG(rng_method)
    
    # Initializes array of proper size
    Bone = np.ones((option.ArraySize, option.ArraySize, option.ArraySize), dtype=np.float32)
    # Sets up indexing of Bone array for use in calculations
    [ii,ij,ik] = np.unravel_index(np.arange(option.ArraySize**3), [option.ArraySize, option.ArraySize, option.ArraySize], 'F')
    
    PD = PGf.getPD(mu, sigma, weighting, option)
    # PD = PGc.probability_dist(mu, sigma, weighting, option)
    
    # Chooses target_porosity Value from experimental distribution
    if option.experimental_porosity is True:
        target_porosity = PD.porosity.rvs(1)[0]
    
    # readjusts target_porosity to account for loss when smoothPores is used
    if option.smoothPores is True:
        target_porosity = target_porosity/option.TP_CORRECTION_FACTOR
    
    # creates log for use in pore networking
    valueslog = np.zeros([12, int(round(7_000_000*target_porosity/mu.osteonlength, ndigits=-3))])
        # [R; C; theta; phi; x; y; minz; z; maxz; isfilled; A; B]
    iteration = 0
    
    # Pre-sets some variables for choosing pore location
    XYprimer = PGc.XYprimer(option)
    
    
    # %% Main Body
    if export.tiff or export.stl:
        while ((1-np.mean(Bone) < target_porosity) and (XYprimer.ignore_target_porosity == 0)) or (XYprimer.grid_complete == 0):
        
            [R, C] = PGf.getRC(option, PD)
            z = random.random() * option.ArraySize
            # theta is angle of trajectory, phi is angle of depression
            theta = PD.theta.rvs(1)
            phi = PD.phi.rvs(1)
            
            # prevent divide-by-zero errors for horizontal pores
            if phi == np.pi/2:
                phi = np.pi/2 - 0.1
        
            # Correction factor (phi is also correcting things in eqn)
            C = C/cos(phi)
            # cuts off cylinder
            minz = z - cos(phi)*0.5*PD.osteonlength.rvs(1)[0]
            maxz = z + cos(phi)*0.5*PD.osteonlength.rvs(1)[0]
        
            if sum(valueslog[10,:]) > params.pores_before_networking and random.random() > params.sealed_osteon_chance:
                [x,y,minz,z,maxz,valueslog,iteration] = PGf.networkPore(valueslog,minz,z,maxz,iteration)
            else:
                [x, y, XYprimer] = PGf.getXY(option, XYprimer)
        
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
                Bone = PGf.poreBlast(Bone)
            else:
                Bone = PGf.poreClast(Bone)
    
        # reverts target_porosity for bookkeeping
        target_porosity = target_porosity*option.TP_CORRECTION_FACTOR
    
    Bone = ~(Bone.astype(bool))
    porosity = np.mean(Bone)
    
    option, mu, sigma = PGf.revertParameters(option, mu, sigma)
    
    
    # %% Save Results
    
    if os.path.isdir(fpath) is False:
        os.makedirs(fpath)
    
    fullcell, sheetprep, sheetcell = PGf.getTextOutput(option, mu, sigma, weighting, params, target_porosity, porosity, RNGkey, fname)
    if export.xlsx is True:
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
        PGf.make3DModel(fpath, fname, Bone)
    
    return Bone

#%% Run from Terminal ############################################################################
if __name__ == "__main__":
    import sys
    if (len(sys.argv) > 1): # At least 1 command line parameter
        param_file = str(sys.argv[1])
        
        if (len(sys.argv) > 2): # At least 2
            namestyle = int(sys.argv[2])
            
            if (len(sys.argv) > 3):
                exports = int(sys.argv[3])
                
                if (len(sys.argv) > 4):
                    rng_method = int(sys.argv[4])
                    
                else: rng_method = False
                
            else: exports = ['tiff']
            
        else: namestyle = 'Timestamp'
        
    else: param_file = None

    cortical_bone_generator(param_file, namestyle, exports, rng_method)
    
    if BEEP:
        from winsound import Beep
        Beep(500,500)

