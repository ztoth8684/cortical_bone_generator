# Cortical Bone Generator

# Components
## PoreGenerator_MK6.py

The main program. Generates a 3-dimensional array intended to replicate cortical bone at the micron scale. Each voxel of the array generated represents a 10 micron cube. Uses functions defined in PoreGenerator_funcs.py
## PoreGenerator_funcs.py
### nameFig(option)
Function that takes in option.namestyle and outputs the name to be used for the files.
### getPD(mu, sigma, weighting, option)
Function that takes in parameters and outputs probability distribution objects, which are later used to create variance within the program.
### getRC(option, PD)
Function that generates radius and circularity values for each pore generated. If DiscreteCircularities or DiscreteDiameters options are used, the values used are selected from that discrete list.
### networkPore(valueslog, minz, z, maxz, iteration)
Function that takes in the list of pores previously generated in order to generate new pores branching off of previous ones. Used to better replicate the networked structure of the canal system in cortical bone.
### getXY(option, XYprimer)
Function that generates x and y coordinate values for each pore generated. Depending on option.LocationType, each pore is either generated with a random center, or according to a square or radial grid pattern.
### erodePores(Bone)
Functions used in conjunction to roughen the image generated and to break up geometric artefacts such as straight lines.
### mergePores(Bone)
Functions used in conjunction to roughen the image generated and to break up geometric artefacts such as straight lines.
### getTextOutput(option, mu, sigma, weighting, TargetPorosity, pores_before_networking, sealed_osteon_chance, transverse_flag_onset, shape_proportions, RNGkey, fname)
Function that takes in all important variables and formats them to be saved as text and excel files along with the image generated.
### make3DModel(fpath, fname, Bone)
Lorem ipsum dolor sit amet...
## GetPoreData.py
### GetPoreInfo(file_directory)
...
### SavePoreData(filename)
...
### ReadPoreData(targetfile)
...
### PoreHist(values, title)
...
### PlotDiameterGraphs(titles, varnames)
...
### UnbalancedTTest(values_test, values_real)
...
### PercentBinHeightChange(Best_values, Literature_values, Real_values)
...
### importBone(fpath, fname)
...
### RuntimeProfiler(filename)
...

## SetupBoneImg.py
Python script for use in Seg3D which adds a function, `showPorosity(filename)`, to be used in the python console. This function takes the name of a file to import and processes the image before displaying it in 3D.
Also adds function `iso()`, which rotates the figure to an isometric viewing angle.

To load the script in Seg3D, you can use the following command in the console:
```py
exec(open('<file_location>/SetupBoneImg.py').read())
```
## GetImagePorosity3D.ijm
ImageJ macro that calculates the porosity of a 3D TIF stack. Can be used to get the bone porosity from a CT scan for use in the TargetPorosity parameter.
# Settings / Parameters
Below is a list of all settings and parameters that can be changed to vary the image generated.  Each section is headed with the name of the MATLAB struct object that contains those variables. Each variable's entry begins with the default value and which program(s) it is present in.
## `option` struct
### `debug`
###### Default: false     (2D | 3D)
Whether to use the same RNG values used the last time the program was run. If all parameters and settings are kept the same, this will result in a duplicate image being generated. This can be used to vary settings and see their effect on the same image, though this might not work for all settings.

Can also be set to an integer value or string in order to use that as the RNG seed.
### `quickgen`
###### Default: true     (2D | 3D)
In the 3D program, whether to open the generated file in ImageJ and to minimize MATLAB when the program finishes running. For this to work, ImageJ must be set as the default application for opening .tif files.

In the 2D program, setting this to true also bypasses the dialog popups and determines whether to save the file.
### `namestyle`
###### Default: 'Trial.tif'     (2D | 3D)
What to name the file(s) generated. If this is set as 'Timestamp,' the file will be named YYYY_MM_DD_HH_MM_N. For any other string, the file will be named the same as the string. ".tif" can be appended or omitted from the name. 

A text file will also be generated with the same name, containing all parameters used.
### `varlink`
###### Default: false     (2D | 3D)
If true, pore radius and circularity will be linked such that radius values significantly different from the mean will be paired with circularity values significantly different from the mean. In effect, larger pores will also tend to be more oblong.
### `mindiameter`
###### Default: 0     (2D | 3D)
Sets a lower bound for the diameter of pores generated.
### `LocationType`
###### Default: 0     (2D | 3D)
Defines the method that determines pore location. If set to 'Circle,' 'Radial,' or 1, pores are generated in a radial pattern. If set to 'Square' or 2, pores are generated in a square grid. Otherwise, x and y coordinates are selected randomly. Z coordinates for each pore are random regardless of this setting.

This setting is more dramatic when used in 2D.
Note: TargetPorosity may be overshot in order to construct a complete grid.
### `Spacing`
###### Default: 16     (2D | 3D)
If option.LocationType is not random, the spacing beween pores. For a square grid, this is in pixels. 
### `location_err`
###### Default: 6     (2D | 3D)
If option.LocationType is not random, the magnitude of the random offset between grid lines and where pores are actualyy generated. If this value is less than 6, TargetPorosity will be ignored.
### `ignoreborder`
###### Default: false     (2D | 3D)
If option.LocationType is 2 (square grid), whether to generate pores along the border of the image.
### `LinearDiscreteDiameters`
###### Default: []     (2D | 3D)
If this array is not empty, the diameter of each pore is instead chosen randomly from this list.

Mutually exclusie with option.WeightedDiscreteDiameters.
### `WeightedDiscreteDiameters`
###### Default: []     (2D | 3D)
If this array is not empty, the diameter of each pore is generated from distributions as normal, then rounded to the nearest value present in this list.

Mutually exclusive with option.LinearDiscreteDiameters.
### `LinearDiscreteCircularities`
###### Default: []     (2D | 3D)
If this array is not empty, the circularity of each pore is instead chosen randomly from this list.

Mutually exclusive with option.WeightedDiscreteCircularities.
### `WeightedDiscreteCircularities`
###### Default: []     (2D | 3D)
If this array is not empty, the circularity of each pore is generated from distributions as normal, then rounded to the nearest value present in this list.

Mutually exclusive with option.LinearDiscreteCircularities.
### `mergePores`
###### Default: true     (2D | 3D)
In the 2D program, calls mergePores.m once in order to merge nearby pores, adding porosity to spaces adjacent to highly porous areas.

In the 3D program, calls mergePores.m and erodePores.m multiple times in order to smooth geometric boundaries and make shapes generated look more natural.
### `variedPoreShape`
###### Default: false     (3D)
Adds variation in the shape of pores, using ellipsoids, hyperboloids, and cones, in addition to cylinders. Ratios to be used for each shape is stored in the shape_proportions variable.
### `ArraySize`
###### Default: 200     (2D | 3D)
Number of pixels/voxels generated on each dimension of the image. At the default, the 3D program will generate a 200 by 200 by 200 voxel image. 
### `maxosteonlength`
###### Default: 220/3     (3D)
The maximum length for a pore to generate in the z direction. Works in tandem with mu.osteonlength and sigma.osteonlength to determine how long a given pore should be.
## `mu` / `sigma` structs
The first default value listed is the average (mu), while the second value is the deviation (sigma).
Default values for all diameter and circularity-raleted parameters from https://doi.org/10.1002/jbmr.3561
### `SED`
###### Defaults: 8, 4.5     (2D | 3D)
Parameters for the frequency at which pores generate with extreme values for diameter and circularity. This is used while option.varlink is true, whereas weighting.SED is used when option.varlink is false.
### `Ndiameter`
###### Defaults: 4.8, 0.4     (2D | 3D)
Parameters for diameter when more moderate values are chosen by SED processes.
### `Ncircularity`
###### Defaults: 0.66, 0.03     (2D | 3D)
Parameters for circularity when more moderate values are chosen by SED processes.
### `Hdiameter`
###### Defaults: 15.7, 3.75     (2D | 3D)
Parameters for diameter when more extreme values are chosen by SED processes.
### `Hcircularity`
###### Defaults: 0.42, 0.045     (2D | 3D)
Parameters for circularity when more extreme values are chosen by SED processes.
### `osteonlength`
###### Defaults: 100, 75/2     (3D)
Parameters for the length of a given pore in the z direction


### `porosity`
###### Defaults: 0.075046512, 0.036744908     (3D)
Parameters for experimental distribution of porosity values.
## `weighting` struct
### `SED`
###### Default: [0.9774 0.0226]     (2D | 3D)
Parameters for the frequency at which pores generate with extreme values for diameter and circularity. This is used while option.varlink is false, whereas mu.SED and sigma.SED are used when option.varlink is true.
### `phi_values`
###### Default: [0, pi/12, pi/2]     (3D)
Parameters for the angles at which pores are generated. With default `phi_probs`, most pores generate at an angle less than the second value of this array. 

Note: use [0,10^-161] with [0,1] for all vertical pores.
### `phi_probs`
###### Default: [0, 0.75, 1]     (3D)
Parameters for the angles at which pores are generated. The second value of this array is the chance for a given pore to generate at an angle less than the second value of weighting.phi_values. All other pores genrate between the second and third values of weighting.phi_values.
## Other Variables
### `TargetPorosity`
###### Default: 0.026     (2D | 3D)
The porosity at which the image generates, represented as a percent of pixels/voxels with a value of 0 (porous in the image).
### `shape_proportions`
###### Default: [0.392, 0.094, 0.351, 0.122, 0.041]     (3D)
Used with option.variedPoreShape to supply the chances of each pore shape generating. The values of this array are as follows: cylinder, proximal opening cone, distal opening cone, ellipsoid, hyperboloid.
The sum of this array should be 1.

Default values from https://doi.org/10.1111/j.1439-0264.2009.00973.x
### `pores_before_networking`
###### Default: 75     (3D)
The number of pores to generate with random x and y coordinates before networkPore.m begins being used to form a networked structure.
### `sealed_osteon_chance`
###### Default: 0.068     (3D)
The chance for a given pore to generate with random x and y coordinates after the pores_before_networking threshold has been reached.

Default values from https://doi.org/10.1002/ar.21309
### `transverse_flag_onset`
###### Default: pi/4     (3D)
The maximum phi value, after which a pores length in the x and y directions is limited.
***
# Using the Program
The program was written with MATLAB R2021a, but should work for all versions after R2019a. 

To setup the program for use, all matlab files should be placed in your `Documents/MATLAB` folder. A folder named "PoreFigs" should also be created in this folder (`Documents/MATLAB/PoreFigs`). This folder is where all files generated by the program will be saved. If you wish to rename this folder, you will also need to change line 5 of nameFig.m and lines 35, 37, and 53 of call.m which specify this folder.

