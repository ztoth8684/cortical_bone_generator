# Cortical Bone Generator

# Components
## PoreGenerator_main.py

The main program. Generates a 3-dimensional array intended to replicate cortical bone at the micron scale. Each voxel of the array generated represents a 10 micron cube. Uses functions defined in PoreGenerator_funcs.py
## PoreGenerator_funcs.py
### nameFig(option)
Function that outputs the name to save the program output files as. If `Timestamp`, saves the file in `YYYY_MM_DD_HH_MM_SS` format.
### setRNG(option)
Sets RNG seed value used by the program. Uses a random seed id set to False, uses the last used seed if set to True, otherwise uses passed string as the seed.
### getRC(option, PD)
Function that generates radius and circularity values for each pore generated. If `LinearDiscrete*` or `WeightedDiscrete*` options are used, the values used are chosen from the list or rounded to the nearest value on the list, respectively.
### networkPore(valueslog, minz, z, maxz, iteration)
Function that uses a list of pores previously generated in order to generate new pores branching off of them. Used to better replicate the networked structure of the canal system in cortical bone.
### getXY(option, XYprimer)
Function that generates x and y coordinate values for each pore generated. Depending on `option.LocationType`, each pore is either generated with a random center, or according to a square (`'Square'`) or radial grid pattern (`'Radial'`). 
### poreBlast(Bone)
Function used to roughen the image generated and to break up geometric artefacts such as straight lines by filling in pores at the edges.
### poreClast(Bone)
Function used to roughen the image generated and to break up geometric artefacts such as straight lines by widening pores and merging nearby ones.
### getTextOutput(...)
Function that takes in all input variables and formats them to be saved as text and excel files along with the image generated.
### make3DModel(fpath, fname, Bone)
Saves generated array as an .STL file.
## PoreGenerator_classes.py
### MixtureModel
Extension of `scipy.stats.rv_continuous` that allows for combined distributions.
### XYprimer(option)
Framework for `option.LocationType` grid generation.
### probability_dist(mu, sigma, weighting, option)
Probability distribution objects for for circularity, diameter, porosity, osteon length, azimuthal angle, and radial angle, which are later used to create variance within the program.
## GetPoreInfo.py
### importBone(fpath, fname)
Reads .TIFF file to retrieve the saved array.
### GetPoreData(file_directory)
Imports folder of .TIFF files, segments them into two-dimensional layers, and labels all pixels in each pore via a connected component labeling algorithm. The diameter of each pore is then found assuming a circular shape. Both the lists for individual samples and the aggregate data is saved.
### SavePoreData(filename)
Saves `counts` and `values` variables from `GetPoreData()` as a .pkl file.
### ReadPoreData(targetfile)
Opens .pkl file from `SavePoreData()` and retrieves the `values` data.
### PoreHist(values, title)
Creates histogram with standardized format from pore data.
### PlotDiameterGraphs(titles, varnames)
Batch creates histograms with `PoreHist()`.
### UnbalancedTTest(values_test, values_real)
Performs an unbalanced t-test on two sets of values.
### PercentBinHeightChange(Best_values, Literature_values, Real_values)
Compares the histogram bins for two samples, returning the average percent change in bin height between the two.
### RuntimeProfiler(filename)
Runs the program with data on how long different code segments take to run.

## SetupBoneImg.py
Python script for use in Seg3D which adds a function, `showPorosity(filename)`, to be used in the python console. This function takes the name of a file to import and processes the image before displaying it in 3D.
Also adds function `iso()`, which rotates the figure to an isometric viewing angle and `top()` which rotates the figure to a top-down view.

To load the script in Seg3D, you can use the following command in the console:
```py
exec(open('<file_location>/SetupBoneImg.py').read())
```
## GetImagePorosity3D.ijm
ImageJ macro that calculates the porosity of a 3D TIF stack. Can be used to get the bone porosity from a CT scan for use in the TargetPorosity parameter.
# Settings / Parameters
Below is a list of all settings and parameters that can be changed to vary the image generated.  Each section is headed with the name of the class structure that contains those variables. Each variable's entry begins with the default value and which program(s) it is present in.
### `option.rng_method`
###### Default: false
Whether to use the same RNG values used the last time the program was run. If all parameters and settings are kept the same, this will result in a duplicate image being generated. This can be used to vary settings and see their effect on the same image, though this might not work for all settings.

Can also be set to an integer value or string in order to use that as the RNG seed.
### `option.namestyle`
###### Default: 'Timestamp'
What to name the file(s) generated. If this is set as 'Timestamp,' the file will be named YYYY_MM_DD_HH_MM_N. For any other string, the file will be named the same as the string. ".tif" can be appended or omitted from the name. 

A text file will also be generated with the same name, containing all parameters used.
### `option.varlink`
###### Default: false
If true, pore radius and circularity will be linked such that radius values significantly different from the mean will be paired with circularity values significantly different from the mean. In effect, larger pores will also tend to be more oblong.
### `option.mindiameter`
###### Default: 0
Sets a lower bound for the diameter of pores generated.
### `option.LocationType`
###### Default: 0
Defines the method that determines pore location. If set to 'Circle,' 'Radial,' or 1, pores are generated in a radial pattern. If set to 'Square' or 2, pores are generated in a square grid. Otherwise, x and y coordinates are selected randomly. Z coordinates for each pore are random regardless of this setting.

This setting is more dramatic when used in 2D.
Note: TargetPorosity may be overshot in order to construct a complete grid.
### `option.Spacing`
###### Default: 16
If option.LocationType is not random, the spacing beween pores. For a square grid, this is in pixels. 
### `option.location_err`
###### Default: 6
If option.LocationType is not random, the magnitude of the random offset between grid lines and where pores are actualyy generated. If this value is less than 6, TargetPorosity will be ignored.
### `option.ignore_target_porosity`
...
### `option.ignoreborder`
###### Default: false
If option.LocationType is 2 (square grid), whether to generate pores along the border of the image.
### `option.LinearDiscreteDiameters`
###### Default: []
If this array is not empty, the diameter of each pore is instead chosen randomly from this list.

Mutually exclusie with `option.WeightedDiscreteDiameters`.
### `option.WeightedDiscreteDiameters`
###### Default: []
If this array is not empty, the diameter of each pore is generated from distributions as normal, then rounded to the nearest value present in this list.

Mutually exclusive with `option.LinearDiscreteDiameters`.
### `option.LinearDiscreteCircularities`
###### Default: []
If this array is not empty, the circularity of each pore is instead chosen randomly from this list.

Mutually exclusive with `option.WeightedDiscreteCircularities`.
### `option.WeightedDiscreteCircularities`
###### Default: []
If this array is not empty, the circularity of each pore is generated from distributions as normal, then rounded to the nearest value present in this list.

Mutually exclusive with `option.LinearDiscreteCircularities`.
### `option.smoothPores`
###### Default: true
Calls `poreBlast` and `poreClast` multiple times in order to smooth geometric boundaries and make shapes generated look more natural.
### `option.variedPoreShape`
###### Default: false
Adds variation in the shape of pores, using ellipsoids, hyperboloids, and cones, in addition to cylinders. Ratios to be used for each shape is stored in the `param.shape_proportions` variable.
### `option.ArraySize`
###### Default: 200
Number of pixels/voxels generated on each dimension of the image. At the default, the 3D program will generate a 200 by 200 by 200 voxel image. 
### `option.maxosteonlength`
###### Default: 220/3
The maximum length for a pore to generate in the z direction. Works in tandem with mu.osteonlength and sigma.osteonlength to determine how long a given pore should be.
### `option.SED_limit`
...
### `mu/sigma .SED`
###### Defaults: `mu`=8, `sigma`=4.5
Parameters for the frequency at which pores generate with extreme values for diameter and circularity. This is used while option.varlink is true, whereas weighting.SED is used when option.varlink is false.
### `mu/sigma .Ndiameter`
###### Defaults: `mu`=4.8, `sigma`=0.4
Parameters for diameter when more moderate values are chosen by SED processes.
### `mu/sigma .Ncircularity`
###### Defaults: `mu`=0.66, `sigma`=0.03
Parameters for circularity when more moderate values are chosen by SED processes.
### `mu/sigma .Hdiameter`
###### Defaults: `mu`=15.7, `sigma`=3.75
Parameters for diameter when more extreme values are chosen by SED processes.
### `mu/sigma .Hcircularity`
###### Defaults: `mu`=0.42, `sigma`=0.045
Parameters for circularity when more extreme values are chosen by SED processes.
### `mu/sigma .osteonlength`
###### Defaults: `mu`=100, `sigma`=37.5
Parameters for the length of a given pore in the z direction
### `mu/sigma .porosity`
###### Defaults: `mu`=0.075046512, `sigma`=0.036744908
Parameters for experimental distribution of porosity values.
### `weighting.phi_values`
###### Default: [0, pi/12, pi/2]
Parameters for the angles at which pores are generated. With default `phi_probs`, most pores generate at an angle less than the second value of this array. 

Note: use [0,10^-161] with [0,1] for all vertical pores.
### `weighting.phi_probs`
###### Default: [0, 0.75, 1]
Parameters for the angles at which pores are generated. The second value of this array is the chance for a given pore to generate at an angle less than the second value of weighting.phi_values. All other pores genrate between the second and third values of weighting.phi_values.
### `weighting.theta_values`
###### Default: [0, 2pi]
Parameters for the angles at which pores are generated. With default `phi_probs`, most pores generate at an angle less than the second value of this array. 

Note: use [0,10^-161] with [0,1] for all vertical pores.
### `weighting.theta_probs`
###### Default: [1]
Parameters for the angles at which pores are generated. The second value of this array is the chance for a given pore to generate at an angle less than the second value of weighting.phi_values. All other pores genrate between the second and third values of weighting.phi_values.
### `params.shape_proportions`
###### Default: [0.392, 0.094, 0.351, 0.122, 0.041]
Used with option.variedPoreShape to supply the chances of each pore shape generating. The values of this array are as follows: cylinder, proximal opening cone, distal opening cone, ellipsoid, hyperboloid.
The sum of this array should be 1.

Default values from https://doi.org/10.1111/j.1439-0264.2009.00973.x
### `params.pores_before_networking`
###### Default: 75
The number of pores to generate with random x and y coordinates before networkPore.m begins being used to form a networked structure.
### `params.top_branches`
...
### `params.bottom_branches`
...
### `params.sealed_osteon_chance`
###### Default: 0.068
The chance for a given pore to generate with random x and y coordinates after the pores_before_networking threshold has been reached.

Default values from https://doi.org/10.1002/ar.21309
### `params.transverse_flag_onset`
###### Default: pi/4
The maximum phi value, after which a pores length in the x and y directions is limited.
### `target_porosity`
###### Default: 'Exp'
The porosity at which the image generates, represented as a percent of pixels/voxels with a value of 0 (porous in the image).
***
# Using the Program
The program was written with Python 3.11
The program uses the following modules:
- numpy
- pandas
- tifffile
- meshlib (https://github.com/MeshInspector/MeshLib)
- 
