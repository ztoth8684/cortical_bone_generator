# Cortical Bone Generator
The program was written with Python 3.11. The program uses the following modules:
- numpy
- pandas
- scipy
- meshlib
- tifffile
- matplotlib
- skimage

Usage:
	Python:
        >>> `from PoreGenerator import cortical_bone_generator`
        >>> `param_file = 'Optimized.txt'`
        >>> `namestyle = 'Timestamp'`
        >>> `exports = ['tiff', 'txt', 'stl', 'xlsx']`
        >>> `rng_method = 'some seed value'`
        >>> `fpath = './pore_files/subfolder/'`
		>>> `Bone = cortical_bone_generator(param_file, namestyle, exports, rng_method, fpath)`
	Terminal:
		$  `python PoreGenerator.py Defaults.txt  Timestamp  ['tiff'] 'seed string' './pore_files/subfolder/'`

Default values are as follows:
	`param_file = None`
		This is the file stating what values of parameters to use. If "None", values from LoadParameters.py are used.
	`namestyle = 'Timestamp'`
		Filename for the generated files. If "Timestamp", the file is named based on the time the core started generating. If "prefix(\*)", the timestamp is prefixed with \*.
	`exports = ['tiff']`
		A list of filetypes to export as. Options are "xlsx", "txt", "tiff", "stl".
	`rng_method = False`
		Sets rng seed value. If "True", re-uses the value from the last time the program was run. If "False", uses a random value. Else, uses the value of "rng_method" as the seed. 
	`fpath = './pore_files/'`
		Directory to save files to.
# Components
## PoreGenerator.py

The main program. Generates a 3-dimensional array intended to replicate cortical bone at the micron scale. Each voxel of the array generated represents a 10 micron cube. Uses functions defined in PoreGenerator_funcs.py, classes defined in PoreGenerator_classes.py, and LoadParameters.py to import parameters.
### Arguments:
param_file
    File that passes set of parameters to use
    Choose a file from `./param_files` folder or pass `None` to generate from the LoadParameters.py file
    
namestyle
    Either 'Timestamp' or the name to be used
    'Timestamp' saves the file in `YYYY_MM_DD_hh_mm_ss.tif` format.
    Format : 'file_name.tif'
    
exports
    A list of filetypes to export
    Use a subset of ['tiff', 'txt', 'stl', 'xlsx']
    txt files exported can be used as parameter files.
    
rng_method
    Set to `False` to use a random seed.
    Set to `True` to keep the same seed used during the last generation.
    Any other value will be used as a seed.
## PoreGenerator_funcs.py
### chooseExports
##### Input: `exports`
Generates program-readable values from `exports` program argument.
### scaleParameters
##### Inputs: `option`, `mu`, `sigma`, `scale`
Scales all length-dimension parameters by `scale` value.

If scale is 1/10: Converts inputs from µm to voxel units
If scale is 10: Converts inputs back from voxel to µm units
### nameFig
##### Input: `namestyle`
Outputs the name to save the program output files as. If `namestyle` is 'Timestamp', saves the file in `YYYY_MM_DD_HH_MM_SS` format.
### setRNG
##### Input: `rng_method`
Sets RNG seed value used by the program. Uses a random seed if set to False, uses the last used seed if set to True, otherwise uses passed string as the seed.
### getPD
##### Inputs: `mu`, `sigma`, `weighting`, `option`
Generates probability distributions used by the program.
### getRC
##### Inputs: `option`, `PD`
Generates radius and circularity values for each pore generated. If `LinearDiscrete*` or `WeightedDiscrete*` options are used, the values used are chosen from the list or rounded to the nearest value on the list, respectively.
### networkPore(valueslog, minz, z, maxz, iteration)
##### Inputs: `valueslog`, `minz`, `z`, `maxz`, `iteration`
Uses a list of pores previously generated in order to generate new pores branching off of them. Used to better replicate the networked structure of the canal system in cortical bone.
### getXY
##### Inputs: `option`, `XYprimer`
Generates x and y coordinate values for each pore generated. Depending on `option.LocationType`, each pore is either generated with a random center, or according to a square (`'Square'`) or radial grid pattern (`'Radial'`). 
### poreBlast
##### Input: `Bone`
Roughens the image generated and breaks up geometric artefacts such as straight lines by filling in pores at the edges.
### poreClast
##### Input: `Bone`
Roughens the image generated and breaks up geometric artefacts such as straight lines by widening pores and merging nearby ones.
### getTextOutput
##### Inputs: `option`, `mu`, `sigma`, `weighting`, `params`, `target_porosity`, `porosity`, `RNGkey`, `fname`
Takes in all input variables and formats them to be saved as text and excel files along with the image generated.
### make3DModel
##### Inputs: `fpath`, `fname`, `Bone`
Saves generated array as an .STL file.
## PoreGenerator_classes.py
### MixtureModel
Extension of `scipy.stats.rv_continuous` that allows for bimodal normal distributions.
### XYprimer
Framework for `option.LocationType` grid generation.
## LoadParameters.py
Populates all parameters needed to run the program. Either retrieves these from a parameter file or from the LoadParameters.py file itself if no `param_file` argument is passed. This is useful for testing different parameters. 

This file also contains explanations of all parameters.
## GetPoreInfo.py
### importBone(fpath, fname)
Reads .tif file to retrieve a saved array for further processing.
### rotateBone(fpath, fname)
Rotates bone file so the view plane is along the long axis. For use with ImageJ to better view generated .tif files.
### invertBone(fpath, fname)
Inverts bone file to swap black / white locations. For use with ImageJ to better view generated .tif files.
### GetPoreData(file_directory)
Imports folder of .tif files, segments them into two-dimensional layers, and labels all pixels in each pore via a connected component labeling algorithm. The diameter of each pore is then found assuming a circular shape. Both the lists for individual samples and the aggregate data are outputted.
### SavePoreData(counts, values, filename)
Saves `counts` and `values` variables from `GetPoreData()` as a .pkl file.
### ReadPoreData(targetfile)
Opens .pkl file from `SavePoreData()` and retrieves the `values` data.
### PoreHist(values, title)
Creates histogram with standardized format from pore data. This format includes 20 bins of size 3 µm, spanning 0–60 µm, and axis labels.
### PlotDiameterGraphs(targetfiles, titles)
Batch creates histograms directly from file using `PoreHist()`.
### PercentBinHeightChange(values_1, values_2)
Compares all histogram bins for two samples, returning the average percent change in bin height between the two.
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
Below is a list of all settings and parameters that can be changed to vary the image generated. This does not include the arguments of the main function.
### Radius-Circularity Correlation {`option.varlink`}
###### Type: bool
If true, pore radius and circularity will be linked such that radius values significantly different from the mean will be paired with circularity values significantly different from the mean. In effect, larger pores will also tend to be more oblong.
### Minimum Diameter {`option.mindiameter`}
###### Type: float
Sets a lower bound for the diameter of pores generated. Units are in µm.
### Generation Pattern {`option.LocationType`}
###### Type: string
Defines the method that determines pore location. If set to 'Radial,' pores are generated in a radial pattern. If set to 'Square,' pores are generated in a square grid. Otherwise, x and y coordinates are selected randomly. Z coordinates for each pore are random regardless of this setting.

Note: `target_porosity` may be overshot for 'Radial' or 'Square' methods in order to construct a complete grid for very low porosity values.
### Pore Spacing {`option.Spacing`}
###### Type: float
If `option.LocationType` is not random, the spacing between pores. Units are in µm.
### Grid Random Offset {`option.location_err`}
###### Type: float
If `option.LocationType` is not random, the magnitude of the random offset between grid lines and where pores are actually generated. Units are in µm.
### Repeat Grid Population {`option.ignore_target_porosity`}
###### Type: bool
If `option.LocationType` is not random, whether to stop pore generation after the grid has been completed.
### Grid Edge Population {`option.ignoreborder`}
###### Type: bool
If `option.LocationType` is 'Square', whether to generate pores along the border of the image.
### Uniform Diameter Choice {`option.LinearDiscreteDiameters`}
###### Type: list of floats
If this array is not empty, the diameter of each pore is chosen randomly from this list instead of from statistical distributions.

Mutually exclusive with `option.WeightedDiscreteDiameters`.
### Weighted Diameter Choice {`option.WeightedDiscreteDiameters`}
###### Type: list of floats
If this array is not empty, the diameter of each pore is generated from distributions as normal, then rounded to the nearest value present in this list.

Mutually exclusive with `option.LinearDiscreteDiameters`.
### Uniform Circularity Choice {`option.LinearDiscreteCircularities`}
###### Type: list of floats
If this array is not empty, the circularity of each pore is chosen randomly from this list instead of from statistical distributions.

Mutually exclusive with `option.WeightedDiscreteCircularities`.
### Weighted Circularity Choice {`option.WeightedDiscreteCircularities`}
###### Type: list of floats
If this array is not empty, the circularity of each pore is generated from distributions as normal, then rounded to the nearest value present in this list.

Mutually exclusive with `option.LinearDiscreteCircularities`.
### Pore Smoothing {`option.smoothPores`}
###### Type: bool
Calls `poreBlast` and `poreClast` multiple times in order to smooth geometric boundaries and make shapes generated look more natural.
### Pore Shape Variation {`option.variedPoreShape`}
###### Type: bool
Adds variation in the shape of pores, using ellipsoids, hyperboloids, and cones, in addition to cylinders. Ratios to be used of each shape is stored in the `param.shape_proportions` variable.
### Core Size {`option.ArraySize`}
###### Type: float
Side length in µm of cubic .tif file generated. The scale is 10 micrometers/voxel.
### Porosity from Distribution {`option.experimental_porosity`}
###### Type: bool
Whether to choose the porosity value from the distribution formed by `mu.porosity` and `sigma.porosity` instead of the inputted value.
### Porosity Correction Factor {`option.TP_CORRECTION_FACTOR`}
###### Type: float
`target_porosity` is divided by this value if `option.smoothPores` is True in order to compensate for the variation in final porosity caused by this routine.

Values to use for this parameter at a given porosity can be found by comparing target and actual porosity after generating such samples. 
### Extreme Pore Chance {`option.SED_limit`}
###### Type: float
Proportion of pores that are generated using 'normal' diameter and circularity distributions, rather than 'high' distributions.
### Diameter (Normal SED) {`mu.Ndiameter / sigma.Ndiameter`}
###### Type: float
Parameters for diameter when more moderate values are chosen by SED processes.
### Circularity (Normal SED) {`mu.Ncircularity / sigma.Ncircularity`}
###### Type: float
Parameters for circularity when more moderate values are chosen by SED processes.
### Diameter (High SED) {`mu.Hdiameter / sigma.Hdiameter`}
###### Type: float
Parameters for diameter when more extreme values are chosen by SED processes.
### Circularity (High SED) {`muHcircularity / sigma.Hcircularity`}
###### Type: float
Parameters for circularity when more extreme values are chosen by SED processes.
### Pore Length {`mu.osteonlength / sigma.osteonlength`}
###### Type: float
Parameters for the length of a given pore in the z direction
### Maximum Pore Length {`option.maxosteonlength`}
###### Type: float
The maximum length for a pore to generate in the z direction. Works with `mu.osteonlength` and `sigma.osteonlength` to determine how long a given pore should be.
### Porosity {`mu.porosity / sigma.porosity`}
###### Type: float
Parameters for experimental distribution of porosity values.
### Azimuthal Angles {`weighting.phi_values`}
###### Type: list of floats
List the azimuthal angles between which pores are generated. Units are in radians.

Note: use [0] for all pores to generate vertically.
### Azimuthal Angle Weighting {`weighting.phi_probs`}
###### Type: list of floats that sum to 1
Chance for azimuthal angles to generate in a given range. The Nth value of this list is the chance for a given pore to generate at an angle between the Nth and (N+1)th value of `weighting.phi_values`.
### Polar Angles {`weighting.theta_values`}
###### Type: list of floats
List the radial angles between which pores are generated. Units are in radians.
### Polar Angle Weighting {`weighting.theta_probs`}
###### Type: list of floats that sum to 1
Chance for radial angles to generate in a given range. The Nth value of this list is the chance for a given pore to generate at an angle between the Nth and (N+1)th value of `weighting.theta_values`.
### Shape Proportions {`params.shape_proportions`}
###### Type: length-5 list of floats that sum to 1
Used with `option.variedPoreShape` to supply the chances of each pore shape generating. The values of this array are as follows: cylinder, proximal opening cone, distal opening cone, ellipsoid, hyperboloid.
### Pores Before Networking {`params.pores_before_networking`}
###### Type: int
The number of pores to generate with random x and y coordinates before the `networkPore` function begins being used to form a networked structure.
### Branches (Top) {`params.top_branches`}
###### Type: length-2 list of ints
Range of number of pores that could branch off the top of each pore in a networked structure.
### Branches (Bottom) {`params.bottom_branches`}
###### Type: length-2 list of ints
Range of number of pores that could branch off the bottom of each pore in a networked structure.
### Sealed Osteon Chance {`params.sealed_osteon_chance`}
###### Type: float
The chance for a given pore to generate with random x and y coordinates after the `params.pores_before_networking` threshold has been reached.
### Transverse Pore Azimuth Limit {`params.transverse_flag_onset`}
###### Type: float
The maximum azimuthal angle after which a pores length in the x and y directions is limited. Units are in radians.
### Target Porosity {`target_porosity`}
###### Type: float
The porosity at which the image generates, represented as a percent of voxels with a value of 0 (porous in the image).
