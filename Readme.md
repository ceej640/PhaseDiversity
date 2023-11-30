# Phase Diversity Wavefront Sensing for Microscopy

This package is the companion code to our paper, Phase diversity-based wavefront sensing for fluorescence microscopy.

Given an image contaminated with an unknown aberration and additional images contaminated with the same unknown aberration and additional known aberrations, phase diversity can estimate the wavefront associated with the unknown aberration. This package contains MATLAB code for 2 implementations of Phase Diversity and their dependency functions.

1) ```script_calibration.m``` - Used for Phase Diversity-based calibration of a deformable mirror. Takes defocus images stacks from a series of voltages applied to each deformable mirror (DM) actuator, estimates wavefronts and computes influence functions. Outputs command matrix (mirror voltages required to issue a desired wavefront).

2) ```processPhaseDiversityImages.m``` - function designed for on-line (or off-line) calls to estimate wavefronts from an image stack of diversity aberrations applied with a DM. Returns object and wavefront estimate, which can be used to apply optical correction with the DM.


## Table of Contents

- Requirements
- Installation
- Datasets
- Usage
- References
- License


## Requirements

The code in this package requires MATLAB 23.2 (2023b) with the following toolboxes:

- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox
- Parallel Computing Toolbox

The MATLAB function ```dftregistration.m``` is also required [1], but included in the folder ```PhaseDiversity/MATLAB code/redist```.

A GPU is required to use the ```flagGPU``` option. 
This code has been tested and verified on Intel Core i7, Intel Xeon, and AMD Ryzen Threadripper CPUs, and Nvidia GeForce and RTX series GPUs.


## Installation

To use this code, clone the repository or download and unzip it into a suitable location. Add the folder ```PhaseDiversity``` and its subfolders to the MATLAB path.

## Datasets

### Deformable Mirror Calibration
To use the test datasets, download from [figshare](https://figshare.com/s/774457a804ee0686930d) and unzip the contents so that datasets is a subfolder of the folder ```PhaseDiversity```.

3 full sets of DM calibration data are provided in the data subfolder ```Datasets/DM Calibration```. This folder contains a subfolder for each dataset (e.g.: ```PD_Cal_1```). Each dataset was acquired on a different field of view. Each calibration folder contains a series of subfolders in the form ```A_n```, where ```n``` is the number corresponding to each DM actuator. The actuator subfolder then contains an image stack for each of two voltages in the form ```A_Poke_N.tif``` where ```N``` is the voltage value.

The command matrix can be output using ```script_calibration.m``` and pointing to the ```DM Calibration``` folder located in ```PhaseDiversity/Datasets```. Figure outputs and MAT data of wavefront estimates, influence functions, and intermediate results of the algorithm can be output by changing the SaveScriptFigures, SaveFunctionFigures and SaveMAT to true. Note that the runtime will slow considerably when figure saving is enabled.

The calibration data was acquired under the following system:

- ```pixelSize``` = 0.104 (Image Pixel Size in µm)
- ```lambda``` = 0.532 (Central Emission Wavelength in µm)
- ```NA``` = 1.2 (Numerical Aperture)
- ```numz``` = 5 (Number of slices in defocus image stack corresponding to the raw image + number of diversity images)
- ```zstepsize``` = 1 (Space between defocus slices in µm)


### AO Correction

Two example sets of AO data are provided in the folder ```PhaseDiversity/Datasets/AO```: an example with multiple beads featuring a scaled test aberration, and a U2OS cell sample with a randomly generated test aberration. Each dataset's folder contains a series of subfolders and two files: GT.tif and CoeffsIn.txt

```GT.tif``` - Ground Truth image before test aberration is applied. For reference only, not used in computation.

```CoeffsIn.txt``` - Table of Zernike Coefficients. Columns correspond to Zernike modes 1-20 (ANSI notation). Rows 2-5 correspond to the known Zernike coefficients for diversity phases. Row 1 is the Zernike coefficients corresponding to the test aberration that are present in all stack images before correction.

Each subfolder is labeled ```Iteration n``` where iteration refers to the cycle number/number of wavefront measurements and n corresponds to correction n-1. For example, the subfolder ```Iteration 1``` corresponds to the uncorrected aberrated image and ```Iteration 2``` corresponds to the first correction. In each iteration subfolder is a single file: ```Stack.tif```, which is an image stack corresponding to the raw image and four diversity phases. The first slice in each image stack corresponds to the raw image. In ```Iteration 1```, this slice corresponds to the aberrated image and in ```Iterations 2...N``` the first slice corresponds to correction n-1.

Each data set has an associated wrapper function, ```runPD_TestBeadData``` and ```runPD_TestCellData```, respectively, which are used to replicate the original on-line calls from LabVIEW. To use, simply call the function with the desired cycle number as input. The function will locate the data folder if PhaseDiversity is found in the MATLAB path and call processPhaseDiversityImages using the original set of input parameters. This function will output the estimated Zernike coefficients, wavefront, and object estimate. It will also return the time to estimate the wavefront, and the difference between the measured wavefront and the expected wavefront for the test aberration.

The input number for the wrapper function is the desired correction cycle equal to n-1 corrections. For example, a value of 1 corresponds to the uncorrected aberrated wavefront and 2 corresponds to the first correction. 10 is the maximum cycle number for the bead test data and 8 is the maximum cycle number for the cell data.



## References

[1] Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).


## License

Copyright © 2023 Howard Hughes Medical Institute

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of HHMI nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
