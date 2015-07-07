# Strong Gravitational Lens Model & Detect

## Authors
Sean W. McLaughlin  
Robert J. Brunner

###sgl-model-detect is a python program that models light profiles of galaxies using a Mixture-of-Gaussians model.
It is built to fit these models with MCMC using the emcee package. However, it does have the capability to use scipy's optimize to do a non-linear least squares fit as well. It reads in either a single .fits image or a directory of images. The most important modules in the package are:
* modelImages: main module. Imports the chosen fitting module, runs the modeler several times and determines the best model.
* mcmcFit: contains all functions necessary for an MCMC model, using emcee
* nlsqFit: contains all functions necessary for a non-linear least squares model, using scipy's curve_fit. 

In addition to it's modeling modules, this distribution contains a few helper functions. 
* goodnessOfFit: calculates the Chi2Dof of the calculated model
* imageClass: Contains the class for the imageObject used in the modeler. Includes subclasses for different formats.
* makeCenterFile: given a directory of .fits images and a collection of (ra,dec)'s for objects in them, creates a file that has the id's and pixel coordinates of those objects, which the modeler takes as input.
* makeToyImages: Generates a collection of artificial galaxy images, to see if the modeler can recover their parameters.
* residualID: An optional additional step to the modeler. This module has a rudementary strong lens detection pipeline, if the modeler is being applied to a strong lens candidate. 

###Installation
The current implementation does not support standard python installation. The code can be cloned from github via:

`git clone https://github.com/ProfessorBrunner/sgl-model-detect.git`

And run via:

`python modelImages.py <args>`

###Dependencies
This package has the following dependencies:
* numpy
* scipy
* emcee
* pyfits
* astropy

For strong lens identification, an additional requirement is
* scikit-learn

Though not required, the following are reccomended to support plotting functionality.
* matplotlib
* seaborn

All above packages can be installed with the pip package manager, via the following command:

`[sudo] pip install <package>`

###Documentation
An IPython notebook detailing use of this module is included in the **doc** folder. It can either be downloaded and viewed, or viewed directly on GitHub. 
