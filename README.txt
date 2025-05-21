BLS_Generate README

SUMMARY

BLS_Generate is a Python code used to generate broadband light scattering (BLS) spectra for spherical homogeneous particles. Full details of the theoretical framework used for these calculations can be found in the associated manuscript [insert doi when available]. A brief summary is given below.

Users input 5 parameters:

wl - the wavelengths at which the scattering intensity is to be evaluated
r - the particle radius values to use for calculations
theta_i, phi_i - together these specify the direction of the incident light. See publication for the geometry assumed
NA - the numerical aperture over which the scattered light is collected

Related to r and wl is the refractive index. Typically this is dependent on the wavelength, although the exact form differs depending on the desired simulation. Examples will be added in time to illustrate common situations of interest. 

Using the input parameters, the code calculates the array spectra where each row is the scattering intensity at each input wl for a given value of r. 

REQUIREMENTS

In order to run this code, you must have Python 3 installed on your computer. Those new to Python are recommended to install Anaconda (https://www.anaconda.com/download) as this is, in the author's experience, the easiest way to get started with Python. Code was written in Pyhton 3.11.5, but should runwith any version of Python 3.

In the main folder one will find the code in two formats: .py and .ipynb. The .py file is intended for those who simply wish to run the code. The .ipynb contains additional explanatory notes and plots spectra at the end of the Notebook (not done in the .py file). Only two packages are required: numpy and time. Numpy should be part of any Python IDE, but installation instructions for numpy can be found at https://numpy.org/install/. The time package is part of the Python Standard Library (https://docs.python.org/3/library/time.html). The final import required is BLSfunctions.py, which is found in the main folder alongside BLS_Generate.py and BLS_Generate.ipynb. BLSfunctions contains functions for evaluating parameters using Mie theory and should be downloaded and included in the same directory that you run BLS_Generate from.

AUTHORS

This code is written and maintained by Aidan Rafferty. He can be contacted at aidan.rafferty@chem.ox.ac.uk (until end of September 2025, will update then). 

COPYRIGHT AND LICENSING

This work is licensed under the CC BY 4.0 License (more information in LICENSE.txt or at https://creativecommons.org/licenses/by/4.0/).