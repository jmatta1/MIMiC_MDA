"""This file contains the configuration dictionary. This dictionary stores all
the parameters needed for the code to run using names that are human readable
and is intended to be user editable"""

# This is the dictionary that stores configuration information
# Below are all the user modifiable configuration options
CONFIGURATION = {}

# path to the shared library that contains the compiled routines for fast calc
# of chi^2, log likelihood, and prior distributions
CONFIGURATION["Shared Lib Path"] = "./C/libChiSq.so"

# Information about inputs to the program
# path to the file containing the experimental data
CONFIGURATION["Input File Path"] = "./58Ni_cefatz_en_row_decomp.csv"  # used
# path to the directory that contains the dwba GR distributions for 100% EWSR
CONFIGURATION["Distribution Directory"] = "./dist/"

# Information about outputs from the program
# file to hold the fitted parameters
CONFIGURATION["Parameter File"] = "./output.csv"
# path to the directory that will hold the distribution outputs
CONFIGURATION["Fits Csv Directory"] = "./fits/csv/"
# directory to store triangle plots
CONFIGURATION["Triangle Plots Directory"] = "./fits/triangles/"
# directory to store fit plots
CONFIGURATION["Fit Plots Directory"] = "./fits/fit_plots/"
# directory to store parameter plots
CONFIGURATION["Parameter Plots Directory"] = "./fits/param_plots/"


# Information about the target nucleus
# A of the target nucleus
CONFIGURATION["Target A"] = 58
# Sets if IVGDR subtraction is carried carried out.
CONFIGURATION["Subtract IVGDR"] = True
# Integral across all ex energy of the IVGDR lorientzian, doesn't matter if
# subtraction is off
CONFIGURATION["IVGDR CS Integral"] = 294.0  # used
# Height of the IVGDR lorentzian in millibarns
CONFIGURATION["IVGDR Height"] = 23.6484  # used
# centroid energy (in MeV) of the IVGDR lorentzian
CONFIGURATION["IVGDR Center"] = 19.0984  # used
# width of the IVGDR lorentzian in MeV
CONFIGURATION["IVGDR Width"] = 7.91453  # used

# Information setting what energies and angles are fit
# Maximum Angle (in degrees) to fit with
CONFIGURATION["Max Theta"] = 10.0  # used
# Energy such that all fitted energies are >= to it
CONFIGURATION["Start Energy"] = 11.4  # used
# Energy such that all fitted energies are <= to it
CONFIGURATION["Final Energy"] = 31.6  # used

# Limits on and paramters of the fit
# Maximum L value to fit with
CONFIGURATION["Maximum L"] = 7
# List of sum rule fractions (on the interval [0, 1])
CONFIGURATION["EWSR Fractions"] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# Confidence interval of the error bars
# ((CONFIGURATION["Sample Points"] - 50) * CONFIGURATION["Number of Walkers"])
# must exceed, by a large margin, the following value:
# 1 / ((1 - CONFIGURATION["Confidence Interval"]) / 2)
# 1 sigma -> 0.682689492
# 2 sigma -> 0.954499736
# 3 sigma -> 0.997300204
# 4 sigma -> 0.99993666
CONFIGURATION["Confidence Interval"] = 0.682689492
# holds the number of points to sample for each walker
CONFIGURATION["Sample Points"] = 500
# holds the number of walks to run
CONFIGURATION["Number of Walkers"] = 1000
# holds the number of concurrent threads
CONFIGURATION["Number of Threads"] = 4
