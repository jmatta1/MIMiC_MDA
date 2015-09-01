"""This file contains the configuration dictionary. This dictionary stores all
the parameters needed for the code to run using names that are human readable
and is intended to be user editable"""

# This is the dictionary that stores configuration information
CONFIGURATION = {}

# path to the shared library that contains the compiled routines for fast calc
# of chi^2, log likelihood, and prior distributions
CONFIGURATION["Shared Lib Path"] = "./C/libChiSq.so"

# Information about inputs to the program
# path to the file containing the experimental data
CONFIGURATION["Input File Path"] = "./58Ni_cefatz_en_row_decomp.csv"
# path to the directory that contains the dwba GR distributions for 100% EWSR
CONFIGURATION["Distribution Directory"] = "./dist/"


# Information about the target nucleus
# A of the target nucleus
CONFIGURATION["Target A"] = 58
# Sets if IVGDR subtraction is carried carried out.
CONFIGURATION["Subtract IVGDR"] = True
# Integral across all ex energy of the IVGDR lorientzian, doesn't matter if
# subtraction is off
CONFIGURATION["IVGDR CS Integral"] = 294.0
# Height of the IVGDR lorentzian in millibarns
CONFIGURATION["IVGDR Height"] = 23.6484
# centroid energy (in MeV) of the IVGDR lorentzian
CONFIGURATION["IVGDR Center"] = 19.0984
# width of the IVGDR lorentzian in MeV
CONFIGURATION["IVGDR Width"] = 7.91453

# Limits on and paramters of the fit
# Maximum L value to fit with
CONFIGURATION["Maximum L"] = 7
# Maximum Angle (in degrees) to fit with
CONFIGURATION["Max Theta"] = 10.0
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
