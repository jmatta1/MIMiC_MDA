"""This file contains the configuration dictionary. This dictionary stores all
the parameters needed for the code to run using names that are human readable
and is intended to be user editable"""

# This is the dictionary that stores configuration information
# Below are all the user modifiable configuration options
CONFIG = {}

# path to the shared library that contains the compiled routines for fast calc
# of chi^2, log likelihood, and prior distributions
CONFIG["Shared Lib Path"] = "./C/libChiSq.so"  # complete

# Information about inputs to the program
# path to the file containing the experimental data
CONFIG["Input File Path"] = "./58Ni_cefatz_en_row_decomp.csv"  # complete
# path to the directory that contains the dwba GR distributions for 100% EWSR
CONFIG["Distribution Directory"] = "./dist/"  # complete

# Information about outputs from the program
# file to hold the fitted parameters
CONFIG["Parameter File"] = "./fits/params.csv"
# path to the directory that will hold the distribution outputs
CONFIG["Fits Csv Directory"] = "./fits/csv/"
# directory to store triangle plots this is used by the seperate, make
# triangles code due to how much memory this consumes
CONFIG["Corner Plots Directory"] = "./fits/corner_plots/"
# directory to store probability plots
CONFIG["Prob Plots Directory"] = "./fits/prob_plots/"
# directory to store fit plots
CONFIG["Fit Plots Directory"] = "./fits/fit_plots/"
# directory to store sampled Markov Chains (in the compressed numpy format)
# this is used by teh seperate make triangles code as a data source
CONFIG["Chain Directory"] = "./fits/chains/"  # complete
# directory to store parameter plots
CONFIG["Parameter Plots Directory"] = "./fits/param_plots/"
# toggle the output of corner plots on or off, if the number of sample points
# (see notes below) is large this output should be turned off, later they can
# be generated by a seperate script, for ~1000000 sample points making this
# plot consumes ~840MB, for ~8000000 points it consumes ~4.8GB
CONFIG["Generate Corner Plots"] = True
# sets the number of samples for corner plots if you have a lot of samples but
# wish to generate the corner plot without the memory penalty you can set this
# to use a random selection of however many samples out of the total
CONFIG["Corner Plot Samples"] = 1000000
# Toggle the output of sample chains, if you disable the output of Corner plots
# this needs to be on if you later want to make the corner plots in the
# seperate script if you want to analyze the samples in other ways
CONFIG["Save Chain Data"] = True

# Information about the target nucleus
# A of the target nucleus
CONFIG["Target A"] = 58  # complete
# Sets if IVGDR subtraction is carried carried out.
CONFIG["Subtract IVGDR"] = True  # complete
# Integral across all ex energy of the IVGDR lorientzian, doesn't matter if
# subtraction is off
CONFIG["IVGDR CS Integral"] = 294.0  # complete
# Height of the IVGDR lorentzian in millibarns
CONFIG["IVGDR Height"] = 23.6484  # complete
# centroid energy (in MeV) of the IVGDR lorentzian
CONFIG["IVGDR Center"] = 19.0984  # complete
# width of the IVGDR lorentzian in MeV
CONFIG["IVGDR Width"] = 7.91453  # complete

# Information setting what energies and angles are fit
# Maximum Angle (in degrees) to fit with
CONFIG["Max Theta"] = 10.0  # complete
# Energy such that all fitted energies are >= to it
CONFIG["Start Energy"] = 11.4  # complete
# Energy such that all fitted energies are <= to it
CONFIG["Final Energy"] = 31.6  # complete

# Limits on and paramters of the fit
# Maximum L value to fit with
CONFIG["Maximum L"] = 7  # incomplete
# List of sum rule fractions (on the interval [0, 1])
CONFIG["EWSR Fractions"] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # complete

# Confidence interval of the error bars
# number of sample points (see notes below) must exceed, by a fair margin,
# the following value: 1 / ((1 - CONFIG["Confidence Interval"]) / 2)
# 1 sigma -> 0.682689492
# 2 sigma -> 0.954499736
# 3 sigma -> 0.997300204
# 4 sigma -> 0.99993666
CONFIG["Confidence Interval"] = 0.682689492  # complete
# holds the number of points to sample for each walker
CONFIG["Sample Points"] = 2050  # complete
# holds the number of walks to run
CONFIG["Number of Walkers"] = 4000  # complete
# holds the number of concurrent threads
CONFIG["Number of Threads"] = 4  # complete
# Lists of starting points for each a_L in the initial fits
# (max a_L = 1.0/(EWSR Fraction_L)
# ... This continues for as many Ls as you have
# for example: CONFIG["Start Pts a23"] = [0.5]
CONFIG["Start Pts a0"] = [0.0, 0.05, 0.10, 0.15, 0.2]  # complete
CONFIG["Start Pts a1"] = [0.0, 0.05, 0.10, 0.15, 0.2]  # complete
CONFIG["Start Pts a2"] = [0.0, 0.05, 0.10, 0.15, 0.2]  # complete
CONFIG["Start Pts a3"] = [0.0, 0.10, 0.2]  # complete
CONFIG["Start Pts a4"] = [0.05, 0.15]  # complete
CONFIG["Start Pts a5"] = [0.05, 0.15]  # complete
CONFIG["Start Pts a6"] = [0.1]  # complete
CONFIG["Start Pts a7"] = [0.1]  # complete
# Number of refined points to use to generate initial positions for walkers
CONFIG["Number Walker Generators"] = 200  # complete
# Spread of the initial sampling positions for walkers
CONFIG["Sample Spread"] = 0.15  # complete
# Number of points at the beginning of the Markov Chain to discard as 'burn in'
CONFIG["Burn-in Points"] = 50  # incomplete
# NOTES:
# Note about how to calculate the number of sample points
# the total number of sample points is the product of the following two values
# number of useable sample points per walker:
# (CONFIG["Sample Points"]-CONFIG["Burn-in Points"])
# and the number of walkers:
# CONFIG["Number of Walkers"])
