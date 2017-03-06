"""This file contains the configuration dictionary. This dictionary stores all
the parameters needed for the code to run using names that are human readable
and is intended to be user editable"""

# This is the dictionary that stores configuration information
# Below are all the user modifiable configuration options
CONFIG = {}

# path to the shared library that contains the compiled routines for fast calc
# of chi^2, log likelihood, and prior distributions
CONFIG["Shared Lib Path"] = "./C/libChiSq.so"

# Information about inputs to the program
# path to the file containing the experimental data
CONFIG["Input File Path"] = "./58Ni_tz_fa_ce_mda.csv"
# path to the directory that contains the dwba GR distributions
# the distributions, regardless of sum rule fraction, should all start at 0 deg
# and have identical sampling angles and end points
CONFIG["Distribution Directory"] = "./58Ni_dist/"

# Information about outputs from the program
# file to hold the fitted parameters
CONFIG["Parameter Files Directory"] = "./fits1/"
# path to the directory that will hold the distribution outputs
CONFIG["Fits Csv Directory"] = "./fits1/fit_csv/"
# directory to store triangle plots this is used by the seperate, make
# triangles code due to how much memory this consumes
CONFIG["Corner Plots Directory"] = "./fits1/corner_plots/"
# directory to store probability plots
CONFIG["Prob Plots Directory"] = "./fits1/prob_plots/"
# directory to store fit plots
CONFIG["Fit Plots Directory"] = "./fits1/fit_plots/"
# directory to store sampled Markov Chains (in the compressed numpy format)
# this can then be used for seperate analysis of the sampling
CONFIG["Chain Directory"] = "./fits1/chains/"
# directory to store parameter plots
CONFIG["Parameter Plots Directory"] = "./fits1/param_plots/"
# directory to store the autocorrellation plot sub directories
CONFIG["Time Series Directory"] = "./fits1/walker_plots/"
# format for plot outputs. options are:
# "svg" - scalable vector graphics
# "svgz" - scalable vector graphics (compressed)
# "pdf" - portable document format
# "ps" - post script
# "eps" - encapsulated post script
# "png" - portable netword graphics
CONFIG["Plot Format"] = "png"
# the height of the plot image in inches
CONFIG["Plot Height"] = 12
# the width of the plot image in inches
CONFIG["Plot Width"] = 9
# the dpi of the plot
CONFIG["Plot DPI"] = 254
# toggle the output of corner plots on or off, if the number of sample points
# (see notes below) is large this output should be turned off, later they can
# be generated by a seperate script, for ~1000000 sample points making this
# plot consumes ~840MB, for ~8000000 points it consumes ~4.8GB
CONFIG["Generate Corner Plots"] = True
# sets the number of samples for corner plots if you have a lot of samples but
# wish to generate the corner plot without the memory penalty you can set this
# to use a random selection of however many samples out of the total
CONFIG["Corner Plot Samples"] = 1000000
# sets the number of bins for each parameter axis in the corner plot
CONFIG["Corner Plot Bins"] = [100, 100, 100, 100, 100, 100, 100, 100]
# allow the corner package to set param ranges instead of [0, param_max]
CONFIG["Corner Default Range"] = True
# sets whether or not the user wants to plot the IVGDR in the unlimited fit
# plots
CONFIG["Plot IVGDR in Fits"] = False
# sets the max L the user wants to plot in the limited fit plot
CONFIG["Fit Plot L Limit"] = 3
# Toggle the output of sample chains, if you disable the output of Corner plots
# this needs to be on if you later want to make the corner plots in the
# seperate script if you want to analyze the samples in other ways
CONFIG["Save Chain Data"] = False
# Number of walkers to plot in the time series plots
CONFIG["Walker Plot Count"] = 4000

# Information about the target nucleus
# A of the target nucleus
CONFIG["Target A"] = 58
# Sets if IVGDR subtraction is carried carried out.
CONFIG["Subtract IVGDR"] = True
# Integral across all ex energy of the IVGDR lorientzian, doesn't matter if
# subtraction is off
CONFIG["IVGDR CS Integral"] = 294.0
# Height of the IVGDR lorentzian in millibarns
CONFIG["IVGDR Height"] = 23.6484
# centroid energy (in MeV) of the IVGDR lorentzian
CONFIG["IVGDR Center"] = 19.0984
# width of the IVGDR lorentzian in MeV
CONFIG["IVGDR Width"] = 7.91453

# Information setting what energies and angles are fit
# Maximum Angle (in degrees) to fit with
CONFIG["Max Theta"] = 10.0
# Energy such that all fitted energies are >= to it
CONFIG["Start Energy"] = 8.4
# Energy such that all fitted energies are <= to it
CONFIG["Final Energy"] = 31.6

# Limits on and paramters of the fit
# Maximum L value to fit with
CONFIG["Maximum L"] = 7
# List of sum rule fractions (on the interval [0, 1])
CONFIG["EWSR Fractions"] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# Confidence interval of the error bars
# number of sample points (see notes below) must exceed, by a fair margin,
# the following value: 1 / ((1 - CONFIG["Confidence Interval"]) / 2)
# 1 sigma -> 0.682689492
# 2 sigma -> 0.954499736
# 3 sigma -> 0.997300204
# 4 sigma -> 0.99993666
CONFIG["Confidence Interval"] = 0.682689492
# holds the number of points to sample for each walker, should be about 10
# autocorrellation times, 1 auto correllation time is the time it takes for the
# walkers for all parameters to drift into their designated sampling region
CONFIG["Sample Points"] = 3000
# holds the number of walks to run
CONFIG["Number of Walkers"] = 4000
# holds the number of concurrent threads
CONFIG["Number of Threads"] = 1
# holds the number of bins to use when determining the peak of the probability
# density distribution
CONFIG["Num Bins"] = 1000
# holds the difference required to make floating point comparisons evaluate
# as equal values
CONFIG["Float Epsilon"] = 1.0e-7
# Lists of starting points for each a_L in the initial fits
# (max a_L = 1.0/(EWSR Fraction_L)
# ... This continues for as many Ls as you have
# for example: CONFIG["Start Pts a23"] = [0.5]
CONFIG["Start Pts a0"] = [0.0, 0.05, 0.10, 0.15, 0.2]
CONFIG["Start Pts a1"] = [0.0, 0.05, 0.10, 0.15, 0.2]
CONFIG["Start Pts a2"] = [0.0, 0.05, 0.10, 0.15, 0.2]
CONFIG["Start Pts a3"] = [0.05, 0.10, 0.15]
CONFIG["Start Pts a4"] = [0.05, 0.15]
CONFIG["Start Pts a5"] = [0.05, 0.15]
CONFIG["Start Pts a6"] = [0.05]
CONFIG["Start Pts a7"] = [0.05]
# number of times to force bfgs to rerun for the initial fits
CONFIG["Forced Extra Fits"] = 3
# Number of refined points to use to generate initial positions for walkers
CONFIG["Number Walker Generators"] = 30
# Spread of the initial sampling positions for walkers
CONFIG["Sample Spread"] = 0.75
# Centroid of the offset applied to sampling positions for walkers
CONFIG["Sample Offset Centroid"] = 0.02
# Width of the offset applied to the sampling positions for walkers
CONFIG["Sample Offset Width"] = 0.01
# Number of points at the beginning of the Markov Chain to discard as 'burn in'
# should be about one auto-correllation time
CONFIG["Burn-in Points"] = 300
# NOTES:
# Note about how to calculate the number of sample points
# the total number of sample points is the product of the following two values
# number of useable sample points per walker:
# (CONFIG["Sample Points"]-CONFIG["Burn-in Points"])
# and the number of walkers:
# CONFIG["Number of Walkers"])
