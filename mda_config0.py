"""This file contains the configuration dictionary. This dictionary stores all
the parameters needed for the code to run using names that are human readable
and is intended to be user editable"""

# This is the dictionary that stores configuration information
# Below are all the user modifiable configuration options
# aside from this line needing to be at the top, all the other options are
# order agnostic
CONFIG = {}

###############################################################################
#  MIMiC_MDA Calculation Library
###############################################################################
# path to the shared library that contains the compiled routines for fast calc
# of chi^2, log likelihood, and prior distributions
CONFIG["Shared Lib Path"] = "./C/libChiSq.so"

###############################################################################
#  Input Files and configuration (exp data and calculated distributions)
###############################################################################
# Information about inputs to the program
# path to the file containing the experimental data
CONFIG["Input File Path"] = "./58Ni_tz_fa_ce_mda.csv"
# path to the directory that contains the dwba GR distributions
# the distributions, regardless of sum rule fraction, should all start at 0 deg
# and have identical sampling angles and end points, the file names in the
# folder must have the naming scheme:
# "A{NucA:d}_Ex{ExEn:4.2f}_L{OrbAng:02d}_T{Iso:d}_F{Frac:03d}.csv"
# where:
# NucA is the nucleus's total number of nucleons
# ExEn is the excitation energy which is always written with 2 decimal places
# OrbAng is the two digit orbital angular momentum, padded with a '0' if needed
# Iso is the isospin, which should always be zero except for the L01_T1 IVGDR
# Frac is the sum rule fraction multiplied by 100 and converted to an integer
CONFIG["Distribution Directory"] = "./58Ni_dist/"
# A of the target nucleus, for reading input and writing output
CONFIG["Target A"] = 58
# List of sum rule fractions (on the interval [0, 1]), this can be used if you
# want to perform the fit with varying sum rule fractions for each distribution
# essentially, this affects the parameter limits and it affects which
# distributions are read in for each L value
CONFIG["EWSR Fractions"] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

###############################################################################
#  Output Directory Structure
###############################################################################
# Information about outputs from the program
# directory to hold the fitted parameters, also holds the diagnostic file
CONFIG["Parameter Files Directory"] = "./fits1/"
# path to the directory that will hold the distribution outputs
CONFIG["Fits Csv Directory"] = "./fits1/fit_csv/"
# directory to store triangle plots this is used by the separate
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

###############################################################################
#  IVGDR Settings and information
###############################################################################
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

###############################################################################
#  Exp data limits and parameter limits
###############################################################################
# Information setting what energies and angles are fit
# Maximum Angle (in degrees) to fit with
CONFIG["Max Theta"] = 10.0
# Energy such that all fitted energies are >= to it
CONFIG["Start Energy"] = 8.4
# Energy such that all fitted energies are <= to it
CONFIG["Final Energy"] = 8.6  # 31.6
# Limits on and paramters of the fit
# Maximum L value to fit with
CONFIG["Maximum L"] = 7

###############################################################################
#  Sampling and Concurrency Configuration
###############################################################################
# holds the number of concurrent threads
CONFIG["Number of Threads"] = 1
# holds the number of points to sample for each walker, should be about 10
# autocorrellation times, 1 auto correllation time is the time it takes for the
# walkers for all parameters to drift into their designated sampling region
CONFIG["Sample Points"] = 16385
# holds the number of walks to run
CONFIG["Number of Walkers"] = 2500
# Number of points at the beginning of the Markov Chain to discard as 'burn in'
# should be one to two of the longest auto-correlation times
# if you get errors like "WARNING:root:Too few points to create valid contours"
# then typicaly you did not set burn-in points high enough because it has some
# points far from the proper region that it is trying to make contours for and
# failing. Try increasing burn-in points and run again
CONFIG["Burn-in Points"] = 512
# Notes about sampling accuracy
# the total number of used sample points can be calculated to be:
# Samples = "Number of Walkers"*("Sample Points"-"Burn-in Points")
# the estimate of the error in the sampled distribution can then be obtained by
# looking the total number of independent samples of the distribution obtained.
# The number of independent samples for each parameter is the used samples
# divided by the autocorrelation time for that parameter
# IndSampP_i = Samples / (Parameter AutoCorr Time)
# We then take the minimum number of independent samples of the distribution
# MinIndSamp = Minimum(IndSampP_0, IndSampP_1, IndSampP_2, ...)
# and from this we can estimate the distribution error to be:
# DistErr ~= 1/Sqrt(MinIndSamp)
# Therefore, with a maximum autocorrelation of 300 and Samples = 19200000
# we get a sampling error of 1/Sqrt(63144) ~= 0.35%
# With a maximum autocorrelation of 300 and Samples = 3000000
# we get a sampling error of 1/Sqrt(10000) ~= 1%


###############################################################################
#  Initial Fits Configuration
###############################################################################
# Lists of starting points for each a_L in the initial fits
# (max a_L = 1.0/(EWSR Fraction_L)
# ... This continues for as many Ls as you have
# for example: CONFIG["Start Pts a23"] = [0.5]
CONFIG["Start Pts a0"] = [0.05, 0.15, 0.25, 0.35, 0.45]
CONFIG["Start Pts a1"] = [0.05, 0.15, 0.25, 0.35, 0.45]
CONFIG["Start Pts a2"] = [0.05, 0.15, 0.25, 0.35, 0.45]
CONFIG["Start Pts a3"] = [0.05, 0.25, 0.45]
CONFIG["Start Pts a4"] = [0.05, 0.2]
CONFIG["Start Pts a5"] = [0.05, 0.2]
CONFIG["Start Pts a6"] = [0.05, 0.2]
CONFIG["Start Pts a7"] = [0.05, 0.2]
# Spread of the initial sampling positions for walkers setting this to be large
# is a good idea, it really forces the system to get a good sampling of the 
# broader distribution which in turn fills out the corner plots nicely (as well
# as giving better results)
CONFIG["Sample Spread"] = 2.0
# Centroid of the offset applied to sampling positions for walkers
CONFIG["Sample Offset Centroid"] = 0.04
# Width of the offset applied to the sampling positions for walkers
CONFIG["Sample Offset Width"] = 0.02
# number of times to force bfgs to rerun for the initial fits
CONFIG["Forced Extra Fits"] = 3
# Number of refined points to use to generate initial positions for walkers
CONFIG["Number Walker Generators"] = 50
# Walker start positions are calculated as follows, select a start position
# from the fit points add a random offset with gaussian distribution,
# the width and centroid of which are specified above, to each parameter. Then
# multiply each parameter by a random value with gaussian distribution centered
# at 0 with width "Sample Spread". If a parameter is negative, negate it, if a
# parameter is 0, make it 0.001

###############################################################################
#  Parameter and Error Finding configuration
###############################################################################
# Confidence interval of the error bars
# number of sample points (see notes below) must exceed, by a fair margin,
# the following value: 1 / ((1 - CONFIG["Confidence Interval"]) / 2)
# 1 sigma -> 0.682689492
# 2 sigma -> 0.954499736
# 3 sigma -> 0.997300204
# 4 sigma -> 0.99993666
CONFIG["Confidence Interval"] = 0.682689492
# holds the number of bins to use when determining the peak of the probability
# density distribution
CONFIG["Num Bins"] = 1000

###############################################################################
#  Optional Output Generation Configuration
###############################################################################
# Toggle the output of sample chains, if some of the plots are turned off, you
# may need to chain data to recreate them later.
CONFIG["Save Chain Data"] = True
# holds if the autocorrelation time should be calculated, its nice to know,
# but requires a lot of extra samples to know accurately
CONFIG["Calc AutoCorr"] = False
# toggle the output of corner plots on or off, if the number of sample points
# (see notes below) is large this output should be turned off, later they can
# be generated by a seperate script, for ~1000000 sample points making this
# plot consumes ~840MB, for ~8000000 points it consumes ~4.8GB
CONFIG["Generate Corner Plots"] = False
# sets whether or not to write out walker plots
# if chains are saved, walker plots could be generated at a later date
CONFIG["Generate Walker Plots"] = False
# toggle the output of the univariate probability plots
CONFIG["Generate Probability Plots"] = False
# toggle the output of fit plots
CONFIG["Generate Fit Plots"] = False
# toggle the output of fit CSVs
CONFIG["Generate Fit CSVs"] = False
# toggle the output of parameter plots
CONFIG["Generate Parameter Plots"] = False


###############################################################################
#  Plot Configuration
###############################################################################
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
# sets the number of samples for corner plots if you have a lot of samples but
# wish to generate the corner plot without the memory penalty you can set this
# to use a random selection of however many samples out of the total
CONFIG["Corner Plot Samples"] = 1000000
# sets the number of bins for each parameter axis in the corner plot
CONFIG["Corner Plot Bins"] = [100, 100, 100, 100, 100, 100, 100, 100]
# allow the corner package to set param ranges instead of [0, 1.1*param_max]
CONFIG["Corner Default Range"] = True
# sets whether or not the user wants to plot the IVGDR in the unlimited fit
# plots
CONFIG["Plot IVGDR in Fits"] = False
# sets the max L the user wants to plot in the limited fit plot
CONFIG["Fit Plot L Limit"] = 3
# Number of walkers to plot in the time series plots
CONFIG["Walker Plot Count"] = 500


###############################################################################
#  Autocorrelation Calculation Accuracy Options
###############################################################################
# Autocorrelation check parameters
# holds the number of autocorrellation times required to estimate the
# autocorrelation time of the chains the larger this is the more accurate the
# estimate
# NOTES: a word of warning on autocorrelation checks: they require either a
# very short window size (CONFIG["ACorr WindSize"] or AC_WIND) which can reduce
# the accurace of the autocorrellation time calculation, or they require a lot
# more than the recommended number of samples. In particular you need to
# fulfill 2*AC_WIND*AC_WIND*(Longest ACorr Time) < CONFIG["Sample Points"]
CONFIG["ACorr WindSize"] = 3.0
# This parameter sets if the effective sample size for the autocorrelation time
# calculation should be restricted to 2^n samples (where n is as large as,
# possible such that 2^n <= CONFIG["Sample Points"]), this does significantly
# speed up acorr time calculation, but it does shrink the effective number of
# pointes available for calculating the acorr time
CONFIG["ACorr Use FFT"] = True
