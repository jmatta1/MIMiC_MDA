# MIMiC MDA
Markov chaIn Monte Carlo Multipole Decomposition Analysis.  A system for giant resonanace Multiple Decomposition Analysis (MDA) using python, the python Markov Chain Monte Carlo (MCMC) package 'emcee' and a dynamicly loaded library written in C to accellerate the calculation of the log-likelihood function needed by the MCMC library.

## Overview
Multipole Decomposition Analysis (MDA) is a technique to extract the amplitudes of the superposition of scattering angular distributions for single states that sum to give the observed angular distribution. While this may seem like a straightforward linear fitting problem, the angular distributions from individual states are 'too similar' or 'not orthogonal enough' for traditional fitting methods to work. MIMiC_MDA solves this issue by using Markov Chain Monte Carlo to extract both best fit parameters and error bars of the distributions directly from the posterior probability distributions of the parameters that it extracted via sampling.

## Input 
### Primary Input File
MIMiC_MDA takes one primary input file on the command line. This input file is a python file that defines a dictionary containing a number of string keys. These keys correspond to configuration parameters, whose function well documented in the file itself with python comments.

### Calculated Angular Distributions
MIMiC_MDA needs calculated angular distributions to decompose the experimental angular distribution into. These files, each containing a single angular distribution are stored in the directory given in the `"Distribution Directory"` key of the primary input file. The file has one angle and cross-section per line, in strictly increasing order of angle. The file name must be given as follows `A{mass_num:d}_Ex{excitation:4.2f}_L{angular_momentum:02d}_T{isovector_change:1d}_F{sum_rule_fraction_percent:03d}.csv` where `mass_num` is the mass number of the nucleus for which these distributions were computed, `excitation` is the excitation energy of the state this was calculated for in MeV, `angular_momentum` is the number of units of orbital angular momentum of the state this was calculated for, `isovector_change` is the units of iso-spin transferred to produce this state (normally 0, only the Iso-Vector Giant Dipole Resonance, IVGDR is included in these decomposition usually), and `sum_rule_fraction_percent` is the sum rule fraction the distribution was calculated with.

### Experimental Data
MIMiC_MDA also needs the set of raw exeperimental data that is to be decomposed. This is in the file whose name is given in the primary input file key `Input File Path`. This file has a csv structure that is as follows for each line.
`ex_energy, angle_1, cross_section_1, cs_error_1, angle_2, cross_section_2, cs_error_2, ...,  angle_N, cross_section_N, cs_error_N`
Where `ex_energy` is the excitation energy of the angular distribution given on that line, `angle_1` is the angle of the first point of the angular distribution, `cross_section_1` is the cross-section at `angle_1`, and `cs_error_1` is the absolute error in `cross_section_1`. This pattern is repeated for all `N` points in the experimental angular distribution at that energy.
