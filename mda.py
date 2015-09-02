#!/usr/bin/python
"""This program performs a multipole decomposition analysis of the data and
options given in the configuration file mda_config.py and using the Markov
Chain Monte Carlo method to sample and get error bars and best fits"""
import sys
import multiprocessing
import math
import numpy as np
import copy
from scipy import interpolate

# prevent bytecode generation for the config file
ORIGINAL_SYS_DONT_WRITE_BYTECODE = sys.dont_write_bytecode
sys.dont_write_bytecode = True
# import the configuration and template variables
from mda_config import CONFIGURATION as config
# restore the dont_write_bytecode variable to its original value
sys.dont_write_bytecode = ORIGINAL_SYS_DONT_WRITE_BYTECODE


def main():
    """performs sanity checks on the configuration data and then calls the
    functions that do the work"""
    # check that config parameters have sane values
    # first check that they are not requesting greater concurrency than the
    # system supports
    cpu_count = multiprocessing.cpu_count()
    if config["Number of Threads"] > cpu_count:
        print "\nInvalid number of threads, on this machine it must be: "
        print 'config["Number of Threads"] <= %d\n' % cpu_count
        sys.exit()
    num_samples = ((config["Sample Points"] - 50)
                   * config["Number of Walkers"])
    samples_needed = int(math.ceil(10.0 /
                                   ((1.0 - config["Confidence Interval"]) /
                                    2.0)))
    if samples_needed > num_samples:
        print SAMPLES_ERROR % (num_samples, config["Sample Points"],
                               config["Number of Walkers"], samples_needed,
                               config["Confidence Interval"])
        sys.exit()
    num_dists = (1 + config["Maximum L"])
    num_ewsr = len(config["EWSR Fractions"])
    if num_ewsr > num_dists:
        print "\nToo many EWSR fractions listed, there must be",\
            "(1 + CONFIGURATION[\"Maximum L\"]) = %d\n" % num_dists,\
            "EWSR fractions listed in the variable",\
            "CONFIGURATION[\"EWSR Fractions\"], %d were given\n" % num_ewsr
        sys.exit()
    elif num_ewsr < num_dists:
        print "\nToo few EWSR fractions listed, there must be",\
            "(1 + CONFIGURATION[\"Maximum L\"]) = %d\n" % num_dists,\
            "EWSR fractions listed in the variable",\
            "CONFIGURATION[\"EWSR Fractions\"], %d were given\n" % num_ewsr
        sys.exit()
    # call the function that calls everything else
    initialize_mda()


def initialize_mda():
    """does the work of the program, reads in all the data and distributions
    subtracts all the ivgdr components if needed then calls the functions
    that do the initial fitting and then the sampling"""
    # read the raw data
    exp_data = read_row_cs_data_file(config["Input File Path"],
                                     config["Max Theta"],
                                     config["Start Energy"],
                                     config["Final Energy"])
    ivgdr_dists, ivgdr_ewsr, sub_data = handle_ivgdr(exp_data)
    print ivgdr_dists[0], '\n', ivgdr_ewsr[0], '\n', sub_data[0]


def handle_ivgdr(data):
    """This function handles the ivgdr subtraction if no subtraction is needed
    then this function returns None, None, and a deep copy of the exp data
    that was passed to it, otherwise it returns a list of distributions,
    ewsr fractions and the experimental data minus the interpolated IVGDR data
    times the ewsr fraction"""
    # if needed handle the ivgdr
    if config["Subtract IVGDR"]:
        frac = IVGDRFraction(config["IVGDR Height"], config["IVGDR Center"],
                             config["IVGDR Width"],
                             config["IVGDR CS Integral"])
        # calculate the list of IVGDR fractions
        ewsrs = [frac.get_ivgdr_fraction(elem[0]) for elem in data]
        # read the distributions
        dists = read_ivgdr_dists([elem[0] for elem in data])
        # extract the list of angles for each energy
        angle_lists = [ en_set[1][:,0] for en_set in data]
        # calculate scaled and interpolated points
        ivgdr_values = [ interpolate_and_scale_ivgdr(dist, ewsr, angle_list)
                        for (dist, ewsr, angle_list) in zip(dists, ewsrs,
                                                          angle_lists)]
        sub_data = copy.deepcopy(data)
        for i in range(len(sub_data)):
            sub_data[i][1][:,1] = sub_data[i][1][:,1] - ivgdr_values[i]                           
        print "IVGDR subtraction is done"
        return dists, ewsrs, sub_data
    else:
        print "IVGDR is turned off"
        return None, None, copy.deepcopy(data)


def interpolate_and_scale_ivgdr(dist, ewsr, angle_list):
    """This function takes an ivgdr distribution, interpolates it, calculates
    the values of the distribution at the provided angles, multiplies them by
    ewsr and returns it"""
    interp = interpolate.interp1d(dist[:,0],dist[:,1],kind="cubic")
    values = interp(angle_list)
    return (ewsr*values)


def read_ivgdr_dists(en_list):
    """reads in the csv files with the IVGDR distributions"""
    dist_file_names = \
        ["{0:s}A{1:d}_Ex{2:4.2f}_L01_T1_F100.csv".format(
            config["Distribution Directory"], config["Target A"], energy)
         for energy in en_list]
    dists_list = []
    for name in dist_file_names:
        dist_file = open(name, 'r')
        dist = [ [float(elem) for elem in line.strip().split(',')]
                for line in dist_file]
        dists_list.append(np.array(dist, dtype=np.float32))
    return dists_list
            


# This class allows the calculation of IVGDR fractions easily
class IVGDRFraction(object):
    """Object to hold the information needed to calculate a lorentzian divided
    by its integral"""
    # pythonic constructor!
    def __init__(self, max_sigma, centroid, width, total_sigma):
        # to simplify the calculation of the lorentzian we dont store max_sigma
        # and the total_sigma instead we simply store the ratio because that is
        # all that is needed to compute the IVGDR sumrule fraction
        self.sig_ratio = (max_sigma / total_sigma)
        self.total = total_sigma
        # to simplify the calculation of the lorentzian we dont store the
        # centroid and the width instead we store a couple values that result
        # from the expansion of 1+((en^2-centroid^2)^2/(en^2*width^2)) this
        # expansion becomes
        # 1-2*(centroid/width)^2+centroid^4/(en^2*width^2)+en^2/width^2
        # happily 1-2*(centroid/width)^2 is a constant
        self.denom_const = (1.0 - 2.0 * math.pow(centroid / width, 2.0))
        self.width_sq = math.pow(width, 2.0)
        self.denom_inv = math.pow(centroid, 4.0)/self.width_sq

    def get_ivgdr_fraction(self, excitation_energy):
        """returns the IVGDR EWSR % at energy = excitation_energy"""
        en_sq = math.pow(excitation_energy, 2.0)
        denom = self.denom_const + self.denom_inv/en_sq + en_sq/self.width_sq
        return self.sig_ratio/denom

    def get_ivgdr_cs(self, excitation_energy):
        """returns the IVGDR total cs at energy = excitation_energy"""
        en_sq = math.pow(excitation_energy, 2.0)
        denom = self.denom_const + self.denom_inv/en_sq + en_sq/self.width_sq
        return (self.total * self.sig_ratio) / denom


def read_row_cs_data_file(file_name, max_angle, min_en, max_en):
    """takes a cross-section file in row format and extracts the data,
    points with angle above max_angle are excluded"""
    # open the csv file
    input_file = open(file_name, 'r')
    output = []
    # read the file line by line, each line has an energy and a list of angles
    # and cross-sections and errors
    for line in input_file:
        # break the line by the commas
        data_list = line.split(',')
        # get the energy for the distribution
        energy = float(data_list[0].strip())
        # check if the energy is in the acceptable range
        if energy < max_en and energy > min_en:
            # convert each cell to a float
            distribution_data = [float(x.strip()) for x in data_list[1:]]
            # distibution_data = map(lambda x: float(x.strip()), data_list[1:])
            dist_len = len(distribution_data)
            distribution = []
            i = 0
            # for each trio of elements in the distribution put them in a tuple
            # and put that tuple in the distribution list this makes the
            # distribution variable contain a list of points and errors:
            # [(x1,y1,dy1),(x2,y2,dy2),...]
            for i in xrange(0, dist_len, 3):
                if distribution_data[i] <= max_angle:
                    distribution.append((distribution_data[i],
                                         distribution_data[i+1],
                                         distribution_data[i+2]))
                # if 15.4 < energy and energy < 15.6:
                #    print distData[i],",",distData[i+1],",",distData[i+2]
            # put the energy and its associated distribution in a list
            output.append([energy, np.array(distribution, dtype=np.float32)])
    return output


SAMPLES_ERROR = """
WARNING: the number of samples:
%d = ((%d - 50) * %d)
is not large enough to ensure 10 points outside of each error bar location
for the given confidence interval. The minimum number of samples necessary is:
%d = ceiling[ 10 / ((1.0 - %12.10f) / 2.0))

For more information look at the mda_config.py file.
"""


if __name__ == "__main__":
    main()
