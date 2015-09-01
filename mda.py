#!/usr/bin/python
"""This program performs a multipole decomposition analysis of the data and
options given in the configuration file mda_config.py and using the Markov
Chain Monte Carlo method to sample and get error bars and best fits"""
import sys  # for exit
import multiprocessing  # for cpu_count
import math
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
    # call the function that calls everything else
    initialize_mda()


def initialize_mda():
    """does the work of the program, reads in all the data and distributions
    subtracts all the ivgdr components if needed then calls the functions
    that do the initial fitting and then the sampling"""
    exp_data = read_row_cs_data_file(config["Input File Path"],
                                     config["Max Theta"])
    print "Raw data is read in"
    print exp_data[0][0], "\n", exp_data[0][1]
    ivgdr_frac = None
    if config["Subtract IVGDR"]:
        ivgdr_frac = IVGDRFraction(config["IVGDR Height"],
                                   config["IVGDR Center"],
                                   config["IVGDR Width"],
                                   config["IVGDR CS Integral"])
        print "IVGDR EWSR % calculation is prepared"
    else:
        ivgdr_frac = IVGDRFraction(0.0, 12.0, 4.0, 20.0)
        print "IVGDR EWSR % calculation is turned off"
    print "Ex =", exp_data[0][0], "MeV IVGDR EWSR %=",\
          ivgdr_frac(exp_data[0][0])


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


def read_row_cs_data_file(file_name, max_angle):
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
        # convert each cell to a float
        distibution_data = [x.strip() for x in data_list[1:]]
        # distibution_data = map(lambda x: float(x.strip()), data_list[1:])
        dist_len = len(distibution_data)
        distribution = []
        i = 0
        # for each trio of elements in the distribution put them in a tuple and
        # put that tuple in the distribution list this makes the distribution
        # variable contain a list of points and errors:
        # [(x1,y1,dy1),(x2,y2,dy2),...]
        for i in xrange(0, dist_len, 3):
            if distibution_data[i] <= max_angle:
                distribution.append((distibution_data[i],
                                     distibution_data[i+1],
                                     distibution_data[i+2]))
            # if 15.4 < energy and energy < 15.6:
            #    print distData[i],",",distData[i+1],",",distData[i+2]
        # put the energy and its associated distribution in a list
        output.append([energy, distribution])
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
