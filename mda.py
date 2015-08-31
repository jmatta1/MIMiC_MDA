#!/usr/bin/python
"""This program performs a multipole decomposition analysis of the data and
options given in the configuration file mda_config.py and using the Markov
Chain Monte Carlo method to sample and get error bars and best fits"""
import sys  # for exit
import multiprocessing  # for cpu_count
from math import ceil
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
    samples_needed = int(ceil(10.0 /
                              ((1.0 - config["Confidence Interval"]) / 2.0)))
    if samples_needed > num_samples:
        print SAMPLES_ERROR % (num_samples, config["Sample Points"],
                               config["Number of Walkers"],  samples_needed, 
                               config["Confidence Interval"])
        sys.exit()


SAMPLES_ERROR = """
WARNING: the number of samples:
%d = ((%d - 50) * %d)
is not large enough to ensure 10 points outside of each error bar location
for the given confidence interval.

the minimum number of samples necessary is:

%d = ceiling[ 10 / ((1.0 - %12.10f) / 2.0))
"""



if __name__ == "__main__":
    main()
