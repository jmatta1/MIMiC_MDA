#!/usr/bin/python
"""This program performs a multipole decomposition analysis of the data and
options given in the configuration file. First it finds starting points using
the BFGS algorithm from a variety of starting positions given in the config
file. Then it takes those starting points, refines them and finds parameter
errors using the Markov Chain Monte Carlo technique, specifically the python
implementation of Goodman & Weare's Affine Invariant MCMC ensemble sampler
which is in the emcee package"""
import sys
import multiprocessing
import math
import numpy as np
import copy
import ctypes as ct
from scipy import interpolate
from scipy import optimize

# firsts we check the command line and grab the module name
if len(sys.argv) != 2:
    print "\nUsage:\n\t./mda.py configuration_file\n\t  or"
    print "\tpython mda.py configuration_file\n"
    sys.exit()
# strip off the .py if it exists
CF_FILE_NAME = None
if sys.argv[1][-3:] == ".py":
    CF_FILE_NAME = sys.argv[1][0:-3]
else:
    CF_FILE_NAME = sys.argv[1]
# prevent bytecode generation for the config file
ORIGINAL_SYS_DONT_WRITE_BYTECODE = sys.dont_write_bytecode
sys.dont_write_bytecode = True
# import the config file
CONFIG = __import__(CF_FILE_NAME).CONFIG
# restore the dont_write_bytecode variable to its original value
sys.dont_write_bytecode = ORIGINAL_SYS_DONT_WRITE_BYTECODE


def main():
    """gets the configuration file to import, imports it and then performs
    sanity checks on the configuration data and then calls the functions that
    do the work"""
    # check that the user gave sane information
    # check that they are not requesting greater concurrency than the
    # system supports
    cpu_count = multiprocessing.cpu_count()
    if CONFIG["Number of Threads"] > cpu_count:
        print "\nInvalid number of threads, on this machine it must be: "
        print 'CONFIG["Number of Threads"] <= %d\n' % cpu_count
        sys.exit()
    num_samples = ((CONFIG["Sample Points"] - 50)
                   * CONFIG["Number of Walkers"])
    samples_needed = int(math.ceil(10.0 /
                                   ((1.0 - CONFIG["Confidence Interval"]) /
                                    2.0)))
    if samples_needed > num_samples:
        print SAMPLES_ERROR % (num_samples, CONFIG["Sample Points"],
                               CONFIG["Number of Walkers"], samples_needed,
                               CONFIG["Confidence Interval"])
        sys.exit()
    num_dists = (1 + CONFIG["Maximum L"])
    num_ewsr = len(CONFIG["EWSR Fractions"])
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
    exp_data = read_row_cs_data_file(CONFIG["Input File Path"],
                                     CONFIG["Max Theta"],
                                     CONFIG["Start Energy"],
                                     CONFIG["Final Energy"])
    # now read and subtract the IVGDR data
    ivgdr_dists, ivgdr_ewsr, sub_data = handle_ivgdr(exp_data)
    # print ivgdr_dists[0], '\n', ivgdr_ewsr[0], '\n', sub_data[0]
    # now read the distributions that are used to fit the data
    dists = [[read_dist(elem[0], i) for i in range(CONFIG["Maximum L"] + 1)]
             for elem in exp_data]
    print "Distributions read in"
    # now interpolate the distributions to get the values at the angles in data
    interp_dists = interp_all_dists(dists, exp_data)
    print "Distributions interpolated and prepared for fitting"
    # now get the data divided by errors without angle values
    fit_data = [(exp_en[1][:, 1]/exp_en[1][:, 2]) for exp_en in exp_data]
    print "Experimental data prepared for fitting"
    # calculate the starting parameter sets for initial searches
    start_params = calc_start_params()
    print "Finished calculating parameter starting point list"""
    # now interleave things so we are ready to use pool.map across everything
    interleaved_data = [(fit_data[i], interp_dists[i], start_params)
                        for i in range(len(fit_data))]
    print "Data is interleaved"
    mp_pool = multiprocessing.Pool(processes=CONFIG["Number of Threads"])
    print ("Starting MDA process, working on %d energies simultaneously" %
           CONFIG["Number of Threads"])
    # output = mp_pool.map(fit_and_mcmc, interleaved_data)
    # single threaded version for debugging
    output = map(fit_and_mcmc, interleaved_data)


def fit_and_mcmc(data_tuple):
    """This function is the workhorse function, it loads the shared library,
    generates the structure, performs the BFGS fit from each starting point
    and then performs the MCMC using those fits"""
    # first unpack the tuple
    fit_data = data_tuple[0]
    interp_dists = data_tuple[1]
    start_points = data_tuple[2]
    # now load the shared library
    cs_lib = ct.cdll.LoadLibrary(CONFIG["Shared Lib Path"])
    # set the return types
    cs_lib.makeMdaStruct.restype = ct.c_void_p
    cs_lib.calculateChi.restype = ct.c_float
    cs_lib.calculateLnLiklihood.restype = ct.c_float
    # build the calculation object
    struct = make_calc_struct(cs_lib, fit_data, interp_dists)
    # now map the perform fit function onto the set of starting points
    final_points = map(lambda x: do_init_fit(x, struct, cs_lib), start_points)
    # clusterize the final points
    centers = clusterize_points(final_points)
    # generate the starting points for the walkers
    starts = gen_walker_starts(centers)
    # perform the MCMC
    # write the triangle plots
    # extract the error bars
    # delete the struct
    cs_lib.freeMdaStruct(struct)
    #return the point and the errors
    return None


def gen_walker_starts(centers):
    """This function takes the set of cluster centers and generates the walker
    start positions from them"""
    return centers


def clusterize_points(points):
    """Takes the points given and tries to find clusters"""
    return points


def do_init_fit(start, struct, cs_lib):
    """Performs a fit from the given starting point"""
    # calculate the bounds
    bnds = [(0.0, 1.0/x) for x in CONFIG["EWSR Fractions"]]
    print call_chi_sq(np.array(start, dtype=np.float32), cs_lib, struct)
    init_params = np.array(start, dtype = np.float32)
    min_pos, min_chi, d = optimize.fmin_l_bfgs_b(call_chi_sq,
                                                 init_params,
                                                 bounds = bnds,
                                                 epsilon = 1e-03,
                                                 approx_grad = True,
                                                 args = (cs_lib, struct),
                                                 iprint = 0, factr = 10.0)
    print "bfgs done"
    print start, min_pos, min_chi
    return start


def call_chi_sq(params, cs_lib, struct):
    """calls the chi^2 function in cs_lib"""
    print params
    temp = cs_lib.calculateChi(struct, params.ctypes.data)
    print temp
    return temp


def make_calc_struct(cs_lib, data, dists):
    """This function takes the data and distributions and dumps the information
    into a freshly created struct"""
    # make a struct
    out_struct = cs_lib.makeMdaStruct(len(data), len(dists))
    # load it with the data
    cs_lib.setMdaData(out_struct, data.ctypes.data)
    # iterate through this distributions
    for i in range(len(dists)):
        # load the distributions
        cs_lib.setMdaDist(out_struct, i, dists[i].ctypes.data)
    # return the struct
    return out_struct


def calc_start_params():
    """This function uses the config data to generate a list of starting
    points for the initial fits performed before the MCMC"""
    # first load the sets into an array
    sl_list = []
    for i in range(CONFIG["Maximum L"]+1):
        sl_list.append(CONFIG[("Start Pts a%d" % i)])
    # now compute the total number of starting points
    # also get a list of the lengths of each starting list
    # also make a starting list of indices with everything set to 0
    num_starts = 1
    lengths = []
    indices = []
    for start_list in sl_list:
        lengths.append(len(start_list))
        num_starts *= len(start_list)
        indices.append(0)
    # now compute the list of starting points
    start_list = []
    for _ in range(num_starts):
        # first construct the start list based on the current indices
        curr_start = []
        for i in range(len(indices)):
            curr_start.append(sl_list[i][indices[i]])
        start_list.append(curr_start)
        # now increment the indices
        increment_ind(indices, lengths)
    # now return the list of starting points
    return start_list


def increment_ind(ind, lens):
    """This function takes a list of indices and lengths and increments the
    highest indice, resetting to zero as needed"""
    start = len(ind)-1
    for i in range(start, -1, -1):
        if (ind[i] + 1) < lens[i]:
            ind[i] += 1
            break
        else:
            ind[i] = 0


def interp_all_dists(dists, data):
    """this function takes the grand list of distributions and the experimental
    data and calls another function to interpolate the distributions with the
    exp angles and divide the interpolated distributions by the appropriate
    error for each point"""
    # first make the output variable
    output = []
    # now iterate across the energies
    for i in range(len(dists)):
        # extract the list of angles for this energy
        angle_list = data[i][1][:, 0]
        # extract the list of errorss for this energy
        error_list = data[i][1][:, 2]
        # make the variable to hold the list of interpolated values
        en_output = []
        # iterate across the L values
        for j in range(len(dists[i])):
            # interpolate and divide the distribution
            interp_data = interpolate_dist(dists[i][j], angle_list,
                                           error_list).astype(np.float32)
            en_output.append(interp_data)
        # append the information for the distributions of this energy to output
        output.append(en_output)
    return output


def interpolate_dist(dist, angles, errors):
    """this function takes a distribution, a list of exp angles, and a list of
    exp error at each angle, it then interpolates the distribution at that
    angle and divides that value by the exp error at that angle"""
    # construct the interpolation
    interp = interpolate.interp1d(dist[:, 0], dist[:, 1], kind="cubic")
    # get the interpolated values divided by the errors and return them
    return interp(angles)/errors


def read_dist(energy, l_value):
    """This function reads the distribution described by the passed parameters
    from disk into an np array"""
    # first construct the file name
    dist_file_name = "{0:s}A{1:d}_Ex{2:4.2f}_L{3:02d}_T0_F{4:03d}.csv".format(
        CONFIG["Distribution Directory"], CONFIG["Target A"], energy, l_value,
        int(100.0*CONFIG["EWSR Fractions"][l_value]))
    dist_file = open(dist_file_name, 'r')
    output = []
    # iterate through the file
    for line in dist_file:
        vals = [float(x.strip()) for x in line.strip().split(',')]
        output.append((vals[0], vals[1]))
    return np.array(output, dtype=np.float32)


def handle_ivgdr(data):
    """This function handles the ivgdr subtraction if no subtraction is needed
    then this function returns None, None, and a deep copy of the exp data
    that was passed to it, otherwise it returns a list of distributions,
    ewsr fractions and the experimental data minus the interpolated IVGDR data
    times the ewsr fraction"""
    # if needed handle the ivgdr
    if CONFIG["Subtract IVGDR"]:
        frac = IVGDRFraction(CONFIG["IVGDR Height"], CONFIG["IVGDR Center"],
                             CONFIG["IVGDR Width"],
                             CONFIG["IVGDR CS Integral"])
        # calculate the list of IVGDR fractions
        ewsrs = [frac.get_ivgdr_fraction(elem[0]) for elem in data]
        # read the distributions
        dists = read_ivgdr_dists([elem[0] for elem in data])
        # extract the list of angles for each energy
        angle_lists = [en_set[1][:, 0] for en_set in data]
        # calculate scaled and interpolated points
        ivgdr_values = [interpolate_and_scale_ivgdr(dist, ewsr, angle_list)
                        for (dist, ewsr, angle_list) in zip(dists, ewsrs,
                                                            angle_lists)]
        sub_data = copy.deepcopy(data)
        for i in range(len(sub_data)):
            sub_data[i][1][:, 1] = sub_data[i][1][:, 1] - ivgdr_values[i]
        print "IVGDR subtraction is done"
        return dists, ewsrs, sub_data
    else:
        print "IVGDR is turned off"
        return None, None, copy.deepcopy(data)


def interpolate_and_scale_ivgdr(dist, ewsr, angle_list):
    """This function takes an ivgdr distribution, interpolates it, calculates
    the values of the distribution at the provided angles, multiplies them by
    ewsr and returns it"""
    interp = interpolate.interp1d(dist[:, 0], dist[:, 1], kind="cubic")
    values = interp(angle_list)
    return ewsr*values


def read_ivgdr_dists(en_list):
    """reads in the csv files with the IVGDR distributions"""
    dist_file_names = \
        ["{0:s}A{1:d}_Ex{2:4.2f}_L01_T1_F100.csv".format(
            CONFIG["Distribution Directory"], CONFIG["Target A"], energy)
         for energy in en_list]
    dists_list = []
    for name in dist_file_names:
        dist_file = open(name, 'r')
        dist = [[float(elem) for elem in line.strip().split(',')]
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
    print "Exp Data is read in"
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
