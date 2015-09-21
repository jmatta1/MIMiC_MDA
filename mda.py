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
import copy
import emcee
import os
import corner
import matplotlib.pyplot as plt
import numpy as np
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

# TODO: peak find based parameters
# TODO: implement parameter plot function
# TODO: implement parameter writer

PLOT_FORMAT_LIST = ["svg", "svgz", "pdf", "ps", "eps", "png"]


def main():
    """performs sanity checks on the configuration data and then calls the
    functions that do the work"""
    # check that the user gave sane information
    # check that they are not requesting greater concurrency than the
    # system supports
    cpu_count = multiprocessing.cpu_count()
    if CONFIG["Number of Threads"] > cpu_count:
        print "\nInvalid number of threads, on this machine it must be: "
        print 'CONFIG["Number of Threads"] <= %d\n' % cpu_count
        sys.exit()
    # check if the user set the sampling high enough for the error bars wished
    num_samples = ((CONFIG["Sample Points"] - CONFIG["Burn-in Points"])
                   * CONFIG["Number of Walkers"])
    samples_needed = int(math.ceil(10.0 /
                                   ((1.0 - CONFIG["Confidence Interval"]) /
                                    2.0)))
    if samples_needed > num_samples:
        print SAMPLES_ERROR % (num_samples, CONFIG["Sample Points"],
                               CONFIG["Burn-in Points"],
                               CONFIG["Number of Walkers"], samples_needed,
                               CONFIG["Confidence Interval"])
        sys.exit()
    # make certain the user gave enough EWSR fractions for the max L
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
            "(1 + CONFIGURATION[\"Maximum L\"]) =pn %d\n" % num_dists,\
            "EWSR fractions listed in the variable",\
            "CONFIGURATION[\"EWSR Fractions\"], %d were given\n" % num_ewsr
        sys.exit()
    # check to make certain that there are at least 10 start points
    len_array = [len(CONFIG["Start Pts a%d" % i]) for i in range(num_dists)]
    num_cells = 1
    for size in len_array:
        num_cells *= size
    if num_cells < CONFIG["Number Walker Generators"]:
        out_str = """\nYou must provide enough starting points such that there are at
least %d points (the number of start points is the length of each start list
multiplied together)""" % CONFIG["Number Walker Generators"]
        print out_str
        sys.exit()
    # check to make certain that num sampes is a multiple of min start points
    if (num_samples % CONFIG["Number Walker Generators"]) != 0:
        print '\nThe product ((CONFIG["Sample Points"]-CONFIG["Burn-in Points"])'\
            '*CONFIG["Number of Walkers"])\n must be a multiple of %d' %\
            CONFIG["Number Walker Generators"]
        sys.exit()
    # check to make certain that the given file format is one of the
    # supported formats
    if not CONFIG["Plot Format"] in PLOT_FORMAT_LIST:
        print "\nThe chosen plot output format is not supported."
        print "The supported values for this option are:"
        print PLOT_FORMAT_LIST, "\n"
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
    ivgdr_info = handle_ivgdr(exp_data)
    sub_data = ivgdr_info[2]
    # print ivgdr_dists[0], '\n', ivgdr_ewsr[0], '\n', sub_data[0]
    # now read the distributions that are used to fit the data
    dists = [[read_dist(elem[0], i) for i in range(CONFIG["Maximum L"] + 1)]
             for elem in exp_data]
    print "Distributions read in"
    # now interpolate the distributions to get the values at the angles in data
    interp_dists = interp_all_dists(dists, exp_data)
    print "Distributions interpolated and prepared for fitting"
    # now get the data divided by errors without angle values
    fit_data = [(exp_en[1][:, 1]/exp_en[1][:, 2]) for exp_en in sub_data]
    print "Experimental data prepared for fitting"
    # calculate the starting parameter sets for initial searches
    start_params = calc_start_params()
    print "Finished calculating parameter starting point list"""
    # now interleave things so we are ready to use pool.map across everything
    interleaved_data = [(exp_data[i][0], fit_data[i], interp_dists[i],
                         start_params) for i in range(len(fit_data))]
    print "Data is interleaved"
    generate_output_dirs()
    # mp_pool = multiprocessing.Pool(processes=CONFIG["Number of Threads"])
    print ("Starting MDA process, working on up to %d energies simultaneously"
           % CONFIG["Number of Threads"])
    # parameters = mp_pool.map(fit_and_mcmc, interleaved_data)
    # single threaded version for debugging
    parameters = map(fit_and_mcmc, interleaved_data)
    # write the individual fits to csv files
    print "Writing fit files"
    write_fits(exp_data, dists, parameters, ivgdr_info)
    # make the fit plots
    print "Writing fit plots"
    make_fit_plots(exp_data, dists, parameters, ivgdr_info)
    # write the two parameter sets
    energy_set = [val[0] for val in exp_data]
    print "Writing parameter sets"
    write_param_sets(parameters, energy_set)
    # write the parameter plots
    print "Writing parameter plots"
    write_param_plots(parameters, energy_set)
    # and now we are completely done
    print "MDA and output complete"


def write_param_plots(parameters, energy_set):
    """This function takes the generated parameters and makes the plots for
    each L of the parameters"""
    # first split the data into the two types
    perc_data = [pset[0] for pset in parameters]
    peak_data = [pset[1] for pset in parameters]


def write_param_sets(parameters, energy_set):
    """This function takes the generated parameters and writes the sets from
    percentiles and the sets from peak find to seperate files"""
    # first split the data into the two types
    perc_data = [pset[0] for pset in parameters]
    peak_data = [pset[1] for pset in parameters]
    # generate the file names
    perc_path = CONFIG["Parameter Files Directory"]
    peak_path = CONFIG["Parameter Files Directory"]
    if CONFIG["Parameter Files Directory"][-1] != "/":
        perc_path += "/"
        peak_path += "/"
    perc_path += ("A%d_percentile_parameters.csv" % CONFIG["Target A"])
    peak_path += ("A%d_peak_parameters.csv" % CONFIG["Target A"])
    # now call the function that writes a parameter set
    write_horizontal_param_file(perc_path, perc_data, energy_set)
    write_horizontal_param_file(peak_path, peak_data, energy_set)


def write_horizontal_param_file(path, params, energy_set):
    """This function takes a set of parameters and their corresponding energies
    and writes the data to the specified path"""
    # open the file for writing
    out_file = open(path, 'w')
    # write the header
    out_file.write(generate_param_file_header())
    # write each line of parameters and energies
    for i in range(len(params)):
        out_file.write(gen_param_file_line(energy_set[i], params[i]))
    # close the file we wrote
    out_file.close()


def gen_param_file_line(energy, params):
    """This function takes the excitation energy and the parameters for that
    energy and generates the string to be written to the file"""
    out_str = ("%f, , " % energy)
    for i in range((CONFIG["Maximum L"]+1)):
        out_str += ("%f, %f, %f, , " % params[i])
    out_str += "\n"
    return out_str


def generate_param_file_header():
    """This function uses config information to generate an appropriate header
    for the parameter file"""
    out_str = "MDA parameters extracted with percentile techniques\n"
    hrow1 = "Excitation, , "
    hrow2 = "Energy, , "
    for i in range((CONFIG["Maximum L"]+1)):
        hrow1 += ("a%d, a%d, a%d, , " % (i, i, i))
        hrow2 += "Value, Pos. Err, Neg. Err., , "
    out_str += hrow1
    out_str += "\n"
    out_str += hrow2
    out_str += "\n"
    return out_str


def make_fit_plots(data, dists, parameters, ivgdr_info):
    """This function takes everything and generates the plots for individual
    fits at each energy"""
    # first make the directory names
    legend = [r"$Fit$"]
    for i in range(len(parameters[0][0])):
        legend.append(r"$l_{%d}$" % i)
    # loop through each set of data and distributions
    for i in range(len(data)):
        # extract the things pertinet to this set from the arguments
        energy = data[i][0]
        dist_set = copy.deepcopy(dists[i])
        param_set = list(copy.deepcopy(parameters[i]))
        exp_points = data[i][1]
        # handle ivgdr
        if CONFIG["Subtract IVGDR"]:
            dist_set.append(ivgdr_info[0][i])
            param_set.append((ivgdr_info[1][i], 0.0, 0.0))
            legend.append(r"$l_{-1}$")
        # loop through the plots to generate
        for i in range(2):  # peaks or percentiles?
            # choose our parameter set and generate the three sets of fits
            pset = gen_param_sets_for_fit_plot(param_set[i])
            for k in range(3):  # params-lo_errs, params, params+hi_errs
                sc_dists = gen_fit_dists(pset[k], dist_set)
                for j in range(2):  # full or limitted
                    plot_path = "%sA%d_E%4.1f.%s" % (CONFIG["Fit Plot Dirs"]
                                                     [6*i+3*j+k],
                                                     CONFIG["Target A"],
                                                     energy,
                                                     CONFIG["Plot Format"])
                    # now decide how much of the distributions to call the
                    # gen fit plot function on
                    if j == 1:
                        gen_fit_plot(exp_points,
                                     sc_dists[:(CONFIG["Fit Plot L Limit"]+2)],
                                     legend[:(CONFIG["Fit Plot L Limit"]+2)],
                                     plot_path)
                    elif (CONFIG["Subtract IVGDR"] and
                          not CONFIG["Plot IVGDR in Fits"]):
                        gen_fit_plot(exp_points,
                                     sc_dists[:(len(sc_dists)-1)],
                                     legend[:(len(sc_dists)-1)],
                                     plot_path)
                    else:
                        gen_fit_plot(exp_points,
                                     sc_dists[:len(sc_dists)],
                                     legend[:len(sc_dists)],
                                     plot_path)


def gen_param_sets_for_fit_plot(params):
    """This function takes a single set of parameters and generates three sets
    of parameters, the first, decreased by the lower error bar, the second
    equal to the parameter fit value, and the third, increased by the upper
    error bar"""
    temp = [[(vals[0] - vals[2]) for vals in params],
            [vals[0] for vals in params],
            [(vals[0] + vals[1]) for vals in params]]
    for pset in temp:
        for i in range(len(pset)):
            if pset[i] < 0.0:
                print "Got a negative parameter when subtracting errors from"
                print "parameters for fit plots, this should not be possible"
                sys.exit()
                pset[i] = 0.0
    return temp


def gen_fit_plot(points, dists, legends, plot_name):
    """This function takes the list of points with errors in the points var
    and the list of distributions to plot in the dists var and generates a
    nicely formatted matplotlib plot displaying them"""
    # since the first distribution is always the plot then the first style will
    # always be a solid red line, no other distribution is red and solid
    line_styles = ["r-", "b--", "g-.", "c:", "m--", "y-.", "b:", "g--", "c-.",
                   "m:", "y--", "b-.", "g:", "c--", "m-.", "y:"]
    pt_x_vals = points[:, 0]
    pt_y_vals = points[:, 1]
    pt_e_vals = points[:, 2]
    fig, axes = plt.subplots()
    axes.set_yscale('log', subsy=[2, 3, 4, 5, 6, 7, 8, 9])
    axes.set_xscale('linear')
    # plot the points
    axes.errorbar(pt_x_vals, pt_y_vals, yerr=pt_e_vals, fmt="ko",
                  label=r"$Exp$", markersize=2.0)
    # plot the distributions
    for i in range(len(dists)):
        axes.plot(dists[i][:, 0], dists[i][:, 1],
                  line_styles[i % len(line_styles)], label=legends[i])
    # set the scale of the x axis
    axes.set_xlim(0.0, math.ceil(pt_x_vals.max()))
    # set the scale of the y axis
    (ymin_val, ymax_val) = find_y_extrema(pt_y_vals.max(), dists,
                                          math.ceil(pt_x_vals.max()))
    axes.set_ylim(ymin_val, ymax_val)
    # label the axes
    axes.set_xlabel(r'Lab Angle $(^{\circ{}})$')
    axes.set_ylabel(r'$(\partial^2 \sigma)/(\partial \Omega \partial E)$ ($mb/(sr*MeV)$)')
    # make the legend
    legend = axes.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1)
    legend.get_frame().set_facecolor("white")
    # save and close the figure
    fig.savefig(plot_name, additional_artists=[legend], bbox_inches='tight')
    plt.close(fig)


def find_y_extrema(data_max, dists, xmax):
    """This function scans the distributions provided searching for the minimum
    value with angle less than xmax, it also looks to find the maximum value
    (be it in a distribution or the data maximum it then returns
    (10^(floor(log(ymin))), 10^(ceiling(log(ymax)))"""
    current_min = 100000000000.0
    current_max = data_max
    for dist in dists:
        for point in dist:
            if point[0] > xmax:
                break
            if point[1] < current_min and point[1] > 0.0:
                current_min = point[1]
            elif point[1] > current_max:
                current_max = point[1]
    log_min = math.floor(2.0*math.log10(current_min))/2.0
    log_max = math.ceil(2.0*math.log10(current_max))/2.0
    return (math.pow(10.0, log_min), math.pow(10.0, log_max))


def write_fits(data, dists, parameters, ivgdr_info):
    """This function takes the parameters, distributions, and data, and writes
    them to a nicely formatted csv file for usage later"""
    # first split the data up into individual runs
    for i in range(len(data)):
        energy = copy.deepcopy(data[i][0])
        points = copy.deepcopy(data[i][1])
        dist_set = copy.deepcopy(dists[i])
        perc_set = copy.deepcopy(parameters[i][0])
        peak_set = copy.deepcopy(parameters[i][1])
        # test if we are subtracting the IVGDR
        if CONFIG["Subtract IVGDR"]:
            dist_set.append(ivgdr_info[0][i])
            perc_set.append((ivgdr_info[1][i], 0.0, 0.0))
            peak_set.append((ivgdr_info[1][i], 0.0, 0.0))
        # calculate the file names
        file_paths = ["%sA%d_E%4.1f.csv" % (CONFIG["Fit Csv Dirs"][0],
                                            CONFIG["Target A"], energy),
                      "%sA%d_E%4.1f.csv" % (CONFIG["Fit Csv Dirs"][1],
                                            CONFIG["Target A"], energy)]
        # write the percentile fit
        write_fit_csv(file_paths[0], points, perc_set, dist_set, energy)
        # write the peak find fit
        write_fit_csv(file_paths[1], points, peak_set, dist_set, energy)


def write_fit_csv(path, points, pset, dist_set, energy):
    """This function takes a file path, a set of experimental data, a parameter
    set, a set of distributions and the ex energy of the distribution and
    writes a nicely formated csv file with all the information in it"""
    # first open the file to be written
    out_file = open(path, 'w')
    # next generate the title and column headings for the file
    out_file.write(gen_csv_title_and_headings(energy))
    # construct the list of lines in the file
    num_lines = len(dist_set[0])
    if num_lines < len(points):
        num_lines = len(points)
    if num_lines < len(pset):
        num_lines = len(pset)
    csv_list = []
    for _ in range(num_lines):
        csv_list.append("")
    # append the exp data
    append_exp_data_to_fit_csv(points, csv_list)
    # append a blank column
    append_str_to_fit_csv(", ", csv_list)
    # append the parameter data
    append_parameters_to_fit_csv(pset, csv_list)
    # append a blank column
    append_str_to_fit_csv(", ", csv_list)
    # calculate the scaled distributions
    params = [param[0] for param in pset]
    scaled_dists = gen_fit_dists(params, dist_set)
    # append the distribution angles
    append_data_column_to_fit_csv(dist_set[0][:, 0], csv_list)
    # append each distribution, and scaled distribution
    for i in range(len(dist_set)):
        append_data_column_to_fit_csv(dist_set[i][:, 1], csv_list)
        append_data_column_to_fit_csv(scaled_dists[i+1][:, 1], csv_list)
    # append a blank column
    append_str_to_fit_csv(", ", csv_list)
    # append the fit distribution
    append_data_column_to_fit_csv(scaled_dists[0][:, 1], csv_list)
    # append the newline characters
    append_str_to_fit_csv("\n", csv_list)
    # write the csv list
    for line in csv_list:
        out_file.write(line)
    # close the output file
    out_file.close()


def append_data_column_to_fit_csv(data, csv_list):
    """This function takes a column of data and appends it to the csv list"""
    for i in range(len(csv_list)):
        if i < len(data):
            csv_list[i] += ("%f, " % data[i])
        else:
            csv_list[i] += ", "


def append_parameters_to_fit_csv(pset, csv_list):
    """This function appends a parameter set to the csv list"""
    # first generate the list of names
    name_list = []
    for i in range(len(pset)):
        if CONFIG["Subtract IVGDR"] and i == (len(pset)-1):
            name_list.append("a-1")
        else:
            name_list.append("a%d" % i)
    # now append the names and error bars
    for i in range(len(csv_list)):
        if i < len(pset):
            csv_list[i] += ("%s, %f, %f, %f, " % (name_list[i], pset[i][0],
                                                  pset[i][1], pset[i][2]))
        else:
            csv_list[i] += ", , , , "


def append_str_to_fit_csv(str_to_append, csv_list):
    """This function appends the given string to the csv_list"""
    for i in range(len(csv_list)):
        csv_list[i] += str_to_append


def append_exp_data_to_fit_csv(points, csv_list):
    """This function takes an experimental data set and the csv list and writes
    the data to the csv list along with a set of appropriately blank lines"""
    for i in range(len(csv_list)):
        if i < len(points):
            csv_list[i] += ("%f, %f, %f, " % (points[i][0], points[i][1],
                                              points[i][2]))
        else:
            csv_list[i] += " , , , "


def gen_csv_title_and_headings(energy):
    """This function takes the excitation energy and returns a string with the
    csv title, and column headings"""
    # first generate the title string
    out_str = "Fit Information for Ex =, %4.1f\n" % energy
    # now generate the two rows of column headings
    # first the easy column headings
    row1 = "Exp., Exp., Exp., , Param, Param., Pos., Neg., , Dist., "
    row2 = "Angle, CS, Error, , Name, Value, Error, Error, , Angle, "
    # write the headings for the distributions
    for i in range((CONFIG["Maximum L"]+1)):
        row1 += ", Scaled, "
        row2 += "l%d, l%d, " % (i, i)
    # if needed write the heading for the ivgrd
    if CONFIG["Subtract IVGDR"]:
        row1 += ", Scaled, "
        row2 += "IVGDR, IVGDR, "
    # write the column headings for the fit distribution
    row1 += ", Fit (Sum of, \n"
    row2 += ", Scaled Dists), \n"
    # add the two rows to out_str and return it
    out_str += row1
    out_str += row2
    return out_str


def gen_fit_dists(params, dists):
    """This function, takes a parameter set and a set of distributions, it then
    scales the distributions by the appropriate parameters and returns the
    scaled distributions"""
    # first scale all the passed distributions
    sc_dists = copy.deepcopy(dists)
    for (param, dist) in zip(params, sc_dists):
        dist[:, 1] = param*dist[:, 1]
    # make the fit distribution and initiailize it
    fit_distribution = []
    for i in range(len(sc_dists[0])):
        fit_distribution.append([sc_dists[0][i][0], sc_dists[0][i][1]])
    # add the other components of the fit distribution
    for j in range(1, len(sc_dists)):
        for i in range(len(sc_dists[j])):
            fit_distribution[i][1] += sc_dists[j][i][1]
    output = [np.array(fit_distribution)]
    for dist in sc_dists:
        output.append(dist)
    return output


def fit_and_mcmc(data_tuple):
    """This function is the workhorse function, it loads the shared library,
    generates the structure, performs the BFGS fit from each starting point
    and then performs the MCMC using those fits"""
    # first unpack the tuple
    energy = data_tuple[0]
    fit_data = data_tuple[1]
    interp_dists = data_tuple[2]
    start_points = data_tuple[3]
    print "Starting work on Energy", energy, "MeV"
    # now load the shared library
    cs_lib = ct.cdll.LoadLibrary(CONFIG["Shared Lib Path"])
    # set the return types
    cs_lib.makeMdaStruct.restype = ct.c_void_p
    cs_lib.calculateChi.restype = ct.c_double
    cs_lib.calculateLnLiklihood.restype = ct.c_double
    # build the calculation object
    struct = make_calc_struct(cs_lib, fit_data, interp_dists)
    print "Commencing initial fits for", energy, "MeV"
    # now map the perform fit function onto the set of starting points
    final_points = [do_init_fit(x, struct, cs_lib) for x in start_points]
    print "Finished initial fits for", energy, "MeV"
    # sort the list of starting points in order of chi^2
    final_points.sort(key=lambda x: x[0])
    # then take the first CONFIG["Number Walker Generators"]
    generators = [x[1] for x in
                  final_points[0:CONFIG["Number Walker Generators"]]]
    # generate the starting points for the walkers
    starts = gen_walker_starts(generators)
    # calculate the boundaries of the ensemble
    bnds = [(0.0, 1.0/x) for x in CONFIG["EWSR Fractions"]]
    # create the sampling ensemble
    ndims = (1 + CONFIG["Maximum L"])
    print "Commencing MCMC for", energy, "MeV"
    sampler = emcee.EnsembleSampler(CONFIG["Number of Walkers"], ndims,
                                    ln_post_prob, args=(cs_lib, struct, bnds))
    # perform the MCMC
    sampler.run_mcmc(starts, CONFIG["Sample Points"])
    print "Finished MCMC for", energy, "MeV"
    # delete the struct
    cs_lib.freeMdaStruct(struct)
    # now do the rest, all of which involves calculating things from samples
    return perform_sample_manips(sampler, ndims, energy)


def perform_sample_manips(sampler, ndims, energy):
    """this function takes the mcmc sampler, grabs the chains from it, and
    generates plots and percentiles and most likely values from the sampled
    points"""
    # retrieve the samples
    num_samples = (CONFIG["Number of Walkers"] * (CONFIG["Sample Points"] -
                                                  CONFIG["Burn-in Points"]))
    samples = sampler.chain[:, CONFIG["Burn-in Points"]:, :].reshape((
        num_samples, ndims))
    # save the samples to the disk in the sp
    if CONFIG["Save Chain Data"]:
        print "Saving MCMC samples for", energy, "MeV"
        chain_file_name = ""
        if CONFIG["Chain Directory"][-1] == '/':
            chain_file_name = CONFIG["Chain Directory"] +\
                "A%d_chain_E%4.1f.%s" % (CONFIG["Target A"], energy,
                                         CONFIG["Plot Format"])
        else:
            chain_file_name = CONFIG["Chain Directory"] +\
                "/A%d_chain_E%4.1f.%s" % (CONFIG["Target A"], energy,
                                          CONFIG["Plot Format"])
        np.savez_compressed(chain_file_name, sampler.chain)
        print "Done saving MCMC samples for", energy, "MeV"
    # make the probability plots
    make_prob_plots(samples, energy)
    # extract the error bars
    quantile_list = np.array([(0.5 - CONFIG["Confidence Interval"] / 2.0), 0.5,
                              (0.5 + CONFIG["Confidence Interval"] / 2.0)])
    # make the corner plot
    if CONFIG["Generate Corner Plots"]:
        print "Commencing corner plot creation for", energy, "MeV"
        lbls = [r"$a_{%d}$" % i for i in range(ndims)]
        ranges = [(-0.0001, (0.0001 + samples[:, i].max())) for i in
                  range(ndims)]
        fig = None
        if CONFIG["Corner Plot Samples"] >= num_samples:
            fig = corner.corner(samples, labels=lbls, range=ranges,
                                quantiles=quantile_list, verbose=False)
        else:
            # randomize the sample array and then extract the first chunk of it
            np.random.shuffle(samples)
            temp_samples = samples[0:CONFIG["Corner Plot Samples"]]
            fig = corner.corner(temp_samples, labels=lbls, range=ranges,
                                quantiles=quantile_list, verbose=False)
        # make the corner plot file_name
        if CONFIG["Corner Plots Directory"][-1] == '/':
            fig_file_name = CONFIG["Corner Plots Directory"] +\
                "A%d_corner_E%4.1f.%s" % (CONFIG["Target A"], energy,
                                          CONFIG["Plot Format"])
        else:
            fig_file_name = CONFIG["Corner Plots Directory"] +\
                "/A%d_corner_E%4.1f.%s" % (CONFIG["Target A"], energy,
                                           CONFIG["Plot Format"])
        fig.savefig(fig_file_name, bbox_inches='tight')
        plt.close(fig)
        print "Done creating corner plot for", energy, "MeV"
    # return the point and the errors
    return calc_param_values(samples, quantile_list, ndims)


def calc_param_values(samples, quantile_list, ndims):
    """This function calculates the parameter values and error bars using both
    peak finding and percentiles"""
    points = [(v[1], v[2]-v[1], v[1]-v[0]) for v in
              zip(*np.percentile(samples, (100.0*quantile_list), axis=0))]
    peaks = find_most_likely_values(samples, ndims)
    return (points, peaks)


def find_most_likely_values(samples, ndims):
    """This function finds values by finding the peak value in the probability
    distribution, it also extracts errors by trying to encompass half the
    selected confidence interval on each size"""
    return [(0.7, 0.2, 0.2)]*ndims


def make_prob_plots(samples, energy):
    """this function takes the list of samples and makes histograms of the
    probability distributions of the parameters using matplotlib and writes
    those histograms to the specified directory"""
    ndims = len(samples[0])
    lbls = [r"$a_{%d}$" % i for i in range(ndims)]
    ranges = [(-0.0001, (0.0001 + samples[:, i].max())) for i in range(ndims)]
    quantile_list = np.array([(0.5 - CONFIG["Confidence Interval"] / 2.0), 0.5,
                              (0.5 + CONFIG["Confidence Interval"] / 2.0)])
    for i in range(ndims):
        temp = samples[:, i]
        fig = corner.corner(temp, labels=[lbls[i]], range=[ranges[i]],
                            quantiles=quantile_list, verbose=False)
        # make the probability plot file_name
        fig_file_name = None
        if CONFIG["Prob Plots Directory"][-1] == '/':
            fig_file_name = CONFIG["Prob Plots Directory"] +\
                "A%d_prob_en_%4.1f_a%02d.%s" % (CONFIG["Target A"], energy, i,
                                                CONFIG["Plot Format"])
        else:
            fig_file_name = CONFIG["Prob Plots Directory"] +\
                "/A%d_prob_en_%4.1f_a%02d.%s" % (CONFIG["Target A"], energy,
                                                 i, CONFIG["Plot Format"])
        fig.savefig(fig_file_name, bbox_inches='tight')
        plt.close(fig)


def ln_post_prob(params, cs_lib, struct, bounds):
    """This function calculates the log of the post probability function"""
    # first check if we are outside the resonable parameter range, if so,
    # return negative infinity, which corresponds to a probability of 0
    for i in range(len(bounds)):
        if params[i] < bounds[i][0] or params[i] > bounds[i][1]:
            return -np.inf
    return cs_lib.calculateLnLiklihood(struct, params.ctypes.data)


def gen_walker_starts(gens):
    """This function takes the set of cluster centers and generates the walker
    start positions from them"""
    # store the number of dimensions
    ndims = (1 + CONFIG["Maximum L"])
    # generate the output
    output = [randomize_position(gens[i % CONFIG["Number Walker Generators"]],
                                 ndims)
              for i in range(CONFIG["Number of Walkers"])]
    return output


def randomize_position(gen, ndims):
    """This function takes a generator position and suitable randomizes it"""
    # make an array of randomizers, normally distributed around zero
    rands = CONFIG["Sample Spread"]*np.random.standard_normal(ndims)
    # convert them to fractions of the original value
    randomizer = (np.ones((ndims), dtype=np.float64) - rands)
    # make the randomized positions
    position = (gen*randomizer)
    # now if any of the positions are negative set them to zero
    for i in range(ndims):
        if position[i] < 0.0:
            position[i] = 0.0
    return position


def do_init_fit(start, struct, cs_lib):
    """Performs a fit from the given starting point"""
    # calculate the bounds
    bnds = [(0.0, 1.0/x) for x in CONFIG["EWSR Fractions"]]
    # print call_chi_sq(np.array(start, dtype=np.float64), cs_lib, struct)
    init_params = np.array(start, dtype=np.float64)
    ret_vals = optimize.fmin_l_bfgs_b(call_chi_sq,
                                      init_params, bounds=bnds,
                                      epsilon=1e-03, approx_grad=True,
                                      args=(cs_lib, struct), iprint=0,
                                      factr=10.0)
    return (ret_vals[1], ret_vals[0])


def call_chi_sq(params, cs_lib, struct):
    """calls the chi^2 function in cs_lib"""
    temp = cs_lib.calculateChi(struct, params.ctypes.data)
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


def generate_output_dirs():
    """This function checks for the existence of the directories output is to
    be placed in, if they do not exist, they are created"""
    # test / create the directories for csv files with individial fits
    make_config_fit_csv_dirs()
    # test / create the directories for fit plots
    make_config_fit_plot_dirs()
    # test / create the directory for corner plots
    if not os.path.exists(CONFIG["Corner Plots Directory"]):
        os.makedirs(CONFIG["Corner Plots Directory"])
    # test / create the directory for probability plots
    if not os.path.exists(CONFIG["Prob Plots Directory"]):
        os.makedirs(CONFIG["Prob Plots Directory"])
    # test / create the directory for Markov Chains
    if not os.path.exists(CONFIG["Chain Directory"]):
        os.makedirs(CONFIG["Chain Directory"])
    # test / create the directory for Parameter Plots
    if not os.path.exists(CONFIG["Parameter Plots Directory"]):
        os.makedirs(CONFIG["Parameter Plots Directory"])
    # test / create the directory for the output file
    if not os.path.exists(CONFIG["Parameter Files Directory"]):
        os.makedirs(CONFIG["Parameter Files Directory"])


def make_config_fit_plot_dirs():
    """This function encapsulates the annoying task of making all the sub dir
    names for holding fit plots"""
    temp = CONFIG["Fit Plots Directory"]
    CONFIG["Fit Plot Dirs"] = [copy.deepcopy(temp) for _ in range(12)]
    if temp[-1] != '/':
        for i in range(len(CONFIG["Fit Plot Dirs"])):
            CONFIG["Fit Plot Dirs"][i] += "/"
    # split the dirs first by peak vs percentile
    for i in range(0, 6):
        CONFIG["Fit Plot Dirs"][i] += "Percentiles/"
    for i in range(6, 12):
        CONFIG["Fit Plot Dirs"][i] += "Peaks/"
    # now split them by limited and unlimited
    for i in range(0, 3):
        CONFIG["Fit Plot Dirs"][i] += "Complete/"
        CONFIG["Fit Plot Dirs"][i+6] += "Complete/"
    for i in range(3, 6):
        CONFIG["Fit Plot Dirs"][i] += "Limited/"
        CONFIG["Fit Plot Dirs"][i+6] += "Limited/"
    # now split them by lo_err, fit, hi_err
    for i in range(0, 12, 3):
        CONFIG["Fit Plot Dirs"][i] += "Low_Edge/"
    for i in range(1, 12, 3):
        CONFIG["Fit Plot Dirs"][i] += "Parameters/"
    for i in range(2, 12, 3):
        CONFIG["Fit Plot Dirs"][i] += "High_Edge/"
    for fit_dir in CONFIG["Fit Plot Dirs"]:
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)


def make_config_fit_csv_dirs():
    """This function takes the configuration information and generates the two
    folders that will hold the fit csv files"""
    CONFIG["Fit Csv Dirs"] = [copy.deepcopy(CONFIG["Fits Csv Directory"]),
                              copy.deepcopy(CONFIG["Fits Csv Directory"])]
    if CONFIG["Fits Csv Directory"][-1] != "/":
        for i in range(2):
            CONFIG["Fit Csv Dirs"][i] += "/"
    CONFIG["Fit Csv Dirs"][0] += "percentiles/"
    CONFIG["Fit Csv Dirs"][1] += "peaks/"
    for fit_dir in CONFIG["Fit Csv Dirs"]:
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)


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
                                           error_list).astype(np.float64)
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
    return (interp(angles)/errors).astype(np.float64)


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
    return np.array(output, dtype=np.float64)


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
        dists_list.append(np.array(dist, dtype=np.float64))
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
            output.append([energy, np.array(distribution, dtype=np.float64)])
    print "Exp Data is read in"
    return output


SAMPLES_ERROR = """
WARNING: the number of samples:
%d = ((%d - %d) * %d)
is not large enough to ensure 10 points outside of each error bar location
for the given confidence interval. The minimum number of samples necessary is:
%d = ceiling[ 10 / ((1.0 - %12.10f) / 2.0))

For more information look at the mda_config.py file.
"""


if __name__ == "__main__":
    main()
