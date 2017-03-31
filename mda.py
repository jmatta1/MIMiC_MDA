#!/usr/bin/python
"""This program performs a multipole decomposition analysis of the data and
options given in the configuration file. First it finds starting points using
the BFGS algorithm from a variety of starting positions given in the config
file. Then it takes those starting points, refines them and finds parameter
errors using the Markov Chain Monte Carlo technique, specifically the python
implementation of Goodman & Weare's Affine Invariant MCMC ensemble sampler
which is in the emcee package

This program can be invoked with:
    ./mda.py configuration_file
  or
    python mda.py configuration_file

Parameters
----------
config_file : string
    This is the file that fills out the configuration and run details of the
    MCMC that needs to be run. For details of what needs to be in the config
    file, see config_example.py, which came with the repository, for more
    information

Returns
-------
"""
import os
import ctypes as ct
import multiprocessing
import sys
import math
import copy
import emcee
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import numpy as np
from scipy import interpolate
from scipy import optimize

# firsts we check the command line and grab the module name
if len(sys.argv) != 2:
    print "\nUsage:\n\t./mda.py configuration_file\n\t  or"
    print "\tpython mda.py configuration_file\n"
    sys.exit()
# now we test for file existence
if not os.path.exists(sys.argv[1]):
    print "Error: File {0:s} does not exist".format(sys.argv[1])
    sys.exit()
# now set up and use the execfile function to read the parameters
TEMP_PARAMS = {}
execfile(sys.argv[1], TEMP_PARAMS)
CONFIG = TEMP_PARAMS["CONFIG"]

FLOAT_EPSILON = 1.0e-7
PLOT_FORMAT_LIST = ["svg", "svgz", "pdf", "ps", "eps", "png"]


def main():
    """performs sanity checks on the configuration data and then calls the
    functions that do the work

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Number of Threads', 'Sample Points',
        'Burn-in Points', 'Confidence Interval', 'Number of Walkers',
        'Maximum L', 'EWSR Fractions', 'Start Pts a%d',
        'Number Walker Generators', and 'Plot Format', keys

    Returns
    -------
    """
    # check that the user gave sane information
    # check that they are not requesting greater concurrency than the
    # system supports
    cpu_count = multiprocessing.cpu_count()
    if CONFIG["Number of Threads"] > cpu_count:
        print "\nInvalid number of threads, on this machine it must be: "
        print 'CONFIG["Number of Threads"] <= {0:d}\n'.format(cpu_count)
        sys.exit()
    # check if the user set the sampling high enough for the error bars wished
    num_samples = ((CONFIG["Sample Points"] - CONFIG["Burn-in Points"]) *
                   CONFIG["Number of Walkers"])
    samples_needed = int(math.ceil(10.0 /
                                   ((1.0 - CONFIG["Confidence Interval"]) /
                                    2.0)))
    if samples_needed > num_samples:
        print SAMPLES_ERROR.format(num_samples, CONFIG["Sample Points"],
                                   CONFIG["Burn-in Points"],
                                   CONFIG["Number of Walkers"], samples_needed,
                                   CONFIG["Confidence Interval"])
        sys.exit()
    # make sure that the number of walkers for time series plots does not
    # exceed the number of walkers
    if CONFIG["Number of Walkers"] < CONFIG["Walker Plot Count"]:
        print 'CONFIG["Number of Walkers"] must exceed '\
            'CONFIG["Walker Plot Count"]'
        sys.exit()
    # make certain the user gave enough EWSR fractions for the max L
    num_dists = (1 + CONFIG["Maximum L"])
    if len(CONFIG["EWSR Fractions"]) > num_dists:
        print TOO_MANY_EWSR_ERROR.format(num_dists,
                                         len(CONFIG["EWSR Fractions"]))
        sys.exit()
    elif len(CONFIG["EWSR Fractions"]) < num_dists:
        print TOO_FEW_EWSR_ERROR.format(num_dists,
                                        len(CONFIG["EWSR Fractions"]))
        sys.exit()
    # make sure as many corner plot bins as fit params were supplied
    if len(CONFIG["Corner Plot Bins"]) > num_dists:
        print TOO_MANY_BINS_ERROR.format(num_dists,
                                         len(CONFIG["Corner Plot Bins"]))
        sys.exit()
    elif len(CONFIG["Corner Plot Bins"]) < num_dists:
        print TOO_FEW_BINS_ERROR.format(num_dists,
                                        len(CONFIG["Corner Plot Bins"]))
        sys.exit()
    # check to make certain that there are at least 10 start points
    len_array = [len(CONFIG["Start Pts a{0:d}".format(i)]) for i in
                 range(num_dists)]
    num_cells = 1
    for size in len_array:
        num_cells *= size
    if num_cells < CONFIG["Number Walker Generators"]:
        out_str = NUM_STARTS_ERROR.format(CONFIG["Number Walker Generators"])
        print out_str
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
    that do the initial fitting and then the sampling

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L', and 'Number of Threads' keys

    Returns
    -------
    """
    print STARTUP_MSG
    # read the raw data
    (exp_data, plot_data) = read_row_cs_data_file()
    # now read and subtract the IVGDR data
    ivgdr_info = handle_ivgdr(exp_data)
    # sub_data = ivgdr_info[2]
    # print ivgdr_dists[0], '\n', ivgdr_ewsr[0], '\n', ivgdr_info[2][0]
    # now read the distributions that are used to fit the data
    dists = [[read_dist(elem[0], i) for i in range(CONFIG["Maximum L"] + 1)]
             for elem in exp_data]
    print "Distributions read in"
    # now interpolate the distributions to get the values at the angles in data
    interp_dists = interp_all_dists(dists, exp_data)
    print "Distributions interpolated and prepared for fitting"
    # now get the data divided by errors without angle values
    fit_data = [(exp_en[1][:, 1]/exp_en[1][:, 2]) for exp_en in ivgdr_info[2]]
    print "Experimental data prepared for fitting"
    # calculate the starting parameter sets for initial searches
    start_params = calc_start_params()
    print "Finished calculating parameter starting point list"""
    # now interleave things so we are ready to use pool.map across everything
    interleaved_data = [(exp_data[i][0], fit_data[i], interp_dists[i],
                         start_params) for i in range(len(fit_data))]
    print "Data is interleaved"
    generate_output_dirs()
    mp_pool = multiprocessing.Pool(processes=CONFIG["Number of Threads"])
    print MDA_START_MSG.format(CONFIG["Number of Threads"])
    fitted_data = mp_pool.map(fit_and_mcmc, interleaved_data)
    # single threaded version for debugging
    # fitted_data = map(fit_and_mcmc, interleaved_data)
    # write the individual fits to csv files
    parameters = [dat[0] for dat in fitted_data]
    diag_data = [dat[1] for dat in fitted_data]
    if CONFIG["Generate Fit CSVs"]:
        print "Writing fit files"
        write_fits(plot_data, dists, parameters, ivgdr_info)
    else:
        print "Skipping fit files"
    # make the fit plots
    if CONFIG["Generate Fit Plots"]:
        print "Writing fit plots"
        make_fit_plots(plot_data, dists, parameters, ivgdr_info)
    else:
        print "Skipping fit plots"
    # write the two parameter sets
    energy_set = [val[0] for val in exp_data]
    print "Writing parameter sets"
    write_param_sets(parameters, energy_set)
    # write the parameter plots
    if CONFIG["Generate Parameter Plots"]:
        print "Writing parameter plots"
        write_param_plots(parameters, energy_set)
    else:
        print "Skipping parameter plots"
    # write the diagnostic information
    print "Writing diagnostic information"
    write_diagnostic_csv(diag_data, exp_data)
    # and now we are completely done
    print "MDA and output complete"


def write_diagnostic_csv(diag_data, exp_dat):
    """This function takes the diagnostic data generated in the fit process and
    writes it to a csv for assessment

     Parameters
    ----------
    diag_data : list
        The a list of the diagnostic data for each fit

    exp_dat : list
        The experimental data for each energy that was fitted

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L', 'Calc AutoCorr', and 'Target A' keys

    Returns
    -------
    """
    num_params = CONFIG["Maximum L"] + 1
    # for each angular distribution, calculate how many points there are
    num_pts = [len(dat[1]) for dat in exp_dat]
    # generate the output file
    file_name = "A{0:d}_diagnostics.csv".format(CONFIG["Target A"])
    file_path = os.path.join(CONFIG["Parameter Files Directory"], file_name)
    outfile = open(file_path, 'w')
    # now write the csv header
    outfile.write("Ex Energy, , Percentile Chi^2, Peak Chi^2, Num Fitted ")
    outfile.write("Data Points, Num Params, Num DoF, , Acceptance Fraction")
    if CONFIG['Calc AutoCorr']:
        outfile.write(', , NumIndSamples, Dist Error, ')
        for i in range(num_params):
            outfile.write(", a{0:d} AutoCorr Time".format(i))
    outfile.write('\n')
    num_samps = CONFIG["Number of Walkers"]*(CONFIG["Sample Points"] -
                                             CONFIG["Burn-in Points"])
    # now move through the fit data and write out all the columns
    # (acor_time, chis, accept_frac)
    for i, diag in enumerate(diag_data):
        outfile.write('{0:f}, , {1:f}, {2:f}'.format(exp_dat[i][0], diag[1][0],
                                                     diag[1][1]))
        dof = (num_pts[i] - num_params)
        outfile.write(', {0:d}, {1:d}, {2:d}'.format(num_pts[i], num_params,
                                                     dof))
        outfile.write(', , {0:f}'.format(diag[2]))
        if CONFIG['Calc AutoCorr']:
            indsamps = float(num_samps)/max(diag[0])
            err = 1.0/math.sqrt(indsamps)
            outfile.write(', {0:f},  {1:f}, '.format(indsamps, err))
            for accfrac in diag[0]:
                outfile.write(", {0:f}".format(accfrac))
        outfile.write('\n')
    outfile.close()


def write_param_plots(parameters, energy_set):
    """This function takes the generated parameters and makes the plots for
    each L of the parameters

    Parameters
    ----------
    parameters : list
        The full set of parameter sets, for all the runs, from both parameter
        set finding methodologies

    energy_set : list of floats
        The list of excitation energies for each run

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L', 'Param Plot Dirs', 'Target A', and
        'Plot Format' keys

    Returns
    -------
    """
    # loop through each set of parameters
    for i in range((CONFIG["Maximum L"]+1)):
        # first split the data into the two types
        perc_data = [pset[0][i] for pset in parameters]
        peak_data = [pset[1][i] for pset in parameters]
        # calculate the two file names
        fmt_str = "A{0:d}_L{1:d}_{2:s}_parameters.{3:s}"
        file_name = fmt_str.format(CONFIG["Target A"], i, "percentile",
                                   CONFIG["Plot Format"])
        perc_path = os.path.join(CONFIG["Param Plot Dirs"][0], file_name)
        file_name = fmt_str.format(CONFIG["Target A"], i, "peak",
                                   CONFIG["Plot Format"])
        peak_path = os.path.join(CONFIG["Param Plot Dirs"][1], file_name)
        make_param_plot(perc_path, perc_data, energy_set, i)
        make_param_plot(peak_path, peak_data, energy_set, i)


def make_param_plot(path, params, energy_set, l_value):
    """This takes a set of parameters for a given L, the energies they are from
    and generates a plot of those parameters. It then writes that plot to the
    specified path

    Parameters
    ----------
    path : string
        The file path for this plot

    params : list of floats
        One set of parameters for each run

    energy_set : list of floats
        The list of excitation energies for each run

    l_value : int
        The orbital angular momentum of the GR this plot is for

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Plot Height', 'Plot Width', and 'Plot DPI' keys

    FLOAT_EPSILON : float
        A very small value that is the threshold for "two float are the same"

    Returns
    -------
    """
    pt_x_vals = np.array(energy_set)
    param_array = np.array(params)
    pt_y_vals = param_array[:, 0]
    pt_e_vals = [param_array[:, 1], param_array[:, 2]]
    hi_vals = pt_y_vals + pt_e_vals[0]
    # make the figure and stuff
    fig, axes = plt.subplots()
    # set up the axes
    axes.set_yscale('linear')
    axes.set_xscale('linear')
    # plot the data
    axes.errorbar(pt_x_vals, pt_y_vals, yerr=pt_e_vals, fmt="ko",
                  label=r"$Exp$", markersize=2.0)
    # set the axis limits
    axes.set_xlim((pt_x_vals.min() - 1.0), (pt_x_vals.max() + 1.0))
    y_max = 1.2 * hi_vals.max()
    if y_max < FLOAT_EPSILON:
        y_max = 0.01
    axes.set_ylim(0.0, y_max)
    # label the axes
    axes.set_xlabel('Excitation Energy (MeV)')
    axes.set_ylabel(r'$a_{{{0:d}}}$'.format(l_value))
    fig.suptitle(r'MDA Results for L={0:d}'.format(l_value))
    # make the legend
    # legend = axes.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1)
    legend = axes.legend(loc='upper left', ncol=1)
    legend.get_frame().set_facecolor("white")
    # save and close the figure
    fig.set_size_inches(CONFIG["Plot Height"], CONFIG["Plot Width"])
    # fig.savefig(path, additional_artists=[legend], bbox_inches='tight')
    fig.savefig(path, bbox_inches='tight', dpi=CONFIG["Plot DPI"])
    plt.close(fig)


def write_param_sets(parameters, energy_set):
    """This function takes the generated parameters and writes the sets from
    percentiles and the sets from peak find to seperate files

    Parameters
    ----------
    parameters : list
        Both types of parameter sets from each run and their associated errors

    energy_set : list of floats
        The list of excitation energies for each run

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Parameter Files Directory' and 'Target A' keys

    Returns
    -------
    """
    # first split the data into the two types
    perc_data = [pset[0] for pset in parameters]
    peak_data = [pset[1] for pset in parameters]
    # generate the file names
    fmt_str = "A{0:d}_{1:s}_parameters.csv"
    file_name = fmt_str.format(CONFIG["Target A"], "percentile")
    perc_path = os.path.join(CONFIG["Parameter Files Directory"], file_name)
    file_name = fmt_str.format(CONFIG["Target A"], "peak")
    peak_path = os.path.join(CONFIG["Parameter Files Directory"], file_name)
    # now call the function that writes a parameter set
    write_horizontal_param_file(perc_path, perc_data, energy_set, "percentile")
    write_horizontal_param_file(peak_path, peak_data, energy_set, "peak find")


def write_horizontal_param_file(path, params, energy_set, ptype):
    """This function takes a set of parameters and their corresponding energies
    and writes the data to the specified path

    Parameters
    ----------
    path : string
        The file path for this csv

    params : list of floats
        A parameter set and its errors

    energy_set : list of floats
        The list of excitation energies for each run

    ptype : string
        The type of fit used to extract the parameter values

    Global Parameters
    -----------------

    Returns
    -------
    """
    # open the file for writing
    out_file = open(path, 'w')
    # write the header
    out_file.write(generate_param_file_header(ptype))
    # write each line of parameters and energies
    for i, param in enumerate(params):
        out_file.write(gen_param_file_line(energy_set[i], param))
    # close the file we wrote
    out_file.close()


def gen_param_file_line(energy, params):
    """This function takes the excitation energy and the parameters for that
    energy and generates the string to be written to the file

    Parameters
    ----------
    energy : float
        The excitation energy of the parameter set

    params : list of floats
        A parameter set and its errors

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L' key

    Returns
    -------
    out_str : string
        line of parameters written out in the csv format
    """
    out_str = ("{0:f}, , ".format(energy))
    for i in range((CONFIG["Maximum L"]+1)):
        out_str += "{0:f}, {1:f}, {2:f}, , ".format(*params[i])
    out_str += "\n"
    return out_str


def generate_param_file_header(ptype):
    """This function uses config information to generate an appropriate header
    for the parameter file

    Parameters
    ----------
    ptype : string
        The type of fit used to extract the parameter values

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L' key

    Returns
    -------
    out_str : string
        Header for the CSV file that will be output for the summary of param
        values
    """
    out_str = "MDA parameters extracted with {0:s} techniques\n".format(ptype)
    hrow1 = "Excitation, , "
    hrow2 = "Energy, , "
    for i in range((CONFIG["Maximum L"]+1)):
        hrow1 += "a{0:d}, a{0:d}, a{0:d}, , ".format(i)
        hrow2 += "Value, Pos. Err, Neg. Err., , "
    out_str += hrow1
    out_str += "\n"
    out_str += hrow2
    out_str += "\n"
    return out_str


def make_fit_plots(data, dists, parameters, ivgdr_info):
    """This function takes everything and generates the plots for individual
    fits at each energy

    Parameters
    ----------
    data : list of lists
        List of excitation energies and associated experimental data for each
        run

    dists : list of lists of numpy arrays
        List of DWBA distributions for each run

    parameters : list of lists of floats
        List of parameter sets derived with both methods for each run

    ivgdr_info : list of lists
        List of IVGDR distributions and ewsr fractions for each run

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Subtract IVGDR', 'Fit Plot Dirs', 'Target A',
        'Plot Format', 'Fit Plot L Limit', and 'Plot IVGDR in Fits' keys

    Returns
    -------
    """
    # first make the directory names
    legend = [r"$Fit$"]
    for i in range(len(parameters[0][0])):
        legend.append(r"$l_{{{0:d}}}$".format(i))
    # loop through each set of data and distributions
    for i, dat in enumerate(data):
        # extract the things pertinet to this set from the arguments
        energy = dat[0]
        dist_set = copy.deepcopy(dists[i])
        param_set = list(copy.deepcopy(parameters[i]))
        exp_points = dat[1]
        # handle ivgdr
        if CONFIG["Subtract IVGDR"]:
            dist_set.append(ivgdr_info[0][i])
            param_set.append((ivgdr_info[1][i], 0.0, 0.0))
            legend.append(r"$l_{-1}$")
        # loop through the plots to generate
        for ind in range(2):  # peaks or percentiles?
            # choose our parameter set and generate the three sets of fits
            pset = gen_param_sets_for_fit_plot(param_set[ind])
            for k in range(3):  # params-lo_errs, params, params+hi_errs
                sc_dists = gen_fit_dists(pset[k], dist_set)
                for j in range(2):  # full or limitted
                    fmt_str = "A{0:d}_E{1:05.2f}.{2:s}"
                    file_name = fmt_str.format(CONFIG["Target A"], energy,
                                               CONFIG["Plot Format"])
                    dind = 6 * ind + 3 * j + k
                    plt_path = os.path.join(CONFIG["Fit Plot Dirs"][dind],
                                            file_name)
                    # now decide how much of the distributions to call the
                    # gen fit plot function on
                    if j == 1:
                        gen_fit_plot(exp_points, energy,
                                     sc_dists[:(CONFIG["Fit Plot L Limit"]+2)],
                                     legend[:(CONFIG["Fit Plot L Limit"]+2)],
                                     plt_path)
                    elif (CONFIG["Subtract IVGDR"] and
                          not CONFIG["Plot IVGDR in Fits"]):
                        gen_fit_plot(exp_points, energy,
                                     sc_dists[:(len(sc_dists)-1)],
                                     legend[:(len(sc_dists)-1)],
                                     plt_path)
                    else:
                        gen_fit_plot(exp_points, energy,
                                     sc_dists[:len(sc_dists)],
                                     legend[:len(sc_dists)],
                                     plt_path)


def gen_param_sets_for_fit_plot(params):
    """This function takes a single set of parameters and generates three sets
    of parameters, the first, decreased by the lower error bar, the second
    equal to the parameter fit value, and the third, increased by the upper
    error bar

    Parameters
    ----------
    params : list of tuples
        This is the list of parameter values of their errors bars, either
        derived from the peak method or the quantile method

    Global Parameters
    -----------------

    Returns
    -------
    parameter_bounds : list of tuples
        This is a list of three parameter sets, the first is all params reduced
        by one error bar, the second is the params, the third is the params
        increased by one error bar
    """
    temp = [[(vals[0] - vals[2]) for vals in params],
            [vals[0] for vals in params],
            [(vals[0] + vals[1]) for vals in params]]
    for pset in temp:
        for param in pset:
            if param < 0.0:
                print "Got a negative parameter when subtracting errors from"
                print "parameters for fit plots, this should not be possible"
                sys.exit()
                param = 0.0
    return temp


def gen_fit_plot(points, energy, dists, legends, plot_name):
    """This function takes the list of points with errors in the points var
    and the list of distributions to plot in the dists var and generates a
    nicely formatted matplotlib plot displaying them

    Parameters
    ----------
    points : list of floats
        Experimental data for one energy

    energy : float
        Excitation energy of the data

    dists : list of numpy arrays
        List of dwba distributions

    legends : list of strings
        List of names of the points and distributions

    plot_name : string
        output path of the plot

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Subtract IVGDR', 'Fit Plot Dirs', 'Target A',
        'Plot Format', 'Fit Plot L Limit''Plot Height', 'Plot Width',
        'Plot DPI', and 'Plot IVGDR in Fits' keys

    Returns
    -------
    """
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
    for i, dis in enumerate(dists):
        axes.plot(dis[:, 0], dis[:, 1],
                  line_styles[i % len(line_styles)], label=legends[i])
    # set the scale of the x axis
    axes.set_xlim(0.0, math.ceil(pt_x_vals.max()))
    # set the scale of the y axis
    (ymin_val, ymax_val) = find_y_extrema(pt_y_vals.max(), dists,
                                          math.ceil(pt_x_vals.max()))
    axes.set_ylim(ymin_val, ymax_val)
    # label the axes
    axes.set_xlabel(r'Lab Angle $(^{\circ{}})$')
    axes.set_ylabel(r'$(\partial^2 \sigma)/(\partial \Omega \partial E)$'
                    ' ($mb/(sr*MeV)$)')
    fig.suptitle(r"MDA Fit for E$_x$={0:4.2f} MeV".format(energy))
    # make the legend
    legend = axes.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1)
    legend.get_frame().set_facecolor("white")
    # save and close the figure
    fig.set_size_inches(CONFIG["Plot Height"], CONFIG["Plot Width"])
    fig.savefig(plot_name, additional_artists=[legend], bbox_inches='tight',
                dpi=CONFIG["Plot DPI"])
    plt.close(fig)


def find_y_extrema(data_max, dists, xmax):
    """This function scans the distributions provided searching for the minimum
    value with angle less than xmax, it also looks to find the maximum value
    (be it in a distribution or the data maximum it then returns
    (10^(floor(log(ymin))), 10^(ceiling(log(ymax)))

    Parameters
    ----------
    data_max : float
        The maximum cross-section in the experimental data

    dists : list of numpy arrays
        The set of distributions for this MCMC

    xmax : float
        the maximum angle that will be plotted

    Global Parameters
    -----------------

    Returns
    -------
    plot_min : float
        A suggested minimum value for plot y axes

    plot_max : float
        A suggested maximum value for plot y axes
    """
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
    them to a nicely formatted csv file for usage later

    Parameters
    ----------
    data : list of lists
        Each element of the list is one MCMC run, each sublist has the ex
        energy as the first element, the set of experimental data points as the
        second element

    dists : list of lists
        Each sub list contains the numpy arrays with dwba distributions for one
        run of the MCMC

    parameters : list of lists
        Each sub list contains two lists. The first list in the sublist is the
        parameters and their errors as derived from the pure percentile method,
        the second is the parameters and their errors as derived from the peak
        method

    ivgdr_info : list of lists
        Each sub list is from one run and contains the following: the IVGDR
        distribution as its first element, and the ivgdr EWSR weight as its
        second

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Subtract IVGDR', 'Fit Csv Dirs', and 'Target A'
        keys

    Returns
    -------
    """
    # first split the data up into individual runs
    for i, dat in enumerate(data):
        energy = copy.deepcopy(dat[0])
        points = copy.deepcopy(dat[1])
        dist_set = copy.deepcopy(dists[i])
        perc_set = copy.deepcopy(parameters[i][0])
        peak_set = copy.deepcopy(parameters[i][1])
        # test if we are subtracting the IVGDR
        if CONFIG["Subtract IVGDR"]:
            dist_set.append(ivgdr_info[0][i])
            perc_set.append((ivgdr_info[1][i], 0.0, 0.0))
            peak_set.append((ivgdr_info[1][i], 0.0, 0.0))
        # calculate the file names
        file_name = "A{0:d}_E{1:5.2f}.csv".format(CONFIG["Target A"], energy)
        file_paths = [os.path.join(CONFIG["Fit Csv Dirs"][0], file_name),
                      os.path.join(CONFIG["Fit Csv Dirs"][1], file_name)]
        # write the percentile fit
        write_fit_csv(file_paths[0], points, perc_set, dist_set, energy)
        # write the peak find fit
        write_fit_csv(file_paths[1], points, peak_set, dist_set, energy)


def write_fit_csv(path, points, pset, dist_set, energy):
    """This function takes a file path, a set of experimental data, a parameter
    set, a set of distributions and the ex energy of the distribution and
    writes a nicely formated csv file with all the information in it

    Parameters
    ----------
    path : string
        the path to the file that will hold this csv data

    points : list of floats
        the experimental data the MCMC was performed on

    pset : list of floats
        the parameter set and their errs (be it from peak finding or quantiles)

    dist_set : list of numpy arrays
        The set of dwba distributions that were scaled and fit to the exp data
        when the MCMC was performed

    energy : float
        The excitation energy of the nucleus that the MCMC was run for

    Global Parameters
    -----------------

    Returns
    -------
    """
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
    for i, dis in enumerate(dist_set):
        append_data_column_to_fit_csv(dis[:, 1], csv_list)
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
    """This function takes a column of data and appends it to the csv list

    Parameters
    ----------
    data : list of floats
        the set of data, be it angles, distribution values, etc to append to
        the ends of the csv lines

    csv_list : list of strings
        The list of strings that when complete and written to a file produces
        a csv with the appropriate formatting

    Global Parameters
    -----------------

    Returns
    -------
    """
    for i, sub_csv_list in enumerate(csv_list):
        if i < len(data):
            sub_csv_list += "{0:f}, ".format(data[i])
        else:
            sub_csv_list += ", "


def append_parameters_to_fit_csv(pset, csv_list):
    """This function appends a parameter set to the csv list

    Parameters
    ----------
    pset : list of floats
        the parameter set to append to the lines of the fit csv

    csv_list : list of strings
        The list of strings that when complete and written to a file produces
        a csv with the appropriate formatting

    Global Parameters
    -----------------

    Returns
    -------
    """
    # first generate the list of names
    name_list = []
    for i in range(len(pset)):
        if CONFIG["Subtract IVGDR"] and i == (len(pset)-1):
            name_list.append("a-1")
        else:
            name_list.append("a{0:d}".format(i))
    # now append the names and error bars
    fmt_str = "{0:s}, {1:f}, {2:f}, {3:f}, "
    for i, csv_sublist in enumerate(csv_list):
        if i < len(pset):
            csv_sublist += fmt_str.format(name_list[i], pset[i][0],
                                          pset[i][1], pset[i][2])
        else:
            csv_sublist += ", , , , "


def append_str_to_fit_csv(str_to_append, csv_list):
    """This function appends the given string to the csv_list

    Parameters
    ----------
    str_to_append : string
        the string to append to every string in the csv_list

    csv_list : list of strings
        The list of strings that when complete and written to a file produces
        a csv with the appropriate formatting

    Global Parameters
    -----------------

    Returns
    -------
    """
    for csv_sublist in csv_list:
        csv_sublist += str_to_append


def append_exp_data_to_fit_csv(points, csv_list):
    """This function takes an experimental data set and the csv list and writes
    the data to the csv list along with a set of appropriately blank lines

    Parameters
    ----------
    points : list of floats
        the list of angles, data, and erroros to be appended to each of the
        lines in csv_list

    csv_list : list of strings
        The list of strings that when complete and written to a file produces
        a csv with the appropriate formatting

    Global Parameters
    -----------------

    Returns
    -------
    """
    fmt_str = "{0:f}, {1:f}, {2:f}, "
    for i, csv_sublist in enumerate(csv_list):
        if i < len(points):
            csv_sublist += fmt_str.format(points[i][0], points[i][1],
                                          points[i][2])
        else:
            csv_sublist += " , , , "


def gen_csv_title_and_headings(energy):
    """This function takes the excitation energy and returns a string with the
    csv title, and column headings

    Parameters
    ----------
    energy : float
        The excitation energy at which this decomposition was performed

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L' key

    Returns
    -------
    out_str : string
        The string wth the first few header lines of the csv file output
    """
    # first generate the title string
    out_str = "Fit Information for Ex =, {0:4.1f}\n".format(energy)
    # now generate the two rows of column headings
    # first the easy column headings
    row1 = "Exp., Exp., Exp., , Param, Param., Pos., Neg., , Dist., "
    row2 = "Angle, CS, Error, , Name, Value, Error, Error, , Angle, "
    # write the headings for the distributions
    for i in range((CONFIG["Maximum L"]+1)):
        row1 += ", Scaled, "
        row2 += "l{0:d}, l{0:d}, ".format(i)
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
    scaled distributions

    Parameters
    ----------
    params : list
        This is the parameter set to be used to scale the dwba distributions

    dists : list of numpy arrays
        This is the list of dwba distributions to be scaled

    Global Parameters
    -----------------

    Returns
    -------
    output : list of lists
        array of NxP size where N is the number of points in the distribution
        and P is the number of distributions used in the MCMC plus two, then
        column 0 is the list of angles, column 1 is the sum of the scaled
        distributions and columns 2 through P-1 are the scaled distributions
    """
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
    and then performs the MCMC using those fits

    Parameters
    ----------
    data_tuple : tuple
        The first element of the tuple is the energy the mcmc is being
        performed on, the second element of the tuple is the experimental data,
        the third element of the tuple is the interpolated distributions to
        use, the final element of the tuple is the set of start point for the
        initial fitting

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Shared Lib Path', 'Number Walker Generators',
        'EWSR Fractions', 'Number of Walkers', and 'Sample Points' keys

    Returns
    -------
    values : list of lists
        The first sub list contains the parameters and their errors as
        extracted via the pure quantiles method, the second sublist is the
        parameters and their errors as extracted by the peak finding then
        quantiles method, this is the output as yielded by the
        perform_sample_manips function
    """
    # first unpack the tuple
    energy = data_tuple[0]
    fit_data = data_tuple[1]
    interp_dists = data_tuple[2]
    start_points = data_tuple[3]
    print "Starting work on Energy", energy, "MeV"
    # now load the shared library
    cs_lib = prepare_shared_object()
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
    # now do the rest, all of which involves calculating things from samples
    # or the things calculated from samples
    ret_vals = perform_sample_manips(sampler, ndims, energy, cs_lib, struct)
    # delete the struct
    cs_lib.freeMdaStruct(struct)
    # return our return values
    return ret_vals


def prepare_shared_object():
    """This function loads an instance of the shared object and sets the
    restypes and argtypes for the functions of interest

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Shared Lib Path' key

    Returns
    -------
    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines
    """
    cs_lib = ct.cdll.LoadLibrary(CONFIG["Shared Lib Path"])
    # tell python the function argument and return types for all the functions
    # that the library exports, even if we do not use all of them
    # restype for makeMdaStruct
    cs_lib.makeMdaStruct.restype = ct.c_void_p
    # argtypes for makeMdaStruct
    cs_lib.makeMdaStruct.argtypes = [ct.c_int, ct.c_int]
    # arg types for freeMdaStruct (it returns void)
    cs_lib.freeMdaStruct.argtypes = [ct.c_void_p]
    # arg types for setMdaData (it returns void)
    cs_lib.setMdaData.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double)]
    # arg types for setMdaDist (it returns void)
    cs_lib.setMdaDist.argtypes = [ct.c_void_p, ct.c_int,
                                  ct.POINTER(ct.c_double)]
    # restype for calculateChi
    cs_lib.calculateChi.restype = ct.c_double
    # argtypes for calculateChi
    cs_lib.calculateChi.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double)]
    # restype for calculateLnLiklihood
    cs_lib.calculateLnLiklihood.restype = ct.c_double
    # argtypes for calculateLnLiklihood
    cs_lib.calculateLnLiklihood.argtypes = [ct.c_void_p,
                                            ct.POINTER(ct.c_double)]
    # restype for calculateLnLiklihoodResids
    cs_lib.calculateLnLiklihoodResids.restype = ct.c_double
    # argtypes for calculateLnLiklihoodResids
    cs_lib.calculateLnLiklihoodResids.argtypes = [ct.c_void_p,
                                                  ct.POINTER(ct.c_double),
                                                  ct.POINTER(ct.c_double)]
    return cs_lib


def perform_sample_manips(sampler, ndims, energy, cs_lib, struct):
    """this function takes the mcmc sampler, grabs the chains from it, and
    generates plots and percentiles and most likely values from the sampled
    points

    Parameters
    ----------
    sampler : mcmc ensemble sampler from emcee
        the ensemble sampler used to perform the mcmc

    ndims : int
        The number of parameters that were sampled

    energy : float
        the excitation energy that this mcmc was carried out for

    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines

    struct : ctypes.c_void_ptr
        Pointer to the struct for calculating chi^2 and liklihoods, etc

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Number of Walkers', 'Sample Points',
        'Burn-in Points', 'Save Chain Data', 'Chain Directory', 'Target A',
        'Confidence Interval', 'Generate Corner Plots', 'Corner Plot Samples',
        'Corner Plots Directory', 'Plot Height', 'Plot Width', 'Plot DPI',
        'Corner Default Range', 'Calc AutoCorr', 'ACorr WindSize',
        'ACorr Use FFT', 'Generate Walker Plots',  and 'Plot Format' keys

    Returns
    -------
    values : list of lists
        The first sub list contains the parameters and their errors as
        extracted via the pure quantiles method, the second sublist is the
        parameters and their errors as extracted by the peak finding then
        quantiles method
    """
    # retrieve the samples
    num_samples = (CONFIG["Number of Walkers"] * (CONFIG["Sample Points"] -
                                                  CONFIG["Burn-in Points"]))
    if CONFIG["Generate Walker Plots"]:
        print "Making Time Series plots for", energy, "MeV"
        gen_time_series_plots(sampler, ndims, energy)
    else:
        print "Skipping Time Series plots for", energy, "MeV"
    # save the samples to the disk in the sp
    if CONFIG["Save Chain Data"]:
        print "Saving MCMC samples for", energy, "MeV"
        base_name = "A{0:d}_chain_E{1:5.2f}.npz".format(CONFIG["Target A"],
                                                        energy)
        chain_file_name = os.path.join(CONFIG["Chain Directory"], base_name)
        np.savez_compressed(chain_file_name, sampler.chain)
        print "Done saving MCMC samples for", energy, "MeV"
    # get a copy of the samples reshaped to all walkers merged together
    print "Reshaping sample chain for", energy, "MeV"
    samples = sampler.chain[:, CONFIG["Burn-in Points"]:, :].reshape((
        num_samples, ndims))
    # extract the error bars
    print "Extracting parameter values for", energy, "MeV"
    quantile_list = [(0.5 - CONFIG["Confidence Interval"] / 2.0), 0.5,
                     (0.5 + CONFIG["Confidence Interval"] / 2.0)]
    values = calc_param_values(samples, quantile_list, ndims)
    peak_vals = [param[0] for param in values[1]]
    # make the probability plots
    if CONFIG["Generate Probability Plots"]:
        print "Generating probability plots for", energy, "MeV"
        make_prob_plots(samples, energy, peak_vals)
    else:
        print "Skipping probability plots for", energy, "MeV"
    # make the corner plot
    if CONFIG["Generate Corner Plots"]:
        print "Commencing corner plot creation for", energy, "MeV"
        lbls = [r"$a_{{{0:d}}}$".format(i) for i in range(ndims)]
        ranges = [(0.00, (0.0001 + samples[:, i].max())) for i in
                  range(ndims)]
        fig = None
        sample_set = None
        if CONFIG["Corner Plot Samples"] >= num_samples:
            sample_set = samples
        else:
            # randomize the sample array and then extract the first chunk of it
            np.random.shuffle(samples)
            sample_set = samples[0:CONFIG["Corner Plot Samples"]]
        if CONFIG["Corner Default Range"]:
            fig = corner.corner(sample_set, labels=lbls,
                                bins=CONFIG["Corner Plot Bins"],
                                quantiles=quantile_list, truths=peak_vals,
                                verbose=False)
        else:
            fig = corner.corner(sample_set, labels=lbls, range=ranges,
                                bins=CONFIG["Corner Plot Bins"],
                                quantiles=quantile_list, truths=peak_vals,
                                verbose=False)
        # make the corner plot file_name
        fmt_str = "A{0:d}_corner_E{1:05.2f}.{2:s}"
        file_name = fmt_str.format(CONFIG["Target A"], energy,
                                   CONFIG["Plot Format"])
        fig_file_name = os.path.join(CONFIG["Corner Plots Directory"],
                                     file_name)
        fig.set_size_inches(CONFIG["Plot Height"], CONFIG["Plot Width"])
        fig.savefig(fig_file_name, bbox_inches='tight', dpi=CONFIG["Plot DPI"])
        plt.close(fig)
        print "Done creating corner plot for", energy, "MeV"
    else:
        print "Skipping corner plot creation for", energy, "MeV"
    acor_time = []
    if CONFIG["Calc AutoCorr"]:
        # calculate the autocorrelation time
        print "Calculating the autocorrelation for", energy, "MeV"
        acor_time = sampler.get_autocorr_time(c=CONFIG["ACorr WindSize"],
                                              fast=CONFIG["ACorr Use FFT"])
    else:
        print "Skipping calculation of the autocorrelation for", energy, "MeV"
    # calculate the chi^2 for the percentile and peak best fits
    print "Calculating fit chi^2 for", energy, "MeV"
    chis = calculate_fit_chis(cs_lib, struct, values)
    # calculate the average acceptance fraction of the mcmc, should be between
    # 0.2 and 0.5
    accept_frac = np.mean(sampler.acceptance_fraction)
    # return the point and the errors
    return (values, (acor_time, chis, accept_frac))


def calculate_fit_chis(cs_lib, struct, values):
    """This function takes the fitted values, percentile and peak, and
    calculates the sum of squared residuals for those parameter sets

    Parameters
    ----------
    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines

    struct : ctypes.c_void_ptr
        Pointer to the struct for calculating chi^2 and liklihoods, etc

    values : tuple of lists
        Tuple of lists and parameter values for peak and percentile methods

    Global Parameters
    -----------------

    Returns
    -------
    perc_chi : float
        The chi^2 of the point with the parameter values calculated from the
        percentiles of the distribution

    peak_chi : float
        The chi^2 of the point with the parameter values calculated from the
        peaks of the distribution
    """
    params = np.array([param[0] for param in values[0]])
    perc_chi = call_chi_sq(params, cs_lib, struct)
    params = np.array([param[0] for param in values[1]])
    peak_chi = call_chi_sq(params, cs_lib, struct)
    return (perc_chi, peak_chi)


def gen_time_series_plots(sampler, ndims, energy):
    """This function generates the time series plots of the walker parameter
    states

    Parameters
    ----------
    sampler : mcmc ensemble sampler from emcee
        the ensemble sampler used to perform the mcmc

    ndims : int
        The number of parameters that were sampled

    energy : float
        the excitation energy that this mcmc was carried out for

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Time Plot Dirs', 'Plot Format', 'Sample Points',
        'Plot Height', 'Plot Width', 'Plot DPI', and 'Walker Plot Count' keys

    Returns
    -------
    """
    fmt_string = "tSeries_E{0:05.2f}_L{1:d}.{2:s}"
    xvals = np.arange(0, CONFIG["Sample Points"])
    samples = sampler.chain
    for i in range(ndims):
        # make the plot name
        fig_name = fmt_string.format(energy, i, CONFIG["Plot Format"])
        fig_file_name = os.path.join(CONFIG["Time Series Directory"],
                                     fig_name)
        fig, axes = plt.subplots()
        # add the first "Walker Plot Count" walkers to the time series plots
        for j in range(CONFIG["Walker Plot Count"]):
            axes.plot(xvals, samples[j, :, i], color='k')
        fig.set_size_inches(CONFIG["Plot Height"], CONFIG["Plot Width"])
        fig.savefig(fig_file_name, bbox_inches='tight', dpi=CONFIG["Plot DPI"])
        plt.close(fig)


def calc_param_values(samples, quantile_list, ndims):
    """This function calculates the parameter values and error bars using both
    peak finding and percentiles

    Parameters
    ----------
    samples : numpy array
        the full set of parameter samples from the MCMC

    quantile_list : numpy array
        An array with the quantiles for lower error bar, value, and high error
        bar, where value is 0.5 and if you want ~1sigma error bars the low and
        high would be 0.16 and 0.84 respectively

    ndims : int
        the number of parameters

    Global Parameters
    -----------------

    Returns
    -------
    points : list
        A list of the parameter values and their errors derived only from
        the quantiles

    peaks : list
        A list of the parameter values and their errors derived from finding
        the distribution peak
    """
    percentile_list = [100.0*q for q in quantile_list]
    points = [(v[1], v[2]-v[1], v[1]-v[0]) for v in
              zip(*np.percentile(samples, percentile_list, axis=0))]
    peaks = find_most_likely_values(samples, ndims)
    return (points, peaks)


def find_most_likely_values(samples, ndims):
    """This function finds values by finding the peak value in the probability
    distribution, it also extracts errors by trying to encompass half the
    selected confidence interval on each side

    Parameters
    ----------
    samples : numpy array
        the full set of parameter samples from the MCMC

    ndims : int
        the number of parameters

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Num Bins' and 'Confidence Interval'  keys

    FLOAT_EPSILON : float
        A very small value that is the threshold for "two float are the same"

    Returns
    -------
    output : list
        A list of the peak parameter values and their error bars
    """
    output = []
    last_index = (len(samples[:, 0])-1)
    # loop through each dimension
    for i in range(ndims):
        # get a sorted list of the parameters
        values = np.sort(samples[:, i])
        if FLOAT_EPSILON > abs(values[last_index] - values[0]):
            # the max and min are the same then there is no need for more
            output.append((values[0], 0.0, 0.0))
            continue
        hist = np.histogram(values, bins=CONFIG["Num Bins"])
        ind = np.argmax(hist[0])
        # find the centroid of the peak bin by looking at the bin edges
        # peak_centroid = (hist[1][ind] + hist[1][ind+1])/2.0
        quantile = 0.0
        for j in range(ind + 1):
            quantile += float(hist[0][j])
        quantile /= float(len(values))
        # find the percentile of the peak
        peak_percentile = (100.0 * quantile)
        lo_p = peak_percentile - 100.0 * (CONFIG["Confidence Interval"] / 2.0)
        if lo_p < 0.0:
            lo_p = 0.0
        hi_p = peak_percentile + 100.0 * (CONFIG["Confidence Interval"] / 2.0)
        if hi_p > 100.0:
            hi_p = 100.0
        vals = np.percentile(values, (lo_p, peak_percentile, hi_p))
        output.append((vals[1], vals[2]-vals[1], vals[1]-vals[0]))
    return output


def make_prob_plots(samples, energy, peak_vals):
    """this function takes the list of samples and makes histograms of the
    probability distributions of the parameters using matplotlib and writes
    those histograms to the specified directory

    Parameters
    ----------
    samples : numpy array
        the full set of parameter samples from the MCMC

    energy : float
        the excitation energy for which this MCMC was carried out

    peak_vals : array of floats
        The values of the peaks of the probability distribution

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Plot Height', 'Plot Width', and 'Plot DPI' keys

    Returns
    -------
    """
    ndims = len(samples[0])
    lbls = [r"$a_{{{0:d}}}$".format(i) for i in range(ndims)]
    ranges = [(-0.0001, (0.0001 + samples[:, i].max())) for i in range(ndims)]
    quantile_list = [(0.5 - CONFIG["Confidence Interval"] / 2.0), 0.5,
                     (0.5 + CONFIG["Confidence Interval"] / 2.0)]
    fmt_str = "A{0:d}_prob_en_{1:05.2f}_a{2:02d}.{3:s}"
    for i in range(ndims):
        temp = samples[:, i]
        fig = corner.corner(temp, labels=[lbls[i]], range=[ranges[i]],
                            quantiles=quantile_list, truths=[peak_vals[i]],
                            verbose=False, bins=CONFIG["Num Bins"])
        # make the probability plot file_name
        base_file_name = fmt_str.format(CONFIG["Target A"], energy, i,
                                        CONFIG["Plot Format"])
        fig_file_name = os.path.join(CONFIG["Prob Plots Directory"],
                                     base_file_name)
        fig.set_size_inches(CONFIG["Plot Height"], CONFIG["Plot Width"])
        fig.savefig(fig_file_name, bbox_inches='tight', dpi=CONFIG["Plot DPI"])
        plt.close(fig)


def ln_post_prob(params, cs_lib, struct, bounds):
    """This function calculates the log of the post probability function

    Parameters
    ----------
    params : numpy array
        The parameter set to be tested

    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines

    struct : ctypes.void_ptr
        pointer to struct needed by calculation functions

    bounds : array of floats
        the bounds (inclusive) for each parameter

    Global Parameters
    -----------------

    Returns
    -------
    output : float
        the log likelihood of the point, if a parameter is outside the bounds
        for the parameter set, then the likelihood is zero so the log
        likelihood is negative infinity technically we are supposed to return
        the ln of the posterior probability which is the likelihood times the
        prior probability function. However, our prior probability is flat, set
        to 1 for the valid range, so it becomes just the log likelihood plus
        zero
    """
    # first check if we are outside the resonable parameter range, if so,
    # return negative infinity, which corresponds to a probability of 0
    for i, bnd in enumerate(bounds):
        if params[i] < bnd[0] or params[i] > bnd[1]:
            return -np.inf
    return cs_lib.calculateLnLiklihood(
        struct, params.ctypes.data_as(ct.POINTER(ct.c_double)))


def gen_walker_starts(gens):
    """This function takes the set of cluster centers and generates the walker
    start positions from them

    Parameters
    ----------
    gen : list of lists
        This is a list of parameter sets generated by the initial fit process
        it is sorted in order of best fit chi^2 to worst

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Maximum L', 'Number Walker Generators', and
        'Number of Walkers' keys

    Returns
    -------
    output : list of numpy arrays
        There are CONFIG["Number of Walkers"] elements in this list, each is a
        starting point for a walker
    """
    # store the number of dimensions
    ndims = (1 + CONFIG["Maximum L"])
    # generate the output
    output = [randomize_position(gens[i % CONFIG["Number Walker Generators"]],
                                 ndims)
              for i in range(CONFIG["Number of Walkers"])]
    return output


def randomize_position(gen, ndims):
    """This function takes a generator position and suitably randomizes it

    Parameters
    ----------
    gen : numpy array
        initial list of parameters to be randomized

    ndims : int
        number of parameters to be randomized

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Sample Spread', 'Sample Offset Centroid', and
        'Sample Offset Width' keys

    FLOAT_EPSILON : float
        A very small value that is the threshold for "two float are the same"

    Returns
    -------
    position : numpy array
        A randomized starting point for an ensemble walker
    """
    # make an array of randomizers, normally distributed around zero
    rands = CONFIG["Sample Spread"] * np.random.standard_normal(ndims)
    # convert them to fractions of the original value
    randomizer = (np.ones((ndims), dtype=np.float64) - rands)
    # generate the randomized offsets
    base_offsets = (CONFIG["Sample Offset Width"] *
                    np.random.standard_normal(ndims))
    # move the offsets to the appropriate centroid
    offsets = ((CONFIG["Sample Offset Centroid"] *
                np.ones((ndims), dtype=np.float64)) +
               base_offsets)
    # make the randomized positions
    position = (randomizer*(gen + offsets))
    # now if any of the positions are negative, negate them, if they are too
    # close to zero, set them to 0.001
    for i in range(ndims):
        if position[i] < 0.0:
            position[i] *= -1.0
        if position[i] < FLOAT_EPSILON:
            position[i] = 0.001
        if position[i] > (1.0/CONFIG["EWSR Fractions"][i]):
            position[i] = ((1.0/CONFIG["EWSR Fractions"][i]) - 0.001)
    return position


def do_init_fit(start, struct, cs_lib):
    """Performs a fit from the given starting point

    Parameters
    ----------
    start : numpy array
        numpy array containing the start parameters set [a0 to aN] where N is
        the Max L

    struct : ctypes.void_ptr
        pointer to struct needed by calculation functions

    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'EWSR Fractions', and 'Forced Extra Fits' keys

    Returns
    -------
    min_chi^2 : float
        The chi^2 at the best fit point

    best_fit_params : numpy array
        The set of parameters for the best fit
    """
    # calculate the bounds
    bnds = [(0.0, 1.0/x) for x in CONFIG["EWSR Fractions"]]
    # print call_chi_sq(np.array(start, dtype=np.float64), cs_lib, struct)
    init_params = np.array(start, dtype=np.float64)
    ret_vals = optimize.fmin_l_bfgs_b(call_chi_sq,
                                      init_params, bounds=bnds,
                                      epsilon=1e-05, approx_grad=True,
                                      args=(cs_lib, struct), iprint=0,
                                      factr=10.0)
    for _ in range(CONFIG["Forced Extra Fits"]):
        ret_vals = optimize.fmin_l_bfgs_b(call_chi_sq,
                                          ret_vals[0], bounds=bnds,
                                          epsilon=1e-05, approx_grad=True,
                                          args=(cs_lib, struct), iprint=0,
                                          factr=10.0)
    return (ret_vals[1], ret_vals[0])


def call_chi_sq(params, cs_lib, struct):
    """calls the chi^2 function in cs_libParameters

    Parameters
    ----------
    params : numpy array
        numpy array containing the parameters [a0 to aN] where N is the Max L
        to be fit

    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines

    struct : ctypes.void_ptr
        pointer to struct needed by calculation functions

    Global Parameters
    -----------------

    Returns
    -------
    chi^2 : float
        The chi^2 of the data given the data and distributions in the struct
        and the parameters passed to calculateChi
    """
    temp = cs_lib.calculateChi(struct,
                               params.ctypes.data_as(ct.POINTER(ct.c_double)))
    return temp


def make_calc_struct(cs_lib, data, dists):
    """This function takes the data and distributions and dumps the information
    into a freshly created struct

    Parameters
    ----------
    cs_lib : ctypes.cdll.LoadLibrary object
        This is the object representing the loaded dll containing the fast c
        routines

    data : numpy array
        this is the array of experimental data, scaled by errors, that the dll
        will calculate chi^2 values with

    dists : list of numpy arrays
        this is the array of interpolated distributions that have been divided
        by the experimental errors

    Global Parameters
    -----------------

    Returns
    -------
    out_struct : ctypes.void_ptr
        While it is listed as a void pointer this is a pointer to the data
        needed by the dll library routines to perform the log likelihood
        calculations
    """
    # make a struct
    out_struct = cs_lib.makeMdaStruct(len(data), len(dists))
    # load it with the data
    cs_lib.setMdaData(out_struct, data.ctypes.data_as(ct.POINTER(ct.c_double)))
    # iterate through this distributions
    for i, dis in enumerate(dists):
        # load the distributions
        cs_lib.setMdaDist(out_struct, i,
                          dis.ctypes.data_as(ct.POINTER(ct.c_double)))
    # return the struct
    return out_struct


def generate_output_dirs():
    """This function checks for the existence of the directories output is to
    be placed in, if they do not exist, they are created

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Corner Plots Directory', 'Prob Plots Directory',
        'Chain Directory', 'Generate Walker Plots', and
        'Parameter Files Directory' keys

    Returns
    -------
    """
    # test / create the directory for the output file
    if not os.path.exists(CONFIG["Parameter Files Directory"]):
        os.makedirs(CONFIG["Parameter Files Directory"])
    # test / create the directories for csv files with individial fits
    if CONFIG["Generate Fit CSVs"]:
        make_config_fit_csv_dirs()
    # test / create the directories for fit plots
    if CONFIG["Generate Fit Plots"]:
        make_config_fit_plot_dirs()
    # test / create the directory for corner plots
    if CONFIG["Generate Corner Plots"]:
        if not os.path.exists(CONFIG["Corner Plots Directory"]):
            os.makedirs(CONFIG["Corner Plots Directory"])
    # test / create the directory for probability plots
    if CONFIG["Generate Probability Plots"]:
        if not os.path.exists(CONFIG["Prob Plots Directory"]):
            os.makedirs(CONFIG["Prob Plots Directory"])
    # test / create the directory for Markov Chains
    if CONFIG["Save Chain Data"]:
        if not os.path.exists(CONFIG["Chain Directory"]):
            os.makedirs(CONFIG["Chain Directory"])
    # test / create the directory for Parameter Plots
    if CONFIG["Generate Parameter Plots"]:
        make_config_param_plot_dirs()
    # test / create the main directory for time series plots
    if CONFIG["Generate Walker Plots"]:
        if not os.path.exists(CONFIG["Time Series Directory"]):
            os.makedirs(CONFIG["Time Series Directory"])


def make_config_param_plot_dirs():
    """This function makes the parameter plot directories

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Parameter Plots Directory' key and creates (and
        uses) the "Param Plot Dirs' key

    Returns
    -------
    """
    temp = CONFIG["Parameter Plots Directory"]
    CONFIG["Param Plot Dirs"] = [os.path.join(temp, "Percentiles"),
                                 os.path.join(temp, "Peaks")]
    for subdir in CONFIG["Param Plot Dirs"]:
        if not os.path.exists(subdir):
            os.makedirs(subdir)


def make_config_fit_plot_dirs():
    """This function encapsulates the annoying task of making all the sub dir
    names for holding fit plots

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Fit Plots Directory' key and creates (and uses)
        the "Fit Plot Dirs' key

    Returns
    -------
    """
    temp = CONFIG["Fit Plots Directory"]
    subdir1 = ["Percentiles", "Peaks"]
    subdir2 = ["Complete", "Limited"]
    subdir3 = ["Low_Edge", "Parameters", "High_Edge"]
    CONFIG["Fit Plot Dirs"] = []
    # make all combinations of subdir1, subdir2, and subdir3 elements
    for dir1 in subdir1:
        for dir2 in subdir2:
            for dir3 in subdir3:
                CONFIG["Fit Plot Dirs"].append(os.path.join(temp, dir1,
                                                            dir2, dir3))
    for fit_dir in CONFIG["Fit Plot Dirs"]:
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)


def make_config_fit_csv_dirs():
    """This function takes the configuration information and generates the two
    folders that will hold the fit csv files

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Fits Csv Directory' key and creates (and uses)
        the "Fit Csv Dirs' key

    Returns
    -------
    """
    CONFIG["Fit Csv Dirs"] = [os.path.join(CONFIG["Fits Csv Directory"], x) for
                              x in ["percentiles", "peaks"]]
    for fit_dir in CONFIG["Fit Csv Dirs"]:
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)


def calc_start_params():
    """This function uses the config data to generate a list of starting
    points for the initial fits performed before the MCMC

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Start Pts a%d' key where %d is some integer from 0
        to Maximum L, it also uses the 'Maximum L' key

    Returns
    -------
    start_list: list of lists
        Here each sub list is a set of starting values for all the parameters
        in the fit
    """
    # first load the sets into an array
    sl_list = []
    for i in range(CONFIG["Maximum L"]+1):
        sl_list.append(CONFIG["Start Pts a{0:d}".format(i)])
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
        for i, ind in enumerate(indices):
            curr_start.append(sl_list[i][ind])
        start_list.append(curr_start)
        # now increment the indices
        increment_ind(indices, lengths)
    # now return the list of starting points
    return start_list


def increment_ind(ind, lens):
    """This function takes a list of indices and lengths and increments the
    highest indice, resetting to zero as needed

    Parameters
    ----------
    ind : list of ints
        this is a list of indices that needs the lowest one to be incremented
        (with propogation of that increment through the higher order indices
        if the increment pushes the index past a boundary)

    lens : list of ints
        the list of lengths of each dimension of the arrays

    Global Parameters
    -----------------

    Returns
    -------
    """
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
    error for each point

    Parameters
    ----------
    dists: list of lists of numpy arrays
        each sub list contains each of the dwba distributions for a given
        excitation energy with the index of the numpy array in the sub list
        corresponding to its L-Value. The index of the sub lists corresponds to
        the excitation energy of the data at the same index

    data: list of lists
        the first element of the list is the excitation energy of that dataset,
        the second element is a numpy containing the dataset in the format
        (angle, cs, cs-err)

    Global Parameters
    -----------------

    Returns
    -------
    output: list of lists of numpy arrays
        similar to dists in ordering and indexing, however the numpy arrays now
        only contain the values of the input dwba distribution interpolated
        to the angles at which there is data and divided by the error in the
        data
    """
    # first make the output variable
    output = []
    # now iterate across the energies
    for i, dis_set in enumerate(dists):
        # extract the list of angles for this energy
        angle_list = data[i][1][:, 0]
        # extract the list of errorss for this energy
        error_list = data[i][1][:, 2]
        # make the variable to hold the list of interpolated values
        en_output = []
        # iterate across the L values
        for dist in dis_set:
            # interpolate and divide the distribution
            interp_data = interpolate_dist(dist, angle_list,
                                           error_list).astype(np.float64)
            en_output.append(interp_data)
        # append the information for the distributions of this energy to output
        output.append(en_output)
    return output


def interpolate_dist(dist, angles, errors):
    """this function takes a distribution, a list of exp angles, and a list of
    exp error at each angle, it then interpolates the distribution at that
    angle and divides that value by the exp error at that angle

    Parameters
    ----------
    dist: numpy array
        Nx2 array of the read in giant resonance DWBA distribution with the
        format (angle, cs)

    angles: numpy array or list
        Array of angles at which the distributions need to be interpolated

    errors: numpy array or list
        Array of errors for the data at the angles given

    Global Parameters
    -----------------

    Returns
    -------
    scaled_interpolated_dist: numpy array
        Array of the distribution values at the provided angles, scaled by the
        inverse of the errors at those points
    """
    # construct the interpolation
    interp = interpolate.interp1d(dist[:, 0], dist[:, 1], kind="cubic")
    # get the interpolated values divided by the errors and return them
    return (interp(angles)/errors).astype(np.float64)


def read_dist(energy, l_value):
    """This function reads the distribution described by the passed parameters
    from disk into an np array

    Parameters
    ----------
    energy: float
        excitation energy of the giant resonance dwba anglar distribution to be
        read in by this function

    l_value: int
        orbital angular momentum of the giant resonance dwba angular
        distribution to be read in by this function

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Distribution Directory', 'Target A', and
        'EWSR Fractions' keys in the dictionary

    Global Parameters
    -----------------

    Returns
    -------
    gr_dwba_dist: numpy array
        numpy array with shape Nx2 where N is the number of data points in the
        dwba angular distribution which is in the format (angle, dwba_cs)
    """
    # first construct the file name
    fmt_str = "A{0:d}_Ex{1:4.2f}_L{2:02d}_T0_F{3:03d}.csv"
    ewsr_frac = int(100.0*CONFIG["EWSR Fractions"][l_value])
    file_name = fmt_str.format(CONFIG["Target A"], energy, l_value, ewsr_frac)
    dist_file_name = os.path.join(CONFIG["Distribution Directory"], file_name)
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
    times the ewsr fraction

    Parameters
    ----------
    data : list of lists
        the first element of the list is the excitation energy of that dataset,
        the second element is a numpy containing the dataset in the format
        (angle, cs, cs-err)

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Subtract IVGDR', 'IVGDR Height', 'IVGDR Center',
        'IVGDR Width', and 'IVGDR CS Integral' keys in the dictionary

    Returns
    -------
    dists: list of numpy arrays
        these arrays containing the IVGDR distributions at each excitation
        energy that there is data for. If CONFIG["Subtract IVGDR"] is false,
        this is set to None

    ewsrs: list of floats
        each float is the fraction of strength of the IVGDR at the
        corresponding excitation energy in data list. If
        CONFIG["Subtract IVGDR"] is false, this is set to None

    sub_data: list of lists
        same format as the input parameter data, except that each cross-section
        has been reduced by the inerpolated and ewsr-scaled IVGDR cross-section
        at that angle and excitation energy. If CONFIG["Subtract IVGDR"] is
        false, this is merely a deep copy of the input data
    """
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
        for i, sub_set in enumerate(sub_data):
            sub_set[1][:, 1] = sub_set[1][:, 1] - ivgdr_values[i]
        print "IVGDR subtraction is done"
        return dists, ewsrs, sub_data
    else:
        print "IVGDR is turned off"
        return None, None, copy.deepcopy(data)


def interpolate_and_scale_ivgdr(dist, ewsr, angle_list):
    """This function takes an ivgdr distribution, interpolates it, calculates
    the values of the distribution at the provided angles, multiplies them by
    ewsr and returns it

    Parameters
    ----------
    dist : numpy array (Nx2)
        array containing the IVGDR angular distribution in the format angle,
        cross-section

    ewsr : float
        fraction of the IVGDR strength present at that energy

    angle_list : list of floats
        list of angles at which there are data points for this energy

    Global Parameters
    -----------------

    Returns
    -------
    scaled_cross-sections: numpy array of floats
        A list of the calculated IVGDR cross-sections for the given angles
        scaled to the appropriate EWSR fraction
    """
    interp = interpolate.interp1d(dist[:, 0], dist[:, 1], kind="cubic")
    values = interp(angle_list)
    return ewsr*values


def read_ivgdr_dists(en_list):
    """reads in the csv files with the IVGDR distributions

    Parameters
    ----------
    en_list : list
        list of excitation energies (in MeV) whose IVGDR distributions need
        to be read in

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Distribution Directory' and 'Target A' keys in the
        dictionary

    Returns
    -------
    dists_list : list
        A list of distributions numpy arrays with dimensions Nx2 where N is the
        number of points in a distribution. The order of the arrays is
        identical to the ordering of energies in the passed en_list
    """
    fmt_str = "A{0:d}_Ex{1:4.2f}_L01_T1_F100.csv"
    file_names = [fmt_str.format(CONFIG["Target A"], energy) 
                  for energy in en_list]
    dist_file_names = [os.path.join(CONFIG["Distribution Directory"], fname)
                       for fname in file_names]
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
    by its integral

    Parameters
    ----------
    max_sigma: float
        The height (in millibarns) of the lorentzian

    centroid: float
        The centroid energy (in MeV) of the Lorentzian

    width: float
        The width (in MeV) of the Lorentzian

    total_sigma: float
        The integral from 0 to infinity of the lorentizian in millibarns

    Global Parameters
    -----------------

    Returns
    -------
    """
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
        """returns the IVGDR EWSR % at energy = excitation_energy

        Parameters
        ----------
        excitation_energy: float
            The energy to calculate the value of the IVGDR fraction at

        Global Parameters
        -----------------

        Returns
        -------
        ivgdr_percentage: float
            The fraction of the ivgdr strength present at that energy
        """
        en_sq = math.pow(excitation_energy, 2.0)
        denom = self.denom_const + self.denom_inv/en_sq + en_sq/self.width_sq
        return self.sig_ratio/denom

    def get_ivgdr_cs(self, excitation_energy):
        """returns the IVGDR total cs at energy = excitation_energy

        Parameters
        ----------
        excitation_energy: float
            The energy to calculate the value of the IVGDR fraction at

        Global Parameters
        -----------------

        Returns
        -------
        ivgdr_cs: float
            The cross-section at that energy (d_sigma/d_Energy)
        """
        en_sq = math.pow(excitation_energy, 2.0)
        denom = self.denom_const + self.denom_inv/en_sq + en_sq/self.width_sq
        return (self.total * self.sig_ratio) / denom


def read_row_cs_data_file():
    """takes a cross-section file in row format and extracts the data,
    points with angle above max_angle are excluded

    Parameters
    ----------

    Global Parameters
    -----------------
    CONFIG : dictionary
        This uses the CONFIG global dictionary that was read in at program
        start. It uses the 'Input File Path', 'Final Energy', 'Start Energy',
        and 'Max Theta' keys in the dictionary

    Returns
    -------
    fit_output: list of lists
        each sub list contains the excitation energy of the data as the first
        element. The second element is a numpy array with dimensions Nx3, where
        N is the number of data points. Each contains the angle, the cross-
        section and the statistical error in the cross-section, the points in
        this array are limited to only those with angle less than
        CONFIG["Max Theta"]

    plot_output: list of lists
        each sub list contains the excitation energy of the data as the first
        element. The second element is a numpy array with dimensions Nx3, where
        N is the number of data points. Each contains the angle, the cross-
        section and the statistical error in the cross-section, all of the data
        read in is contained in this array, not just the data with small enough
        theta
    """
    # open the csv file
    input_file = open(CONFIG["Input File Path"], 'r')
    fit_output = []
    plot_output = []
    # read the file line by line, each line has an energy and a list of angles
    # and cross-sections and errors
    for line in input_file:
        # break the line by the commas
        data_list = line.split(',')
        # get the energy for the distribution
        energy = float(data_list[0].strip())
        # check if the energy is in the acceptable range
        if energy < CONFIG["Final Energy"] and energy > CONFIG["Start Energy"]:
            # convert each cell to a float
            distribution_data = [float(x.strip()) for x in data_list[1:]]
            # distibution_data = map(lambda x: float(x.strip()), data_list[1:])
            dist_len = len(distribution_data)
            fit_distribution = []
            plot_distribution = []
            i = 0
            # for each trio of elements in the distribution put them in a tuple
            # and put that tuple in the distribution list this makes the
            # distribution variable contain a list of points and errors:
            # [(x1,y1,dy1),(x2,y2,dy2),...]
            for i in xrange(0, dist_len, 3):
                if distribution_data[i] <= CONFIG["Max Theta"]:
                    fit_distribution.append((distribution_data[i],
                                             distribution_data[i+1],
                                             distribution_data[i+2]))
                plot_distribution.append((distribution_data[i],
                                          distribution_data[i+1],
                                          distribution_data[i+2]))
            # put the energy and its associated distribution in a list
            fit_output.append([energy, np.array(fit_distribution,
                                                dtype=np.float64)])
            plot_output.append([energy, np.array(plot_distribution,
                                                 dtype=np.float64)])
    print "Exp Data is read in"
    return (fit_output, plot_output)


SAMPLES_ERROR = """
WARNING: the number of samples:
{0:d} = (({1:d} - {2:d}) * {3:d})
is not large enough to ensure 10 points outside of each error bar location
for the given confidence interval. The minimum number of samples necessary is:
{4:d} = ceiling[ 10 / ((1.0 - {5:12.10f}) / 2.0))

For more information look at the mda_config.py file.
"""


MDA_START_MSG = "Starting MDA, working on up to {0:d} energies simultaneously"


TOO_MANY_EWSR_ERROR = """
Too many EWSR fractions listed, there must be
            (1 + CONFIGURATION[\"Maximum L\"]) = {0:d}
EWSR fractions listed in the variable CONFIGURATION[\"EWSR Fractions\"],
            {1:d} were given"""


TOO_FEW_EWSR_ERROR = """
Too few EWSR fractions listed, there must be
            (1 + CONFIGURATION[\"Maximum L\"]) = {0:d}
EWSR fractions listed in the variable CONFIGURATION[\"EWSR Fractions\"],
            {1:d} were given"""


TOO_MANY_BINS_ERROR = """
Too many corner bins listed, there must be
            (1 + CONFIGURATION[\"Maximum L\"]) = {0:d}
corner bins listed in the variable CONFIGURATION[\"Corner Plot Bins\"],
            {1:d} were given"""


TOO_FEW_BINS_ERROR = """
Too few corner bins listed, there must be
            (1 + CONFIGURATION[\"Maximum L\"]) = {0:d}
corner bins listed in the variable CONFIGURATION[\"Corner Plot Bins\"],
            {1:d} were given"""


NUM_STARTS_ERROR = """
You must provide enough starting points such that there are at
least {0:d} points (the number of start points is the length of each start list
multiplied together)"""


STARTUP_MSG = """
The best way to abort this program while it is using the multiprocessing module
(which it does even when 'Number of Threads' is set to 1) is not to use 
'Ctrl+C' but instead 'Ctrl+Z' followed by 'kill %Job Number' (Job Number is
usually 1) This avoids any insanity from the interruption of the
multiproccessing module
"""

if __name__ == "__main__":
    main()
