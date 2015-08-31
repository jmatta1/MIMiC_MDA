config={}

# path to the file containing the experimental data
config["Input File Path"] = ./test.csv
# path to the shared library that contains the compiled routines for fast calc
# of chi^2, log likelihood, and prior distributions
config["Shared Lib Path"] = "./C/libChiSq.so"

# Sets if IVGDR subtraction is carried carried out.
config["Subtract IVGDR"] = True
# Integral across all ex energy of the IVGDR lorientzian, doesn't matter if
# subtraction is off
config["IVGDR CS Integral"] = 1200.0
# Height of the IVGDR lorentzian in millibarns
config["IVGDR Height"] = 200.0
# centroid energy (in MeV) of the IVGDR lorentzian
config["IVGDR Center"] = 12.0
# width of the IVGDR lorentzian in MeV
config["IVGDR Width"] = 3.0
