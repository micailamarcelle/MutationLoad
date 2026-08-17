#!/usr/bin/python

'''
File: tskitAnalysis_bottlegrowth.py
Author: Micaila Marcelle
Purpose: Analyzes a given run of the simulation to construct the SFS, then feeds
this SFS into dadi and aims to fit a bottleGrowth model to it. Ultimately not used 
within the paper, since this model was not found to provide a good fit. 
'''

# Imports all necessary libraries/packages
import math
import tskit
import io
import numpy as np
import matplotlib.pyplot as plt
import msprime
import pandas as pd
from scipy.signal import savgol_filter
# import scipy
# import nlopt
import dadi

print("Successfully imported necessary modules")

# Reads in the necessary tskit tables for all generations of interest, adding
# the results to the corresponding list for each type of table
node_list = []
edge_list = []
site_list = []
mut_list = []

# For convenience, we specify the starting generation, ending generation, and
# tskit-writing frequency to simplify the task of changing these values for
# different simulation runs
start_gen = 100000
end_gen = 200000
tskit_freq = 5000

# We iterate through all of the tskit files that we want to include in our
# analysis, reading in the tables corresponding to the generations of interest
# and appending to the appropriate list
'''
i = start_gen
while (i < end_gen):
    # Reads in the data from the tables for the current generation
    with open("nodetable_gen" + str(i) + ".txt") as file:
        node_list.append(file.read())
    with open("edgetable_gen" + str(i) + ".txt") as file:
        edge_list.append(file.read())
    with open("sitetable_gen" + str(i) + ".txt") as file:
        site_list.append(file.read())
    with open("mutationtable_gen" + str(i) + ".txt") as file:
        mut_list.append(file.read())
        
    # Increments i to access the tables for the next generation of interest
    i += tskit_freq
'''

# The names of the files for the final generation, however, do not follow the
# same pattern as those for the other generations, so we handle reading in these
# final files separately

with open("nodetable.txt") as file:
    node_list.append(file.read())
with open("edgetable.txt") as file:
    edge_list.append(file.read())
with open("sitetable.txt") as file:
    site_list.append(file.read())
with open("mutationtable.txt") as file:
    mut_list.append(file.read())

    
print("Successfully opened table files")
    
# Uses these tables in order to construct the corresponding tree sequences
# Note that the process of constructing tree sequences is repeated for each
# checkpoint, leading to a number of tree sequences equal to the number of
# checkpoints
tree_sequence_list = []
for i in range(len(node_list)):
     tree_sequence_list.append(tskit.load_text(io.StringIO(node_list[i]), 
                               io.StringIO(edge_list[i]), 
                               io.StringIO(site_list[i]), 
                               io.StringIO(mut_list[i]), 
                               strict = False))

print("Successfully constructed tree sequences from tables")

# Adds neutral mutations onto the resulting tree sequence
# This is done for each tree sequence within the list generated above, all with 
# the same neutral mutation rate and model
rate = 0.00001
for i in range(len(tree_sequence_list)):
    tree_sequence_list[i] = msprime.sim_mutations(tree_sequence_list[i], rate = rate, model = msprime.InfiniteAlleles(), discrete_genome = False, keep = False)

print("Successfully simulated neutral mutations")

# Generates the unfolded site frequency spectrum for each tree within the list
site_frequency_spectrum_list = []
for i in range(len(tree_sequence_list)):
    site_frequency_spectrum_list.append(tree_sequence_list[i].allele_frequency_spectrum(span_normalise = False, polarised = True))

# Finally, we average these site frequency spectra to form a single time-averaged
# site frequency spectrum for the run currently being considered 
siteFrequencySpectrum = []
for i in range(len(site_frequency_spectrum_list[0])):
     cur_sum = 0
     for j in range(len(site_frequency_spectrum_list)):
          cur_sum += site_frequency_spectrum_list[j][i]
     cur_avg = cur_sum / len(site_frequency_spectrum_list)
     siteFrequencySpectrum.append(cur_avg)
     
nonNormal_SFS = siteFrequencySpectrum
print("Successfully generated site frequency spectrum")
    

# Defines the inverse function f(x) = 1 / x for convenience in plotting
def invFunc(x):
    return 1 / x

# Defines a helper function for the Cvijovic et al. distortions
def g(x, n, ud, sd):
    return 1 / (1 + np.exp(n * sd * x * np.exp(-ud / sd)))

# Defines a piecewise function according to the distortions outlined by Cvijovic et al. (2018)
def distortions(x, n, sigma, ud, sd, un, expon):
      if (x < 1 / (N * sigma)):
            return (2 * n * un) / x
      elif (x > 1 / (n * sigma) and x < expon / (N * sd)):
            return (n * un) / (n * sd * x * x * np.sqrt((ud / sd) * np.emath.logn((ud / sd), (expon / (n * sd * x)))))
      elif (x > expon / (n * sd) and x < 1 - (expon / (n * sd))):
            return (2 * n * un) / (expon * x)
      elif (x > 1 - (expon / (n * sd)) and x < 1 - (1 / (n * sigma))):
            return (n * un) / (n * sd * (1 - x) * np.emath.logn((ud / sd), (expon / (n * sd * (1 - x)))))
      elif (x > 1 - (1 / (n * sigma))):
            return (2 * n * un) / x
      else:
            return np.nan

    
# Sets up any necessary parameters, which are determined via keyboard input for ease in updating
N = 1000 # int(input("Enter N: "))
L = 23000 # int(input("Enter L: "))
Ud_perLinkage = 2.5 / L # float(input("Enter Ud (whole-genome): ")) / L
Sd = 0.004 # float(input("Enter Sd: "))
sigma = np.sqrt(Ud_perLinkage * Sd)
pos_expon = np.exp(Ud_perLinkage / Sd)
# To find theta, we average the values for our tree sequences
theta_sum = 0
for i in range(len(tree_sequence_list)):
     theta_sum += tree_sequence_list[i].diversity(mode = "branch")
theta = theta_sum / len(tree_sequence_list)
cvijovic_ne = N * np.exp(-Ud_perLinkage / Sd) 
joseph_ne = N * np.exp(-8 * ((Ud_perLinkage * L) - Ud_perLinkage) * Sd * 0.5) * np.exp(-1 * (Ud_perLinkage * L / 23) / (2) ) 
coal_ne = theta / (4 * N)

print("Ne: " + str(coal_ne * 2))



# Uses the site frequency spectrum to generate the associated histogram
# Note that this code is heavily derived from http://sesame.uoregon.edu/~adkern/stdpopsim/doc/tutorial.html 

# Sets the number of bars to consider for the histogram (in case only a partial SFS is desired)
bars = len(siteFrequencySpectrum)

# Constructs an appropriate x-axis based on Cvijovic et al.
X = [x / (2 * N) for x in range(1, bars)]

# Sets appropriate axis labels for the overall plot
plt.xlabel("Allele count / 2N", fontweight = "bold", fontsize = 20)
plt.ylabel("Normalized number of sites", fontweight = "bold", fontsize = 20)

# Constructs the arrays for both of the neutral curves, along with for the Cvijovic et al. distortions
inverses_n = []
inverses_cvijovic_ne = []
inverses_joseph_ne = []
inverses_coal_ne = []
distortion_list = []
for i in range(0, len(X)):
    inverses_n.append((2 * N * rate) / X[i])
    inverses_cvijovic_ne.append((2 * cvijovic_ne * rate) / X[i])
    inverses_joseph_ne.append((2 * joseph_ne * rate) / X[i])
    inverses_coal_ne.append((2 * coal_ne * rate) / X[i])
    distortion_list.append(distortions(X[i], N, sigma, Ud_perLinkage, Sd, rate, pos_expon))
    
# Cleans the distortions data, removing any significant spikes and replacing 
# them simply with NaN values 
# Note that the adjustments to the upper and lower bound will likely need to be adjusted
# for different parameter ranges

'''
upper_bound = 1 - (pos_expon / (N * Sd)) + 0.1
lower_bound = pos_expon / (N * Sd) - 0.11
for i in range(0, len(X)):
	if (distortion_list[i] > inverses_n[i]):
		distortion_list[i] = np.nan
	elif (X[i] > lower_bound and X[i] < upper_bound and distortion_list[i] != inverses_cvijovic_ne[i]):
		distortion_list[i] = np.nan
'''


# Normalizing the neutral curves & Cvijovic et al. distortions
norm_val = inverses_n[0]
for i in range(0, len(X)):
    inverses_n[i] /= norm_val
    inverses_cvijovic_ne[i] /= norm_val 
    inverses_joseph_ne[i] /= norm_val
    inverses_coal_ne[i] /= norm_val
    distortion_list[i] /= norm_val
    
print("Normalization value: " + str(norm_val))
    
# Normalizes the site frequency spectrum to better match Cvijovic et al. figures
norm_SFS = []
total_sum = sum(siteFrequencySpectrum)
for i in range(0, len(siteFrequencySpectrum)):
    norm_SFS.append(siteFrequencySpectrum[i] / (norm_val * L))
    
# Uses a moving average to try and "smooth" the actual SFS
clean_SFS = pd.DataFrame({"data" : norm_SFS[1:]})
clean_SFS = clean_SFS.rolling(window = 50).mean()
clean_SFS = clean_SFS['data'].to_list()
i = 0
while (np.isnan(clean_SFS[i])):
     clean_SFS[i] = norm_SFS[i + 1]
     i += 1
     
for j in range(0, 210):
     clean_SFS[i + j] = norm_SFS[i + j]
     
# Alternatively, uses a Savitzky-Golay filter in order to smooth particularly noisy data
golay_SFS = savgol_filter(norm_SFS[1:], 200, 3, mode="nearest")
for i in range(0, len(golay_SFS)):
     if golay_SFS[i] <= 10 ** -9:
          golay_SFS[i] = np.nan
          
# Plots the actual SFS, ensuring agreement with current axis length
while(len(clean_SFS) != len(X)):
     norm_SFS.append(np.nan)
     clean_SFS.append(np.nan)
# plt.plot(X, norm_SFS[1:], "o", label = "Observed", color = "#ffcaca")
          
# Plotting both neutral curves for N and Ne
plt.plot(X, inverses_n, label = "Neutral with N = " + str(N), color = "#f47a00", linewidth = 5)
plt.plot(X, inverses_coal_ne, label = "Neutral with Ne = " + str(round(coal_ne * 2)), color = "#31aebb", linewidth = 5)

# plt.plot(X, inverses_joseph_ne, label = "Whole-genome background selection Ne", color = "#ff8000", linewidth = 5)
# plt.plot(X, inverses_cvijovic_ne, label = "Independent asexual linkage blocks Ne", color = "#007191", linewidth = 5, linestyle = "dashed", dashes = (3, 3))
# plt.plot(X, inverses_coal_ne, label = "Coalescent Ne", color = "#31aebb", linewidth = 5)
    
# Plots the cleaned SFS 
plt.plot(X, clean_SFS, label = "Observed SFS", color = "#d31f11", linewidth = 5, linestyle = "dashed", dashes = (3, 3))
# plt.plot(X, golay_SFS, label = "Savitzky-Golay Smoothed SFS", color = "r", linewidth = 2)

# Writes the cleaned SFS into a file to allow for later re-use
with open("cleanedSFS_output.txt", "w") as file:
    for value in clean_SFS:
        file.write(str(value) + "\n")
        
# Writes the reduced Ne curve into a file to allow for later re-use
with open("reducedNe_" + str(round(coal_ne * 2)) + "_output.txt", "w") as file:
    for value in inverses_coal_ne:
        file.write(str(value) + "\n")
        
# Writes the neutral N curve into a file to allow for later re-use
with open("neutralN_" + str(N) + "_output.txt", "w") as file:
    for value in inverses_n:
        file.write(str(value) + "\n")
        
# Writes the x-axis into a file to allow for later re-use
with open("xAxis_output.txt", "w") as file:
    for value in X:
        file.write(str(value) + "\n")


# Plotting the Cvijovic et al. distortions
# plt.plot(X, distortion_list, label = "Independent linkage block prediction", color = "#007191", linewidth = 5, linestyle = "dashed", dashes = (3, 3))

# Sets the scaling, legend, and ensures appropriate structure of the plot
plt.xscale("logit")
plt.yscale("log")
plt.ylim(10 ** -6, 2.0)

# Adjusts the order in the legend to improve readability
handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 2, 1]
plt.legend([handles[i] for i in order], [labels[i] for i in order], fontsize = 16)
plt.legend()

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.tight_layout()

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(13)
plt.savefig("SFS_plus_expectations.png")


# Now, we utilize the functionality of dadi, applying it to our generated SFS in order
# to determine whether the signal that we're seeing is something that dadi is able to 
# find as well

# To do this, we must first convert our SFS into the format required by dadi
dadi_spectrum = dadi.Spectrum(nonNormal_SFS)

# Saves the dadi format of the site frequency spectrum to file
dadi_spectrum.to_file("dadi_SFS.fs")

# Sets the sample size to be the entire sample (1000 individuals * diploid) for now, and sets the
# grid points approximately based on this sample size
sample_size = 2000
grid_points = [sample_size + 120, sample_size + 130, sample_size + 140]
sample_size = (sample_size,)

# We then set the demographic model to be the bottlegrowth_1d model, which is used for
# instantaneous size change followed by exponential growth
demographic_model = dadi.Demographics1D.bottlegrowth_1d
# demographic_model = dadi.Demographics1D.snm_1d
# demographic_model = dadi.Demographics1D.growth

# Since our generated SFS is unfolded, we wrap the demographic model in another function
# to account for this
demographic_model = dadi.Numerics.make_anc_state_misid_func(demographic_model)

# We also wrap the demographic model in a function that uses grid points for greater accuracy
demographic_model_ex = dadi.Numerics.make_extrap_func(demographic_model)

# The starting parameters are defined in an arbitrary manner
params = [1, 1, 0.01, 0.01]
# params = [1, 0.01, 0.01]

# We then define the boundaries of optimization
lower_bounds = [1e-2, 1e-2, 1e-3, 0]
upper_bounds = [100, 100, 10, 0.999]
# lower_bounds = [1e-2, 1e-3, 0]
# upper_bounds = [100, 100, 0.999]

# Note that the optimization is run 20 times in order to consider how well the parameters are converging
for i in range(20):
    # From here, we perturb our starting parameters
    params = dadi.Misc.perturb_params(params, fold = 1, upper_bound = upper_bounds, lower_bound = lower_bounds)

    # The optimization is then run
    print("Beginning optimization")
    popt, ll = dadi.Inference.opt(params, dadi_spectrum, demographic_model_ex, grid_points, lower_bound = lower_bounds, upper_bound = upper_bounds, verbose = len(params))
    print("Ending optimization")

# With these optimized parameters, we generate the model SFS that best fits these parameters, using
# our bottlegrowth_1d demographic model
generated_sfs = demographic_model_ex(popt, sample_size, grid_points)

# This model SFS is then plot
fig = plt.figure()
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(generated_sfs, dadi_spectrum)
fig.savefig("dadi_model_fit_bottlegrowth.png")

