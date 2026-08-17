'''
File: log_likelihood_significance.py
Author: Micaila Marcelle
Purpose: Initially constructed to try and determine the significance of a 
given log likelihood. Ultimately, though, determined that this wasn't the 
best approach for the paper, so not ultimately used
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
N = 1000 
L = 23000 
Ud_perLinkage = 2.5 / L 
Sd = 0.004 
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
    



# For comparing the log likelihoods, we first get the observed SFS
dadi_spectrum = dadi.Spectrum(nonNormal_SFS)

# Sets the sample size to be the entire sample (1000 individuals * diploid) for now, and sets the
# grid points approximately based on this sample size
sample_size = 2000
grid_points = [sample_size + 120, sample_size + 130, sample_size + 140]
sample_size = (sample_size,)

# Next, we set up the models
demographic_model_growth = dadi.Demographics1D.growth
demographic_model_bottle = dadi.Demographics1D.bottlegrowth_1d
demographic_model_neutral = dadi.Demographics1D.snm_1d

demographic_model_growth = dadi.Numerics.make_anc_state_misid_func(demographic_model_growth)
demographic_model_bottle = dadi.Numerics.make_anc_state_misid_func(demographic_model_bottle)
demographic_model_neutral = dadi.Numerics.make_anc_state_misid_func(demographic_model_neutral)

demographic_model_ex_growth = dadi.Numerics.make_extrap_func(demographic_model_growth)
demographic_model_ex_bottle = dadi.Numerics.make_extrap_func(demographic_model_bottle)
demographic_model_ex_neutral = dadi.Numerics.make_extrap_func(demographic_model_neutral)

# The approx. best-fit parameters in each case, found through several optimizations
# Specifically for Ud = 5, Sd = 0.004
popt_growth = np.array([1.25079, 0.335456, 0])
popt_bottle = np.array([1.02829, 1.25105, 0.311103, 0])
popt_neutral = np.array([])

# Gets the associated log likelihoods
model_growth = demographic_model_ex_growth(popt_growth, sample_size, grid_points)
ll_growth = dadi.Inference.ll_multinom(model_growth, dadi_spectrum)

model_bottle = demographic_model_ex_growth(popt_bottle, sample_size, grid_points)
ll_bottle = dadi.Inference.ll_multinom(model_bottle, dadi_spectrum)

# We then need to consider the neutral model, first computing a log likelihood for it
model_neutral = demographic_model_ex_neutral(popt_neutral, sample_size, grid_points)
ll_neutral = dadi.Inference.ll_multinom(model_neutral, dadi_spectrum)




