'''
File: SameNe_DistortionsPlot.py
Author: Micaila Marcelle
Purpose: Constructs a plot showing the SFS curves corresponding to two primary
cases: one for Ud = 10, Sd = 0.006 and one for Ud = 2, Sd = 0.09. These cases
produce approximately the same reduction in Ne, so this plot allows us to better
compare the actual distortions associated with these Ne reductions. Used to construct
Fig. 3 in paper
'''

# Handles the necessary imports
import math
import io
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Reads in the necessary SFS data from files into arrays so that this data can be plotted
SFS_Ud_2_Sd_09 = []
SFS_Ud_10_Sd_006 = []
census_N_neutral_curve = []
reduced_Ne_neutral_curve = []
X = []


with open("cleanedSFS_output_2_09.txt", "r") as file:
	for line in file:
		line = line.strip()
		SFS_Ud_2_Sd_09.append(float(line))
		
with open("cleanedSFS_output_10_006.txt", "r") as file:
	for line in file:
		line = line.strip()
		SFS_Ud_10_Sd_006.append(float(line))
		
with open("neutralN_10000_output.txt", "r") as file:
	for line in file:
		line = line.strip()
		census_N_neutral_curve.append(float(line))
		
with open("reducedNe_5153_output.txt", "r") as file:
	for line in file:
		line = line.strip()
		reduced_Ne_neutral_curve.append(float(line))
		
with open("xAxis_output.txt", "r") as file:
	for line in file:
		line = line.strip()
		X.append(float(line))
		
# Sets appropriate axis labels for the overall plot
plt.xlabel("Allele count / 2N", fontweight = "bold", fontsize = 20)
plt.ylabel("Relative number of sites", fontweight = "bold", fontsize = 20)

	
# Plots the resulting neutral curves
# #377eb8
# #4daf4a
plt.plot(X, census_N_neutral_curve, label = "Neutral N = 10,000", color = "#A0A0A0", linewidth = 3)
plt.plot(X, reduced_Ne_neutral_curve, label = "Neutral Ne = " + str(5153), color = "#ffb366", linewidth = 3)

# Plots all three SFS
# #e41a1c
# #984ea3
plt.plot(X, SFS_Ud_2_Sd_09, label = "Ud = 2, Sd = 0.09, Ne = 5153", color = "#00B3E5", linewidth = 6, linestyle = "dashed", dashes = (3, 3))
plt.plot(X, SFS_Ud_10_Sd_006, label = "Ud = 10, Sd = 0.006, Ne = 5204", color = "#d31f11", linewidth = 6, linestyle = "dashed", dashes = (3, 3))


# Sets the scaling, legend, and ensures appropriate structure of the plot
plt.xscale("logit")
plt.yscale("log")
plt.ylim(5 * (10 ** -6), 2.0)

# Constructs the legend and sets up the details of the plot layout
plt.legend()
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.tight_layout()

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(13)
plt.savefig("DistortionsComparisonPlot.png")

