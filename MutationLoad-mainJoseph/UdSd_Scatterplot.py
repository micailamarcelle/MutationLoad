'''
File: UdSd_Scatterplot.py
Author: Micaila Marcelle
Purpose: Constructs several scatterplots, used to compare variance in log
fitness, U * s, and Ne/N. Though not all of the scatterplots constructed by
this file were ultimately used in the paper, two were used to construct Fig. 1
'''

# Handles the necessary imports
import math
import io
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Sets up basic parameters
N = 1000

# Sets up our x-axis, representing values of Ud * Sd
# x_ud = [0.0001, 0.0002, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1]
# x_ud_sd = [0.0005, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.025, 0.1, 0.145, 0.25, 0.5]
x_ud_sd = [0.001, 0.004, 0.005, 0.008, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.29]

# We then compute values for the Charlesworth 2012 Ne/N, which is e^-8Ush, with h being taken as 0.5
charlesworth_ne_over_n = []
for value in x_ud_sd:
	charlesworth_ne_over_n.append(math.exp(-8 * value * 0.5))
	
# From here, we want to compare this curve to the results of our simulations, specifically looking at
# a constant Ud of 10
# simulation_ne_over_n_ud_5 = [1.000, 0.983, 0.939, 0.872, 0.816, 0.784, 0.754, 0.652, 0.19, 0.326, 0.135]
simulation_ne_over_n_ud_10 = [1, np.nan, 0.966, np.nan, 0.902, 0.788, 0.712, 0.65, 0.592, 0.47, 0.252]

# We'll also compare this to the result of holding sd constant at 0.004
simulation_ne_over_n_sd_004 = [np.nan, 0.868, np.nan, 0.826, np.nan, 0.784, np.nan, 0.65, np.nan, np.nan, np.nan]

# Sets appropriate axis labels for the overall plot
plt.xlabel("Ud * s", fontweight = "bold", fontsize = 30)
plt.ylabel("Ne/N", fontweight = "bold", fontsize = 30)

# Plots the curve for Charlesworth and the points for our simulations
plt.plot(x_ud_sd, charlesworth_ne_over_n, label = "Charlesworth 2012", color = "#007191", linewidth = 3)
plt.plot(x_ud_sd, simulation_ne_over_n_ud_10, 'o', label = "Observed Ud = 10", color = "#d31f11", markersize = 8)
plt.plot(x_ud_sd, simulation_ne_over_n_sd_004, 's', label = "Observed s = 0.004", color = "#0ea503", markersize = 7)

# Sets the scaling, limits, and legend
plt.xscale("log")
plt.ylim(0, 1.05)
plt.legend()

# Sets font sizes
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(10)
plt.savefig("Charlesworth_vs_simulations_ne.png")



# Constructs a new figure, which will be used to compare various constant values of
# Ud to the Charlesworth 2012 curve
plt.figure()

# Sets up the x-axis (Ud * sd values) and plotted points for Ud = 1
x_ud_1 = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
simulation_ne_over_n_const_ud_1 = [1.00, 0.999, 0.99, 0.947, 0.927, 0.81, 0.704, 0.412]

# Also, for points with sd < 1/Ne, we plot a larger, black circle to create a border
x_ud_1_borders = [0.0001, 0.0005]
ud_1_border_circles = [1.00, 0.999]

# Sets up the x-axis and plotted points for Ud = 5
x_ud_5 = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]
simulation_ne_over_n_const_ud_5 = [1.00, 1.00, 1.00, 0.94, 0.876, 0.619, 0.546, 0.192, 0.087]

# Again, for points with sd < 1/Ne, plot a larger, black circle to create a border
x_ud_5_borders = [0.0001, 0.0005, 0.001, 0.005]
ud_5_border_circles = [1.00, 1.00, 1.00, 0.94]

# Sets up the x-axis and plotted points for Ud = 10
x_ud_10 = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 1]
simulation_ne_over_n_const_ud_10 = [1, 1, 1, 0.965, 0.9, 0.617, 0.489, 0.066]

# Again, for points with sd < 1/Ne, plot a larger, black circle to create a border
x_ud_10_borders = [0.0001, 0.0005, 0.001, 0.005, 0.01]
ud_10_border_circles = [1, 1, 1, 0.965, 0.9]

# Sets up the x-axis and plotted values for Charlesworth
x_charlesworth = [0.0001, 0.0002, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1]
charlesworth_ne_over_n_ud_plot = []
for value in x_charlesworth:
	charlesworth_ne_over_n_ud_plot.append(math.exp(-8 * value * 0.5))
	# charlesworth_ne_over_n_ud_plot.append(1.0 / (1 + 4 * value))
	
# For the sake of considering whether our hypothesis of the combined effects of 
# linked and unlinked BGS holds, we also plot curves representing the Matheson &
# Masel combined equation for each Ud
matheson_masel_ud_1 = []
matheson_masel_ud_5 = []
matheson_masel_ud_10 = []
for value in x_charlesworth:
	sd_for_ud_1 = value / 1
	sd_for_ud_5 = value / 5
	sd_for_ud_10 = value / 10
	matheson_masel_ud_1.append(math.exp(-4 * (1 - (1 / 23)) * sd_for_ud_1 * 0.5) * math.exp(- (1 / 23) / 1))
	matheson_masel_ud_5.append(math.exp(-4 * (5 - (5 / 23)) * sd_for_ud_5 * 0.5) * math.exp(- (5 / 23) / 1))
	matheson_masel_ud_10.append(math.exp(-4 * (10 - (10 / 23)) * sd_for_ud_10 * 0.5) * math.exp(- (10 / 23) / 1))
	
# Sets appropriate axis labels for the overall plot
plt.xlabel("Ud * sh", fontweight = "bold", fontsize = 30)
plt.ylabel("Ne/N", fontweight = "bold", fontsize = 30)

# Plots the Charlesworth curve
plt.plot(x_charlesworth, charlesworth_ne_over_n_ud_plot, label = "Charlesworth 2012", color = "#007191", linewidth = 3)

# Plots the Matheson & Masel curves
plt.plot(x_charlesworth, matheson_masel_ud_1, label = "Matheson & Masel Ud = 1", color = "#fff132", linewidth = 1)
plt.plot(x_charlesworth, matheson_masel_ud_5, label = "Matheson & Masel Ud = 5", color = "#fd8d3c", linewidth = 1)
plt.plot(x_charlesworth, matheson_masel_ud_10, label = "Matheson & Masel Ud = 10", color = "#e31a1c", linewidth = 1)

# Plots the curve for the constant Ud simulations
plt.plot(x_ud_10_borders, ud_10_border_circles, 's', color = "black", markersize = 11)
plt.plot(x_ud_10, simulation_ne_over_n_const_ud_10, 's', label = "Observed Ud = 10", color = "#e31a1c", markersize = 9)
plt.plot(x_ud_1_borders, ud_1_border_circles, 'o', color = "black", markersize = 11)
plt.plot(x_ud_1, simulation_ne_over_n_const_ud_1, 'o', label = "Observed Ud = 1", color = "#fff132", markersize = 9)
plt.plot(x_ud_5_borders, ud_5_border_circles, '^', color = "black", markersize = 11)
plt.plot(x_ud_5, simulation_ne_over_n_const_ud_5, '^', label = "Observed Ud = 5", color = "#fd8d3c", markersize = 9)


# Sets the scaling, limits, and legend
plt.xscale("log")
plt.ylim(0, 1.05)

# Sets font sizes
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)

handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 2, 3, 1]
plt.legend([handles[i] for i in order], [labels[i] for i in order], fontsize = 16)
# plt.legend()

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(8.5)
fig.set_figwidth(10)
plt.savefig("Charlesworth_vs_simulations_ne_const_ud.png")



# Constructs the final new figure, which will be used to compare constant
# values of Sd to the Charlesworth curve
plt.figure()

# Sets up the x-axis (again, Ud * sd values) and plotted points for Sd = 0.001
# Desired additional Ud * Sd values: 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1
x_sd_001 = [0.001, 0.002, 0.005, 0.01, 0.015]
simulation_ne_over_n_const_sd_001 = [0.99, 0.97, 0.94, 0.9, 0.859]

# To indicate the points for which sd < 1/Ne, we also plot another mark behind these
# points to create a clear border
x_sd_001_borders = [0.001, 0.002, 0.005, 0.01, 0.015]
sd_001_border_circles = [0.99, 0.97, 0.94, 0.9, 0.859]

# Sets up the x-axis and plotted points for Sd = 0.005
# Desired additional Ud * Sd values: 0.001, 0.075, 0.1, 0.25, 0.5, 0.75, 1
x_sd_005 = [0.001, 0.005, 0.01, 0.025, 0.0375, 0.05, 0.075]
simulation_ne_over_n_const_sd_005 = [0.991, 0.947, 0.896, 0.76, 0.683, 0.617, 0.522]

# Again, uses additional plotted marks to indicate points for which sd < 1/Ne
# In this case, no such points!

# Sets up the x-axis and plotted points for Sd = 0.01
# Desired additional Ud * Sd values: 0.001, 0.005, 0.025, 0.075, 0.25, 0.5, 0.75, 1
x_sd_01 = [0.001, 0.005, 0.01, 0.02, 0.05, 0.075, 0.1]
simulation_ne_over_n_const_sd_01 = [0.992, 0.964, 0.927, 0.856, 0.619, 0.56, 0.489]

# Again, uses additional plotted marks to indicate points for which sd < 1/Ne
# In this case, no such points!

# Sets up the x-axis and plotted points for Sd = 0.05
# Desired additional Ud * Sd values: 0.001, 0.005, 0.01, 0.025, 0.075, 0.75, 1
x_sd_05 = [0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.375, 0.5, 0.75]
simulation_ne_over_n_const_sd_05 = [0.997, 0.979, 0.96, 0.899, 0.81, 0.649, 0.349, 0.227, 0.156, 0.088]

# Sets appropriate axis labels for the overall plot
plt.xlabel(r'$U_{d} \cdot sh$' + '\n', fontweight = "bold", fontsize = 40)
plt.ylabel(r'$N_{e} / N$', fontweight = "bold", fontsize = 40)

# Plots the curve for Charleworth and the constant Ud simulations
plt.plot(x_charlesworth, charlesworth_ne_over_n_ud_plot, label = "Unlinked BGS expectations", color = "#007191", linewidth = 3)
plt.plot(x_sd_001_borders, sd_001_border_circles, 'o', color = "black", markersize = 18)
plt.plot(x_sd_001, simulation_ne_over_n_const_sd_001, 'o', label = "sh = 0.001", color = "#fed976", markersize = 14)
plt.plot(x_sd_005, simulation_ne_over_n_const_sd_005, 's', label = "sh = 0.005", color = "#fd8d3c", markersize = 14)
plt.plot(x_sd_01, simulation_ne_over_n_const_sd_01, '^', label = "sh = 0.01", color = "#fc4e2a", markersize = 14)
plt.plot(x_sd_05, simulation_ne_over_n_const_sd_05, 'v', label = "sh = 0.05", color = "#b10026", markersize = 14)

# Sets the scaling, limits, and legend
plt.xscale("log")
plt.ylim(0, 1.05)
plt.legend(fontsize=23, markerscale=2.1)

# Sets font sizes
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(10)
fig.set_figwidth(12)
plt.savefig("Charlesworth_vs_simulations_ne_const_sd.png")






# To compare using Ud * sd for the variance in log fitness against directly using the variance
# in log fitness output by the simulations, we also have a second version of the above figure,
# in which the x-position of each point is the observed variance in log fitness
plt.figure()

# Sets up the x-axis and plotted points for sd = 0.001, with the x-axis now representing the 
# observed variance in log fitness rather than Ud * sd
x_sd_001_observed = [0.000688461041, 0.00135541954, 0.0032539794, 0.006208274279, 0.008957101505]
simulation_ne_over_n_const_sd_001_observed = [0.989, 0.978, 0.947, 0.904, 0.868]

# To indicate the points for which sd < 1/Ne, we also plot another mark behind these
# points to create a clear border
x_sd_001_borders_observed = [0.000688461041, 0.00135541954, 0.0032539794, 0.006208274279, 0.008957101505]
sd_001_border_circles_observed = [0.989, 0.978, 0.947, 0.904, 0.868]

# Sets up the x-axis and plotted points for sd = 0.005
x_sd_005_observed = [0.001004157315, 0.005040064996, 0.01008160062, 0.02514041267, 0.03765324535, 0.04991276221, 0.07289986029]
simulation_ne_over_n_const_sd_005_observed = [0.99, 0.945, 0.9, 0.771, 0.687, 0.624, 0.526]

# Sets up the x-axis and plotted points for sd = 0.01
x_sd_01_observed = [0.001008991421, 0.005058350642, 0.01010844362, 0.02031052185, 0.0511885836, 0.07722211446, 0.1035505997]
simulation_ne_over_n_const_sd_01_observed = [0.994, 0.964, 0.928, 0.861, 0.688, 0.577, 0.492]

# Sets up the x-axis and plotted points for sd = 0.05
x_sd_05_observed = [0.001049504606, 0.005278162653, 0.01057457997, 0.02651468127, 0.05336638886, 0.1079965781, 0.2806646988, 0.4345654024, 0.5974998559, 0.9500965962]
simulation_ne_over_n_const_sd_05_observed = [0.998, 0.981, 0.961, 0.9, 0.809, 0.657, 0.36, 0.23, 0.154, 0.087]

# Sets appropriate axis labels for the overall plot
plt.xlabel("Variance in Log Fitness" + '\n', fontweight = "bold", fontsize = 30)
plt.ylabel(r'$N_{e} / N$', fontweight = "bold", fontsize = 40)

# Plots the curve for Charleworth and the constant Sd simulations
plt.plot(x_charlesworth, charlesworth_ne_over_n_ud_plot, label = "Unlinked BGS expectations", color = "#007191", linewidth = 3)
plt.plot(x_sd_001_borders_observed, sd_001_border_circles_observed, 'o', color = "black", markersize = 18)
plt.plot(x_sd_001_observed, simulation_ne_over_n_const_sd_001_observed, 'o', label = "s = 0.001", color = "#fed976", markersize = 14)
plt.plot(x_sd_005_observed, simulation_ne_over_n_const_sd_005_observed, 's', label = "s = 0.005", color = "#fd8d3c", markersize = 14)
plt.plot(x_sd_01_observed, simulation_ne_over_n_const_sd_01_observed, '^', label = "s = 0.01", color = "#ea330e", markersize = 14)
plt.plot(x_sd_05_observed, simulation_ne_over_n_const_sd_05_observed, 'v', label = "s = 0.05", color = "#b10026", markersize = 14)

# Sets the scaling, limits, and legend
plt.xscale("log")
plt.ylim(0, 1.05)
plt.legend(fontsize=23, markerscale=2.1)

# Sets font sizes
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(10)
fig.set_figwidth(12)
plt.savefig("const_sd_observed_variance.png")






# Also, to determine how well Ud * s works as a model for variance in log fitness,
# we plot variance in log fitness vs. Ud * s
plt.figure()

# First, sets up the points corresponding to s = 0.001
x_ud_s_001 = x_sd_001
y_variance_001 = x_sd_001_observed

# Next, sets up the points corresponding to s = 0.005
x_ud_s_005 = x_sd_005
y_variance_005 = x_sd_005_observed

# Next, sets up the points corresponding to s = 0.01
x_ud_s_01 = x_sd_01
y_variance_01 = x_sd_01_observed

# Finally, sets up the points corresponding to s = 0.05
x_ud_s_05 = x_sd_05
y_variance_05 = x_sd_05_observed

# Sets up the x-axis, which corresponds to the variance in log fitness for a sequence
# of simulations
# y_variance_in_log_fitness = [0.0006756698375, 0.00101420272, 0.001343005083, 0.003233647986, 0.005060359339, 0.0061657119, 0.007278035732, 0.008858118524, 0.0100797584, 0.01258478135, 0.02037679527, 0.02515136236, 0.03752580984, 0.05124300504, 0.07269014919, 0.08069059612, 0.1036601665, 0.1364256604, 0.2814049671, 0.4354388166, 0.5988994561, 0.9455273607]

# Sets up the y-axis, which corresponds to the Ud * s values associated with these 
# observed variances in log fitness
# x_ud_s = [0.001, 0.001, 0.002, 0.005, 0.005, 0.01, 0.012, 0.015, 0.01, 0.0125, 0.02, 0.025, 0.0375, 0.05, 0.075, 0.075, 0.1, 0.125, 0.25, 0.375, 0.5, 0.75]

# Sets up the labels for the x- and y-axes
# plt.xlabel("Ud * sh", fontweight = "bold", fontsize = 30)
plt.ylabel("Variance in Log Fitness", fontweight = "bold", fontsize = 30)
plt.xlabel(r'$U_{d} \cdot s$' + '\n', fontweight = "bold", fontsize = 40)

# For reference, we first plot a line corresponding to y = x, representing the case in whcih 
# Ud * s perfectly describes the variance in log fitness 
linear_baseline = np.linspace(0.001, 1, 500)
plt.plot(linear_baseline, linear_baseline, label = "Variance in log fitness = Ud * s", linestyle = "dashed", color = "gray")

# Plots the actual points
plt.plot(x_ud_s_001, y_variance_001, 'o', color = "black", markersize = 18)
plt.plot(x_ud_s_001, y_variance_001, 'o', label = "s = 0.001", color = "#fed976", markersize = 14)
plt.plot(x_ud_s_005, y_variance_005, 's', label = "s = 0.005", color = "#fd8d3c", markersize = 14)
plt.plot(x_ud_s_01, y_variance_01, '^', label = "s = 0.01", color = "#ea330e", markersize = 14)
plt.plot(x_ud_s_05, y_variance_05, 'v', label = "s = 0.05", color = "#b10026", markersize = 14)

# Sets the scaling and legend
plt.xscale("log")
plt.yscale("log")
plt.legend(fontsize=20, markerscale=1.8)

# Sets font sizes
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)

# Finally, saves the resulting figure
fig = plt.gcf()
fig.set_figheight(10)
fig.set_figwidth(11)
plt.savefig("VarianceInLogFitness_vs_UdSd.png")