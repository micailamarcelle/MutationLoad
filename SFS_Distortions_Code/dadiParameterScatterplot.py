'''
File: dadiParameterScatterplot.py
Author: Micaila Marcelle
Purpose: Constructs two scatterplots. The first expresses that while the growth ratio remains
the same when sample size is varied, the number of generations over which the growth occurs
decreases. The second scatterplot expresses that as Ud gets smaller, the inferred growth
ratio decreases (i.e. dadi infers less of a growth). Used to construct Fig. 5 within paper
'''

# Handles the necessary imports
import math
import io
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Sets up both the x-axis and the points being plotted for plot 1
x_sample_size = [100, 500, 1000, 2000]
y_growth_ratio = [1.52, 1.52, 1.53, 1.54]
y_generations = [745, 587, 517, 450]

# Constructs the corresponding figure, which will be built with two separate y-axes to avoid
# an overload of figures within the paper
fig1, ax1 = plt.subplots()

ax1.set_xlabel("Sample size", fontsize=16, fontweight='bold')
ax1.set_ylabel("Contemporary / ancient population size", color="#e31a1c", fontsize=16, fontweight='bold')
ax1.set_ylim(0, 2)

ax2 = ax1.twinx()

ax2.set_ylabel("Length of exponential growth \n(in generations)", color="#202fd8", fontsize=16, fontweight='bold')
ax2.set_ylim(0, 1000)
ax2.plot(x_sample_size, y_generations, 's', color = "#202fd8", markersize = 9)
ax1.plot(x_sample_size, y_growth_ratio, 'o', color = "#e31a1c", markersize = 9)

ax1.tick_params(axis='both', labelsize=12)
ax2.tick_params(axis='both', labelsize=12)

plt.show()




# Sets up the x- and y-axes for the second figure
x_ud = [2.5, 5, 10]
y_generations = [555, 504, 420]
y_growth_ratio = [1.13, 1.27, 1.54]

# Constructs the corresponding figure, which again is built with two separate y-axes
fig1, ax1 = plt.subplots()

ax1.set_xlabel("Ud", fontsize=16, fontweight='bold')
ax1.set_ylabel("Contemporary / ancient population size", color="#e31a1c", fontsize=16, fontweight='bold')
ax1.set_ylim(0, 2)

ax2 = ax1.twinx()

ax2.set_ylabel("Length of exponential growth \n(in generations)", color="#202fd8", fontsize=16, fontweight='bold')
ax2.set_ylim(0, 1000)
ax2.plot(x_ud, y_generations, 's', color = "#202fd8", markersize = 9)
ax1.plot(x_ud, y_growth_ratio, 'o', color = "#e31a1c", markersize = 9)

ax1.tick_params(axis='both', labelsize=16)
ax2.tick_params(axis='both', labelsize=16)

plt.show()