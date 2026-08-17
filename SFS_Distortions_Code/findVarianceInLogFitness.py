'''
File: findVarianceInLogFitness.py
Author: Micaila Marcelle
Purpose: Finds mean variance in log fitness for a given run of the simulation
using output from the corresponding simulation, and excluding burn-in phase (i.e.,
furst 100,000 generations)
'''

import math
import tskit
import io
import numpy as np
import matplotlib.pyplot as plt
import msprime
import pandas as pd

# Gets the file name from the user
file_name = input("Enter CSV file name: ")

# Reads this CSV file into a data frame
df = pd.read_csv(file_name)

# Gets the mean, excluding any values prior to generation
# 100,000, since all burn-ins end before this point
df = df[df["Nxtimesteps"] > 100000]
# df = df[df["Variance.in.log.fitness"] <= 0.01]
mean_val = df["Variance.in.log.fitness"].mean()

# Prints the resulting mean value
print(mean_val)