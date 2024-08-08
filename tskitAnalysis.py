#!/usr/bin/python

'''
The following program builds on the tskit library in order to analyze the succinct
tree sequences constructed by the simulation and to generate a site frequency spectrum
from these. Note that this code requires numpy, matplotlib, and msprime to be included
in the current working environment.

Author: Micaila Marcelle
'''

# Imports all necessary libraries/packages
import tskit
import io
import numpy as np
import matplotlib.pyplot as plt
import msprime

print("Successfully imported necessary modules")

# Reads in the necessary tskit tables 
with open("nodetable.txt") as file:
    node_table = file.read()
with open("edgetable.txt") as file:
    edge_table = file.read()
with open("sitetable.txt") as file:
    site_table = file.read()
with open("mutationtable.txt") as file:
    mut_table = file.read()
    
print("Successfully opened table files")
    
# Uses these tables in order to construct the corresponding tree sequence
# Note: once this overall file is working, as Ulises mentioned, it would be
# a good idea to do some experimentation for whether the sites/mutations 
# tables really need to be included.
treeSequence = tskit.load_text(io.StringIO(node_table), 
                               io.StringIO(edge_table), 
                               io.StringIO(site_table),
                               io.StringIO(mut_table),
                               strict = False)

print("Successfully constructed tree sequence from tables")

# Adds neutral mutations onto the resulting tree sequence
# Note: this rate should be changed later! Currently just a filler value- should
# be adjusted for greater accuracy related to program needs
treeSequence = msprime.sim_mutations(treeSequence, rate = 0.0001, model = msprime.InfiniteAlleles(), discrete_genome = True, keep = False)

print("Successfully simulated neutral mutations")

# Generates the site frequency spectrum
# Note: currently getting the folded site frequency spectrum- should this be changed, 
# or is it appropriate in this case? 
# Potentially do some experimentation with folded/unfolded once this code is actually functioning
siteFrequencySpectrum = treeSequence.allele_frequency_spectrum(span_normalise = False, polarised = True)
print(siteFrequencySpectrum)

# Uses this site frequency spectrum to generate the associated histogram
# Likely will need to change the labels on the histogram? Not sure exactly what they should be
# Note that this code is heavily derived from http://sesame.uoregon.edu/~adkern/stdpopsim/doc/tutorial.html 
plt.xlabel("Allele count", fontweight = "bold")
plt.ylabel("Number of sites", fontweight = "bold")
plt.bar(x=[x for x in range(0, 21)], height = siteFrequencySpectrum[:21], width = 0.6)
plt.tight_layout()
plt.savefig("alleleFrequencySpectrum.png")
