This is a fork of the Masel Lab's forward time simulation, designed for studying questions
related to mutation load and demographic inference within population genetics. It simulates
the genomes of individuals in a population over time by splitting the genome into discrete
linkage blocks, with each linkage block containing a value that summarizes the effects of 
all of the mutations that have taken place in that part of the genome. 

Within my fork of the code, functionality was added to allow for runs involving recessivity
to be made, rather than the current model of multiplicative co-dominance. This functionality
will then be used in order to investigate the effects of bottlenecks on demographic inference.
