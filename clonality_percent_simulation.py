#! /usr/bin/env python

""" Estimates the fraction of a given transformation that is composed of sister colonies, 
based on the percent of cells that had time to undergo one doubling during the recovery stage 
(it's assumed that none had time to double more than once) and the percent of cells plated 
(plating a lower percentage obviously cuts down on the number of repeats.)

N_samples is just the number of samples to run the simulation on - the default is 100000, 
which is large enough to give consistent results (within 1%) and small enough to run in 
a couple of seconds.

USAGE: clonality_percent_simulation.py percent_doubled percent_plated [N_samples] """

import random, sys

try:
    if len(sys.argv)==3:
        percent_doubled, percent_plated = [int(x) for x in sys.argv[1:]]
        N_samples = 100000
    elif len(sys.argv)==4:
        percent_doubled, percent_plated, N_samples = [int(x) for x in sys.argv[1:]]
    else:
        raise ValueError
except ValueError:
    print __doc__
    sys.exit("Error: Two or three arguments required!  All arguments must be integers.")

if not 0<=percent_doubled<=100 and 0<=percent_plated<=100:
    print __doc__
    sys.exit("Error: The first two arguments must be between 0 and 100% (percentages)")

if N_samples<0:
    print __doc__
    sys.exit("Error: The third argument must be higher than 0!")

sample_ID_list = []
for i in range(N_samples):
    sample_ID_list.append(i)
    x = random.randint(1,100)
    if x<=percent_doubled:
        sample_ID_list.append(i)

N_cells = len(sample_ID_list)
print "%s starting transformants, after %s%% cell doublings, yielded %s cells"%(N_samples, percent_doubled, N_cells)
cells_plated = []

for i in range(int(float(percent_plated)/100*N_cells)):
    x = random.randint(0,len(sample_ID_list)-1)
    cells_plated.append(sample_ID_list[x])
    del sample_ID_list[x]
    
total_cells = len(cells_plated)
unique_cells = len(set(cells_plated))
N_repeats = (total_cells - unique_cells) * 2
percent_repeats = 100*float(N_repeats)/total_cells
print "Plating %s%% of the cells yielded %s total cells, %s unique - %.1f%% of the cells are repeats."%(percent_plated, total_cells, unique_cells, percent_repeats)
