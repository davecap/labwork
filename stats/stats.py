#!/usr/bin/python

import fileinput
import numpy
import scipy.stats.stats
import random

d = []

for line in fileinput.input():
    if line:
        d.append(float(line))

#n = int(len(d)/3)
#block_averages = []
#for i in range(0,50):
#    random.shuffle(d)
#    block_averages.append(numpy.mean(d[0:n]))
#    #print numpy.mean(d[0:n])

#print numpy.std(block_averages)
print numpy.std(d)
#nd = numpy.array(d)
#nba = numpy.array(block_averages)

#print scipy.stats.stats.sem(nba)
