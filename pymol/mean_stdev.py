#! /usr/bin/python
import string, sys, math, os

i = 0  # line counter

sum = [0.] * 100
sum2 = [0.] * 100
max = [-99999999.] * 100
min = [999999999.] * 100
mean = [0.] * 100
stdev = [0.] * 100

for line in sys.stdin.readlines():
	if line[0] != "#":
		cols = string.split(line)
		x = map(float, cols)
		n = len(x)
		#print n, range(n)

####################################################################
#        # subtract column value from previous column to get difference
#        for j in range(n-1,0,-1):
#                x[j] = x[j] - x[j-1]
####################################################################

		for j in range(n):
			sum[j] = sum[j] + x[j]
			sum2[j] = sum2[j] + x[j]*x[j]
			if x[j] > max[j]: max[j] = x[j]
			if x[j] < min[j]: min[j] = x[j]

		i = i + 1

num = i

print "Num pts: ", num, "Num cols: ", n

#print "col:    Mean:    stdev:      Max:     Min:"
for j in range(n):
	mean[j] = sum[j]/num
	stdev[j] = math.sqrt( (sum2[j] - mean[j]*sum[j])/num)
#	print "%4d  %8.4f  %8.4f  %8.4f  %8.4f" % \
#	(j, mean[j], stdev[j], max[j], min[j])

# comment out the following, and uncomment the two commented lines above
#in order to change the output direction.

print 'Column',
for j in range(n):
	print "%8d " % j,

print '\nMean   ',
for j in range(n):
	print "%8.4f " % mean[j],
print '\nStd dev',
for j in range(n):
	print "%8.4f " % stdev[j],
print '\nMax    ',
for j in range(n):
	print "%8.4f " % max[j],
print '\nMin    ',
for j in range(n):
	print "%8.4f " % min[j],
