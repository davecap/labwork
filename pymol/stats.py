#!/usr/bin/env python
# Copyright (c) 2003-2004 Robert L. Campbell
import sys, os, getopt
import scipy

import numpy as N

def stats(x):
  """
  returns an array containing:
  mean(x),std(x),median(x),max(x),min(x)
  """
  x = N.asarray(x)
  if len(x) != 0:
    return N.mean(x,0),N.std(x,0,ddof=1),N.median(x,0),N.max(x,0),N.min(x,0)
  else:
    return None,None,None,None,None

def avg_dev(x):
  """
  calculate the average deviation of the data in array x
  see Numerical Recipes, Chapter 13

  returns a new array with one less dimension than [x]
  """
  x = N.asarray(x)
  my_mean = N.mean(x,0)
  return N.multiply(1./len(x),N.add.reduce(N.fabs(N.subtract(x,my_mean))))

def var(x):
  """
  calculate the variance (2nd moment) of the data in array x
  see Numerical Recipes, Chapter 13

  returns a new array with one less dimension than [x]

  Could just use numpy.var(x,0,ddof=1)
  """
  x = N.asarray(x)
  diff = N.subtract(x,N.mean(x,0))
  #return N.multiply(std(x),std(x))
# could have used the above, but the following is faster for lots of data
  return N.multiply(1./(len(x)-1.0),N.add.reduce(diff**2))

def skew(x):
  """
  calculate the skewness (3rd moment) of the data in array x
  see Numerical Recipes, Chapter 13

  returns a new array with one less dimension than [x]
  """
  x = N.asarray(x)
  st = N.divide(N.subtract(x,N.mean(x,0)),N.std(x,0,ddof=1))
  return N.add.reduce(st**3)/len(x)

def kurtosis(x):
  """
  calculate the kurtosis (4th moment) of the data in array x
  see Numerical Recipes, Chapter 13

  returns a new array with one less dimension than [x]
  """
  x = N.asarray(x)
  st = N.divide(N.subtract(x,N.mean(x,0)),N.std(x,0,ddof=1))
  return N.subtract(N.add.reduce(st**4)/len(x),3)

#def kurt(x):
#  """ formula found on wikipedia
#  doesn't seem to work for me
#  """
#  x = N.asarray(x)
#  n = len(x)
#  k = N.add.reduce(N.subtract(x,N.mean(x,0))**4)/(N.var(x,0,ddof=1)**2)
#  G2 = (((n+1)*n) / ((n-1)*(n-2)*(n-3)))  * k - (3*(n-1)**2)/((n-2)*(n-3)) 
#  return G2
#
def mode(x,j):
  """
  calculate the mode for continuous data in array x
  see Numerical Recipes, Chapter 13

  usage:  index_list, probability_list = mode(array_of_data,window)

  returns two lists: 
    1) the index {i.e. the value from the data calculated as (x[i]+x[i+window])/2}
    2) the probability of finding that value

  """
# make sure data is in an array and make sure it is sorted 
# (will not maintain synchronicity between columns though, but that shouldn't matter
# for the mode calculation!)
  x = N.asarray(x)
  x = N.msort(x)

# create the index array
  ind = N.zeros((len(x)-j,x.shape[1]),float)
# create the probability array

  p = N.zeros((len(x)-j,x.shape[1]),float)
  n=len(x)
  for i in range(n-j):
    ind[i] = N.multiply(0.5,add(x[i],x[i+j]))
    p[i] = N.divide(j,N.multiply(n,N.subtract(x[i+j],x[i])))
  return ind, p

def histogram(x,numbin=10,binwidth=0, binmin=0,binmax=0):
  mean,stdev,median,max,min = stats(x)

  try:
    num_data_sets = len(x[0])
  except TypeError:
    num_data_sets = 1

# don't bother with arrays here
  if num_data_sets == 1:
    if binwidth == 0:
      interval = ((max-min)/float(numbin))
    else:
      interval = binwidth
      if binmin ==0 and binmax == 0:
        numbin = int((max-min)/interval)
      else:
        numbin = int((binmax-binmin)/interval)
        min=binmin
        max=binmax

# use arrays to store values for min,max,numbin,intervals, etc.
  else:
    pass
#    if binwidth == 0:
#      interval = zeros((num_data_sets),float)
#      for j in range(num_data_sets):
#        print j
#        interval[j] = (max[j]-min[j])/float(numbin)
#    else:
#      numbin = zeros((num_data_sets),float)
#      interval = zeros((num_data_sets),float)
#      for j in range(num_data_sets):
#        interval[j] = binwidth
#        if binmin ==0 and binmax == 0:
#          numbin[j] = max[j]-min[j]/interval[j]
#        else:
#          numbin[j] = binmax-binmin/interval[j]

  print numbin,num_data_sets
  bin = N.zeros((numbin,num_data_sets),float)


# look in color_b.py for ideas about histrograms
#  for i in range(len(x)):
#    for j in range(len(x[i])):
#      for k in range(numbin):
#        start = min[j] + k*interval[j]
#        end = min[j] + (k+1)*interval[j]
#        if x[i][j] >= start and x[i][j] < end:
#          bin[k] += 1
#
#  for k in range(numbin):
#    for j in range(len(x[i])):
#      print min[j] + k * interval[j], min[j] + (k+1)*interval[j], bin[k]
  for i in range(len(x)):
    for k in range(numbin):
      start = min + k*interval
      end = min + (k+1)*interval
      if x[i] >= start and x[i] < end:
        bin[k] += 1

  histogram = []

  for k in range(numbin):
#    print min + k * interval, min + (k+1)*interval, (min + k * interval + min + (k+1)*interval)/2, bin[k]
    histogram.append((min+k*interval,min+(k+1)*interval,bin[k]))
  return min,max,interval,histogram,bin

def histo(x,numbin=10,binwidth=0, binmin=None,binmax=None):
  return scipy.stats.histogram(x,numbins=numbin,defaultlimits=(binmin,binmax),printextras=True)


# least-squares fit of two arrays, x and y of equal length
def lsq(x,y):
  """
  returns slope,intcpt,cc,err_slope,err_intcpt,sigma,prob_err,sx,sxx
          cc=correlation coeff sigma=RMS of residuals, prob_err=probable error
          sx=SUM(x), sxx=SUM(x**2) (used in calculating errors in calculated y values)
  """

  if len(x) != len(y):
    sys.stderr.write("ERROR: unequal lengths of x and y arrays: %d %d" % (len(x),len(y)))
    sys.exit()
  else:
    n = len(x)
# initialize sums to 0
  sx = 0 
  sxx = 0 
  sy = 0 
  syy = 0 
  sxy = 0

  for i in range(n):
    # set x value to be counter and y value to be coordinate value
    # x[i] = q[i] y[i] = p[i];
    # do sums 
    sx += x[i]
    sy += y[i]
    sxx += x[i]*x[i]
    syy += y[i]*y[i]
    sxy += x[i]*y[i]

  # calculate slope, intercept, corr. coeff., and errors.

  det = n*sxx - sx*sx
  slope = (n*sxy - sx*sy)/det
  intcpt = (-sx*sxy+sxx*sy)/det
  cc = (intcpt*sy + slope*sxy - sy*sy/n) / (syy - sy*sy/n)
  a = slope*slope*sxx + n*intcpt*intcpt + syy + 2*(slope*intcpt*sx-intcpt*sy-slope*sxy)
  err_slope = N.sqrt( (n*a) / ((n-2)*det) )
  err_intcpt = err_slope*N.sqrt(sxx/n)
  sigma = N.sqrt(a/(n-2))
  prob_err = sigma*2./3.

  return slope,intcpt,cc,err_slope,err_intcpt,sigma,prob_err,sx,sxx
# least-squares fit of two arrays, x and y of equal length

def lsq_noerrors(x,y):
  """
  returns slope,intcpt
  """

  if len(x) != len(y):
    sys.stderr.write("ERROR: unequal lengths of x and y arrays: %d %d" % (len(x),len(y)))
    sys.exit()
  else:
    n = len(x)
# initialize sums to 0
  sx = 0 
  sxx = 0 
  sy = 0 
  syy = 0 
  sxy = 0

  for i in range(n):
    # set x value to be counter and y value to be coordinate value
    # x[i] = q[i] y[i] = p[i];
    # do sums 
    sx += x[i]
    sy += y[i]
    sxx += x[i]*x[i]
    syy += y[i]*y[i]
    sxy += x[i]*y[i]

  # calculate slope, intercept, corr. coeff., and errors.

  det = n*sxx - sx*sx
  slope = (n*sxy - sx*sy)/det
  intcpt = (-sx*sxy+sxx*sy)/det
  #cc = (intcpt*sy + slope*sxy - sy*sy/n) / (syy - sy*sy/n)
  #a = slope*slope*sxx + n*intcpt*intcpt + syy + 2*(slope*intcpt*sx-intcpt*sy-slope*sxy)
  #err_slope = math.sqrt( (n*a) / ((n-2)*det) )
  #err_intcpt = err_slope*math.sqrt(sxx/n)
  #sigma = math.sqrt(a/(n-2))
  #prob_err = sigma*2./3.

  return slope,intcpt

def print_stats(x,direction='horizontal'):
  my_mean, my_stdev, my_median, my_max, my_min = stats(x)

  if my_mean != None:
    numpts = len(x)         # normally use 'i' to count this
    numcol = len(x[0])      # normally use 'j' to count this

    print 'Num pts: ', numpts, 'Num cols: ', numcol

    if direction == 'horizontal':
      print 'col:    Mean:       Stdev:      Median:      Min:        Max:'
      for j in range(numcol):
        print '%4d  %10.8g  %10.8g  %10.8g  %10.8g  %10.8g' % \
          (j, my_mean[j], my_stdev[j], my_median[j], my_min[j], my_max[j])

    else:
      print 'Column  ',
      for j in range(numcol):
        print '%8d ' % j,

      print '\nMean   ',
      for j in range(numcol):
        print '%8.8g ' % my_mean[j],
      print '\nStd dev',
      for j in range(numcol):
        print '%8.8g ' % my_stdev[j],
      print '\nMedian ',
      for j in range(numcol):
        print '%8.8g ' % my_median[j],
      print '\nMax    ',
      for j in range(numcol):
        print '%8.8g ' % my_max[j],
      print '\nMin    ',
      for j in range(numcol):
        print '%8.8g ' % my_min[j],
  else:
    print "### No data ###"

def window_avg(x,window):
  print "# Averaged over window size: %3d" % window
  x = N.asarray(x)
  half_win = window/2
  i=0
  num_data = len(x)
  num_cols = len(x[0])
  x_avg = N.zeros((num_data-window+1,num_cols),float)

  for i in range(num_data - window + 1):
    total = N.zeros((1,num_cols),float)
    x_avg[i][0] = x[i+half_win][0]
    for j in range(window):
      for k in range(1,num_cols):
        total[0][k] += x[i+j][k]
#        print i,j,k,x[i+j][k]
#      print total[0]
      for k in range(1,num_cols):
        x_avg[i][k] = total[0][k]/float(window)
  return x_avg

def window_lsq(x,window):
  print "# Slopes over window size: %3d" % window
  x = N.asarray(x)
  half_win = window/2
  i=0
  num_data = len(x)
  num_cols = len(x[0])
  slopes = N.zeros((num_data-window+1,num_cols),float)

  for i in range(num_data - window + 1):
    slopes[i][0] = x[i+half_win][0]
    slopes[i][1] = lsq_noerrors(x[i:i+window,0],x[i:i+window,1])[0]
  return slopes

if __name__ == '__main__':

  try:
    opts,args = getopt.getopt(sys.argv[1:],'f:h',['help','file=','horiz','vert','horizontal','vertical'])
  except:
    sys.stderr.write("\n*********************************\n")
    sys.stderr.write("\n      Unknown options %s\n" % str(sys.argv[1:]))
    sys.stderr.write("\n*********************************\n\n")
    usage()

  if len(args) == 1:
    input = open(args[0])
  else:
    input = sys.stdin

  direction = 'horizontal'

  for o,a in opts:
    if o in ('-h', '--help'):
      usage()
    elif o in ('-f','--file'):
      file_in = a
    elif o in ('--horiz','--horizontal'):
      direction = 'horizontal'
    elif o in ('--vert','--vertical'):
      direction = 'vertical'
  x = []

  for line in input.readlines():
# ignore comments beginning with '#' or blank lines
    if line[0] != '#' and len(line) != 1:
      cols = line.split()
      x.append(map(float, cols))

  y = N.asarray(x)
  print_stats(y,direction=direction)
