#! /usr/bin/env python

import sys, string, getopt
import stats

#-- Default settings as necessary
def usage():

  print 'usage:  window_avg [options...] <file_in >file_out'
  print '(  where [options] are:'
  print '(    -w <window size> for averaging the default window size is 7)'
  print '(    -h      {this help message}'
  sys.exit(1)

def get_options():
  w = 7
  file_in = sys.stdin
  try:
    opts,args = getopt.getopt(sys.argv[1:],'hw:',['window='])
  except :
    print 'Unrecognized Option'
    usage()
    return w

  if len(args) == 1:
    file_in = open(args[0])
  else:
    file_in = sys.stdin

  for o,a in opts:
    if o in ('-w', '--window'):
      w=int(a)
    elif o in ('-h'):
      usage()

  return w,file_in

####################################################################
# Any other options that should be put here?
####################################################################

# use list in the eventuality that other options may be added later
(window,input) = get_options()

if __name__ == '__main__':

  x = []

#  if len(sys.argv) > 1:
#    input = open(sys.argv[1])
#  else:
#    input = sys.stdin

  for line in input.readlines():
# ignore comments beginning with '#', or '@' (xvg files from gromacs) or blank lines
    if line[0] != '#' and line[0] != '@' and len(line) != 1:
      cols = string.split(line)
      x.append(map(float, cols))

  xavg = stats.window_avg(x,window)
  for i in range(len(xavg)):
#    print str(xavg[i])
    print "%12.6g" % xavg[i][0],
    for j in range(1,len(xavg[i])):
      print "%12.6g" % (xavg[i][j]),
    print
