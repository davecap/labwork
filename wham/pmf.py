import sys
import numpy
import subprocess

def run_wham(min, max, bins, metafilepath, temp=315.0, tol=0.0001, **kwargs):
    outfile = "%s_wham.out" % metafilepath
    sys.stderr.write("Running WHAM: %s -> %s\n" % (metafilepath, outfile))
    command = "wham %f %f %d %f %f 0 %s %s" % (min, max, bins, tol, temp, metafilepath, outfile)
    sys.stderr.write("\t%s\n" % (command))
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate()
    sys.stderr.write("\tDone running WHAM -> %s\n" % (outfile))

def process_pmf(filename, shift=False):
    x = []
    y = []
    infile = open(filename, 'r')
    while infile:
        line = infile.readline().strip()
        if len(line) == 0:
            break
        elif line.startswith('#'):
            continue
        split_line = [ s.strip() for s in line.split('\t') ]
        if len(split_line) > 3 and split_line[1] != 'inf':
            x.append(float(split_line[0]))
            y.append(float(split_line[1]))
    
    if shift:
        # Calculate the PMF shift by averaging the first 10 bins (bulk water dG)
        s = 1.0*numpy.array(y)[0:10].mean()
        sys.stderr.write("PMF shift by %0.2f\n" % s)
        y = [ (i-s) for i in y ]
    
    return zip(x,y)
        
def dG_bind(pmf, imin=0, imax=1):
    """ Calculate the relative binding free energy for a PMF given bounds """
    # dGbind = -(1/kB)
    kB = 0.00198721 # kcal/mol/K
    t = 315.0 # K
    Beta = kB*t
    
    x_prev = None
    dG_prev = None
    s = 0.0
    
    for x,dG in pmf:
        if x > imax:
            break
        elif x > imin and dG_prev is not None:
            # e^{-B*[(dG+dGprev)/2.0]}*(x-x_prev)
            s += numpy.exp((-1.0/Beta)*(dG+dG_prev)/2.0)*(x-x_prev)
        dG_prev = dG
        x_prev = x
    
    dGbind = (-1.0*Beta)*numpy.log(s)
    return dGbind


