#!/usr/bin/python

import os
import sys
import optparse
import numpy

from pmf import process_pmf, dG_bind

def main():    
    usage = """
        usage: %prog [options] <PMF>

        Calculate dG_bind
    """
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-bm", "--bound-min", dest="bound_min", type="float", default=35, help="Lower bound for binding [default: %default]")
    parser.add_option("-bx", "--bound-max", dest="bound_max", type="float", default=25, help="Upper bound for binding [default: %default]")
    parser.add_option("-um", "--unbound-min", dest="unbound_min", type="float", default=0, help="Lower bound for standard state [default: %default]")
    parser.add_option("-ux", "--unbound-max", dest="unbound_max", type="float", default=10, help="Upper bound for standard state [default: %default]")
    parser.add_option("--autoshift", dest="autoshift", default=True, action="store_false", help="Auto-shift PMF [default: %default]")
    parser.add_option("-s", "--std", dest="standard", default=True, action="store_false", help="Do standard state correction [default: %default]")
    
    (options, args) = parser.parse_args()
    
    pmf = process_pmf(args[0], shift=options.autoshift)
    
    correction = 0.0
    if options.standard:
        correction = dG_bind(pmf, imin=options.unbound_min, imax=options.unbound_max)
        print "Correction: %0.2f" % correction

    print dG_bind(pmf, imin=options.bound_min, imax=options.bound_max)-correction

if __name__=='__main__':
    main()
