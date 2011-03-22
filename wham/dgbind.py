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
    parser.add_option("--bind-min", dest="bind_min", type="float", default=10, help="Lower bound for binding [default: %default]")
    parser.add_option("--bind-max", dest="bind_max", type="float", default=20, help="Upper bound for binding [default: %default]")
    parser.add_option("--std-min", dest="std_min", type="float", default=0, help="Lower bound for standard state [default: %default]")
    parser.add_option("--std-max", dest="std_max", type="float", default=10, help="Upper bound for standard state [default: %default]")
    parser.add_option("--autoshift", dest="autoshift", default=True, action="store_false", help="Auto-shift PMF [default: %default]")
    
    (options, args) = parser.parse_args()
    
    pmf = process_pmf(args[0], shift=options.autoshift)
    print dG_bind(pmf, imin=options.bind_min, imax=options.bind_max)


if __name__=='__main__':
    main()
