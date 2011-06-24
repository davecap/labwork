#!/usr/bin/python

import os
import sys
import optparse
import numpy
from configobj import ConfigObj

def main():    
    usage = """
        usage: %prog [options] <config.ini>
    """
    
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()
    
    config = ConfigObj(args[0])
    project_path = os.path.abspath(os.path.dirname(config.filename))

    reps = []

    for r_id, r_config in config['replicas'].items():
        c = float(r_config['coordinate'])
        # for each distance, calculate the average and stdev of the difference from the restraining coordinate
        try:
            f = open(os.path.join(project_path,r_id,'distances'), 'r')
        except:
            print "Couldn't open distances file for %s" % r_id
        else:
            coords = []
            for r in f.readlines():
                coords.append(c-float(r))
            d = numpy.array(coords)
            reps.append((r_id, c, d))

    for r in sorted(reps, key=lambda r: r[1]):
        r_id = r[0]
        c = r[1]
        d = r[2]
        sem = numpy.std(d, ddof=1)/numpy.sqrt(d.size)
        if c < -35 and c > -40:
            print "* %f [%s]: %f +- %f" % (c, r_id, d.mean(), sem)
        else:
            print "%f [%s]: %f +- %f" % (c, r_id, d.mean(), sem)
        
if __name__=='__main__':
    main()
