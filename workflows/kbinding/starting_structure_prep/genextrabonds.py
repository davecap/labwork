#!/usr/bin/python

import os
import optparse
#from configobj import ConfigObj, flatten_errors
import subprocess

VMD_PATH = "vmd"

TCL = """
set H102NE2 [[atomselect top "segname PEPA and resid 102 and name NE2"] get index]
set H102FE [[atomselect top "segname PEPA and resid 102 and name FE"] get index]
set H419NE2 [[atomselect top "segname PEPA and resid 419 and name NE2"] get index]
set H419FE [[atomselect top "segname PEPA and resid 419 and name FE"] get index]
set HSD411NE2 [[atomselect top "segname PEPA and resid 411 and name NE2"] get index]
set ASP412OD2 [[atomselect top "segname PEPA and resid 412 and name OD2"] get index]
set HSD333NE2 [[atomselect top "segname PEPA and resid 333 and name NE2"] get index]
set HSE284ND1 [[atomselect top "segname PEPA and resid 284 and name ND1"] get index]
set HSD334NE2 [[atomselect top "segname PEPA and resid 334 and name NE2"] get index]
set MG [[atomselect top "segname I and resid 1 and name MG"] get index]
set CUB [[atomselect top "segname I and resid 3 and name CU"] get index]
set CUA1 [[atomselect top "segname J and resid 1 and name CU"] get index]
set CUA2 [[atomselect top "segname J and resid 2 and name CU"] get index]
set C252SG [[atomselect top "segname PEPB and resid 252 and name SG"] get index]
set C256SG [[atomselect top "segname PEPB and resid 256 and name SG"] get index]
set H217ND1 [[atomselect top "segname PEPB and resid 217 and name ND1"] get index]
set H260ND1 [[atomselect top "segname PEPB and resid 260 and name ND1"] get index]
set M263SD [[atomselect top "segname PEPB and resid 263 and name SD"] get index]
set E254O [[atomselect top "segname PEPB and resid 254 and name O"] get index]
set E254OE1 [[atomselect top "segname PEPB and resid 254 and name OE1"] get index]

puts "PARSE #HDH102:NE2 <-> HDH102:FE 0.25"
puts "PARSE bond $H102NE2 $H102FE 50.0 2.5"
puts "PARSE #HDH419:NE2 <-> HDH419:FE 0.25"
puts "PARSE bond $H419NE2 $H419FE 200.0 2.5"
puts "PARSE #CUB <-> HDH419:FE"
puts "PARSE bond $CUB $H419FE 200.0 3.0"
puts "PARSE #HSD411:NE2 <-> MG 0.25"
puts "PARSE bond $HSD411NE2 $MG 10.0 2.1"
puts "PARSE #ASP412:OD2 <-> MG 0.25"
puts "PARSE bond $ASP412OD2 $MG 10.0 1.8"
puts "PARSE #E254:OE1 <-> MG"
puts "PARSE bond $E254OE1 $MG 10.0 1.9"
puts "PARSE #HSD333:NE2 <-> CuB 5 0.22"
puts "PARSE bond $HSD333NE2 $CUB 100.0 2.2"
puts "PARSE #HSE284:ND1 <-> CuB 5 0.22"
puts "PARSE bond $HSE284ND1 $CUB 100.0 2.2"
puts "PARSE #HSD334:NE2 <-> CuB 5 0.22"
puts "PARSE bond $HSD334NE2 $CUB 100.0 2.2"
puts "PARSE #CuA1 to CuA2 0.45"
puts "PARSE bond $CUA1 $CUA2 100.0 4.5"
puts "PARSE #CYS252:SG to CuA1 0.21"
puts "PARSE bond $C252SG $CUA1 500.0 2.1"
puts "PARSE #CYS252:SG to CuA2 0.21"
puts "PARSE bond $C252SG $CUA2 500.0 2.1"
puts "PARSE #CYS256:SG to CuA1 0.21"
puts "PARSE bond $C256SG $CUA1 500.0 2.1"
puts "PARSE #CYS256:SG to CuA2 0.21"
puts "PARSE bond $C256SG $CUA2 500.0 2.1"
puts "PARSE #HSE217:ND1 to CuA2 0.24"
puts "PARSE bond $H217ND1 $CUA2 500.0 2.4"
puts "PARSE #HSE260:ND1 to CuA1 0.243"
puts "PARSE bond $H260ND1 $CUA1 500.0 2.43"
puts "PARSE #MET263:SD to CuA2 0.255"
puts "PARSE bond $M263SD $CUA2 500.0 2.55"
puts "PARSE #GLU254:O to CuA1 0.235"
puts "PARSE bond $E254O $CUA1 500.0 2.35"
quit
"""

def main():    
    usage = """
        usage: %prog [options] <PSF> <PDB>
    """
    
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()
    
    if len(args) < 2:
        parser.error("Please specify input files")
    
    command = "vmd -dispdev none -psf %s -pdb %s" % (args[0], args[1])
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate(input=TCL)
    
    for line in stdoutdata.split('\n'):
        if line.startswith('PARSE'):
            line_parts = line.split(' ')
            print ' '.join(line_parts[1:])

if __name__=='__main__':
    main()


