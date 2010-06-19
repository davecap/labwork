# Read a topology file and convert to another format
#
# Note: only supports reading CHARMM RTF files and writing to Gromacs RTP files
#
# Author: David Caplan <david@davidacaplan.com>
# Date: December 2009

from pyparsing import *
import string
import optparse
import sys

from molecule import *
from topology import *

def main():
    usage = """
    CHARMM RTF File Parser
    usage: %prog <input file>
    """
    
    parser = optparse.OptionParser(usage)
    # parser.add_option("-x", dest="col_x", default=0, help="X Column [default: %default]")
    options, args = parser.parse_args()
    if len(args) == 0:
       parser.error("No input file specified")
    
    #f = open(args[0])
    
    try:
        print "Parsing RTF File: %s" % args[0]
        topology = CharmmTopology()
        topology.read(args[0])
        topology.summary()

        #gt = GromacsTopology()
        #gt.write(topology, 'test')
        
        #parser = CHARMMRTFParser(args[0])
        #fields = parser.parse()
        #parser.summary()
    except ParseException, err:
        print err.line
        print " "*(err.column-1) + "^"
        print err
    
    print "Done."


test_RTF_data = """
*  Test RTF file
*   David Caplan 2009
*

! this is a test comment

MASS     1 H      1.00800 ! poop
! this is a test regular comment

MASS    11 C     12.01100
! random stuff
! 
MASS    13 CH2E  14.02700
MASS    14 CH3E  15.03500

DECL -C ! test comment
DECL +N

DEFA FIRS NTER LAST CTER

! this is a test comment
RESI ALA     0.00000 ! Alanine residue
GROU
ATOM N    NH1    -0.35
ATOM CA   CH1E    0.10
GROU
ATOM CB   CH3E    0.00
ANGL -C   N    CA             N    CA   C              CA   C    +N
DIHE -C   N    CA   C         N    CA   C    +N        CA   C    +N   +CA
IMPH N    -C   CA   H         C    CA   +N   O         CA   N    C    CB
CMAP -C  N  CA  C   N  CA  C  +N
DONO H    N    -C   CA
ACCE O C
BILD -C   CA   *N   H      0.0000    0.00  180.00    0.00   0.0000  ! test1a fhello (test) :) ! |\ []

RESI OH2     0.00000
GROU
ATOM OH2  OH2    -0.40000
ATOM H2   H       0.20000
BOND OH2  H1        OH2  H2
THET H1   OH2  H2
DONO H2   OH2 -O -O
ACCE OH2
PATC FIRS NONE LAST NONE

END
"""

if __name__ == '__main__':
    main()
