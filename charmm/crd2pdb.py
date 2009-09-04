#!/usr/bin/python

# This script converts a CHARMm crd file to a PDB file. (This has been done 
# several times before - hopefully this is more successful :-)
#
# EPF 27-03-2000
#

import sys,string,struct,os

if len(sys.argv) < 2:
    print "usage: crd2pdb crd-file > pdb-file"
    os._exit(0)

crdfile = sys.argv[1]
f=open(crdfile) 
print 'REMARK   Converted from '+crdfile

buffer = f.readline()
while buffer[0:1] == "*":
    buffer = f.readline()
    numofatoms = float(buffer)
    
    j = 0;
    while buffer != "":
        buffer = f.readline()
        if buffer != "":
            atomnr   = int(buffer[0:5]) 
            resnr    = int(buffer[5:10]) 
            residue  = buffer[11:15] 
            atom     = buffer[16:20] 
            x         = float(buffer[20:30]) 
            y       = float(buffer[30:40]) 
            z       = float(buffer[40:50]) 
            segment  = buffer[51:55]
            reschain = buffer[56:63]     
            bfactor  = float(buffer[63:70])  

            if string.find(string.letters,reschain[0:1]) != -1:
                chain = reschain[0:1]
            else:
                chain = " "
            
            if (resnr<10000):
                atomdata =(atomnr,atom,residue,chain,resnr,x,y,z,0.0,bfactor) 
                pdbline = 'ATOM  %5i  %3s%4s%c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f' % atomdata 
            else:
                atomdata =(atomnr,atom,residue,resnr,x,y,z,0.0,bfactor) 
                pdbline = 'ATOM  %5i  %3s%4s%5i    %8.3f%8.3f%8.3f%6.2f%6.2f' % atomdata 

            print pdbline
            j = j + 1
    
    print "REMARK There are "+str(j)+" atoms in this file"
    if j != numofatoms:
        print "REMARK Actual number of atoms does not match the number given in crd file ("+str(numofatoms)+")"
