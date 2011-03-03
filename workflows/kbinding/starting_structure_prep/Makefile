# Remove ions from initial structure
md.psf md.pdb: initial.psf initial.pdb
	python removeions.py initial.psf initial.pdb

# Generate extrabonds
extrabonds: md.psf md.pdb
	python genextrabonds.py md.psf md.pdb > extrabonds

# run the minimization and equilibration
md.dcd: extrabonds md.psf md.pdb
	mpirun -np 8 ~/bin/namd2 namd2.conf >> md.log

# extract last frame from DCD
final.pdb: md.dcd md.psf
	catdcd -o final.pdb -otype pdb -stype pdb -s md.pdb -first 1 -last 1 md.dcd 2>&1 >> /dev/null

clean:
	rm -rf extrabonds
	rm -rf final.*
	rm -rf md.*
	rm -rf *.log
	rm -rf *.txt
	rm -rf rst.*

.PHONY: clean

default: final.pdb
