# Requires VMD.
#vmd -e set_dihedral.tcl -dispdev none -pdb <initial.pdb> -args <dihedral angle>
# mol load psf 1m56_wt.psf pdb 1m56_wt.pdb

set movesel [ atomselect top "segname PEPA and resid 139 and (name CG or name OD1 or name ND2 or name HD21 or name HD22 or name HB1 or name HB2)" ]
set N139N [ atomselect top "segname PEPA and resid 139 and name N" ]
set N139CA [ atomselect top "segname PEPA and resid 139 and name CA" ]
set N139CB [ atomselect top "segname PEPA and resid 139 and name CB" ]
set N139CG [ atomselect top "segname PEPA and resid 139 and name CG" ]

set ind1 [ $N139N get index ]
set ind2 [ $N139CA get index ]
set ind3 [ $N139CB get index ]
set ind4 [ $N139CG get index ]

set tmpmolid 0

foreach i $argv {
    set NEWDIHED $i
    
    set dihedral [measure dihed [list [list $ind1 $tmpmolid] [list $ind2 $tmpmolid] [list $ind3 $tmpmolid] [list $ind4 $tmpmolid] ]]
    set bsel1 [atomselect $tmpmolid "index $ind2"]
    set bsel2 [atomselect $tmpmolid "index $ind3"]
    set delta [expr -1 * ($dihedral - $NEWDIHED)]
    set mat [trans bond [lindex [$bsel1 get {x y z}] 0] [lindex [$bsel2 get {x y z}] 0] $delta deg]
    $movesel move $mat
    $bsel1 delete
    $bsel2 delete

    set all [ atomselect top all ]
    $all writepdb ./dihed_$NEWDIHED.pdb
}

quit

