# 286
set E286 [ atomid PEPA 286 CA ]
addatom $E286

set N139 [ atomid PEPA 139 CA ]
addatom $N139

set ION [ atomid POT 4712 POT ]
addatom $ION

print "K=$K"
print "COORD=$COORD"
    
# Call procedure 
proc calcforces {} { 
    global HSE10CA
    global HSE7CA
    global HSE549CA
    global ION 
    global E286
    global N139
    
    loadcoords c 
    
    global K
    global COORD

    # current ion coordinates
    set cur_x [lindex $c($ION) 0]
    set cur_y [lindex $c($ION) 1]
    set cur_z [lindex $c($ION) 2]
    #print "CUR XYZ $cur_x $cur_y $cur_z"
    
    set E286CA [ lindex $c($E286) 2 ]
    set ZDIST [ expr $cur_z-$E286CA ]
    #print "ZDIST $ZDIST"

    # only apply cylinder restraint if ion is near or outside the entrance to the D-channel
    if {$ZDIST < -24.0} {
        set N139CAx [ lindex $c($N139) 0 ]
        set N139CAy [ lindex $c($N139) 1 ]
        set N139CAz [ lindex $c($N139) 2 ]

        #print "Applying cylinder restraint!"
        # average X and Y distances for the histidines

        # ion distance from the center
        set ion_dist [ expr sqrt([ expr ($cur_x-$N139CAx)*($cur_x-$N139CAx) + ($cur_y-$N139CAy)*($cur_y-$N139CAy)])]
        
        # only apply the restraint if the ion leaves the cylinder radius
        if {$ion_dist > 4} {
            #print "ION OUTSIDE CYLINDER"
            # calculate force to apply to ion ($cur_x, $cur_y) -> ($avg_x, $avg_y)
            set k 50
            set r12 [ vecsub "$N139CAx $N139CAy $cur_z" "$cur_x $cur_y $cur_z" ]
            set force [ expr $k*($ion_dist) ]
            addforce $ION [ vecscale $force $r12 ]
            #print "ION XY force $force"
        }
    }

    # Calculate the force along Z for the ion
    #COORD is the Z distance from the E286 Calpha
    set force [ expr -$K*($ZDIST-$COORD) ]    
    addforce $ION [vecscale $force { 0 0 1 }]
    #print "ION Z force $force"

}

