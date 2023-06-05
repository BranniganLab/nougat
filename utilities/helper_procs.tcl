# ConvertRadianToDegree (Previously: RtoD)--
#       Converts radians to degrees
#
# Arguments:
#       r         {float}    Angle in radians
#
# Results:
#       Returns angle in degrees
proc RtoD {r} {
    global M_PI
    return [expr $r*180.0/$M_PI]
}

# fitVecsToSel (Previously: tilt_angles)--
#
#       Performs a least squares fit for each lipid tail of N bead length
#       given the x, y and z values in order. This allows you to use 1 instance 
#       of $lsqnormfactor only,rather than calculate it on the fly each time.
#
# Arguments:
#       length      {int}     Number of beads in a single lipid tail.
#       xvals       {list}    Ordered list of x coordinates. 
#       yvals       {list}    Ordered list of y coordinates. 
#       zvals       {list}    Ordered list of z coordinates.
#
# Results:
#       The result is a list of list that contains a normalized fitted vector 
#       and the number of beads in an acyl chain. 
#
#       ex. 
#       {{{0.534522 -0.801783 0.267261} {0.534522 -0.801783 0.267261} {0.534522 -0.801783 0.267261} {0.534522 -0.801783 0.267261}}...}   

proc fitVecsToSel {length xvals yvals zvals} {
    set tiltList []
    set lsqNormFactor [calculateLsqNormFactor $length]
    set xvec [fitTailVectors $length $xvals $lsqNormFactor]
    set yvec [fitTailVectors $length $yvals $lsqNormFactor]
    set zvec [fitTailVectors $length $zvals $lsqNormFactor]
    for {set i 0} {$i < [llength $xvec]} {incr i} {
        set vector "[lindex $xvec $i] [lindex $yvec $i] [lindex $zvec $i]"
        set norm [vecnorm $vector]

        # it is desirable to have a list that is the same length as $xvals
        # for ease of indexing in the binning process; hence lrepeat
        lappend tiltList [lrepeat $length $norm]
    }

    # this is now a list of lists of lists, but we want a list of lists
    set finalTiltList [concatenateList $tiltList "NULL"]

    return $finalTiltList
}

# calculateLsqNormFactor (Previously: calc_lsq_normfactor)-- 
#       
#       calculates the normalization factor for a least square fitting.
#
# Arguments:
#       length      {int}     number of beads in the least square fitting.
#       
# Results: 
#       The result is a list of normalization factors of the form (i-(N-1)/2).
#       
#       ex. 
#       >>> calculateLsqNormFactor 4
#       {-1.5 -0.5 0.5 1.5} 

proc calculateLsqNormFactor { length } {
    set diff [expr $length-1]
    set d [expr 0.5*$diff] ;#normalization factor
    set I []

    #Make a list from 0 to length
    for {set k 0} {$k < $length} {incr k} {
        lappend I $k
    }

    set lsqNormFactor [vecexpr $I $d sub]

    return $lsqNormFactor
}

# fitTailVectors (Previously: fit_all_tails)-- 
#
#       Fits all points x to x = ai + b, i=0...N-1
#
# Arguments: 
#       tailLength        {int}    number of beads in a tail
#       listOfTailCoords  {list}   list of x, y or z coordinates for all beads in a tail
#                                  for all lipids to be evaluated
#       lsqNormFactor     {list}   list of normalization factors (see: calculateLsqNormFactor)  
# 
# Results:
#       returns the value a = sum[ (i-(N-1)/2) * x_i] for each tail ; reference: Bevington

proc fitTailVectors {tailLength listOfTailCoords lsqNormFactor} {
    set fitValues []
    set NumLipids [expr [llength $listOfTailCoords]/$tailLength]
    set startIdx 0

    # iterate through list of all tail coords, separating them into one 
    # lipid tail at a time with lrange
    for {set i 0} {$i < $NumLipids} {incr i} {
        set startIdx [expr $i*$tailLength]
        set endIdx [expr $startIdx+$tailLength-1]
        set coords [lrange $listOfTailCoords $startIdx $endIdx]

        # Perform least squares fitting:
        # Multiply x_i and the differences (i-(N-1)/2). stack=1vec
        # Sum over the vector. stack=1scalar
        # Append to $fit_values
        lappend fitValues [vecexpr [vecexpr $coords $lsqNormFactor mult] sum]
    }
    return $fitValues
}

# concatenateList (Previously: cat_list)--
#
#       Concatenates list elements into one long string, separated by 
#       spaces on either side of the $delimiter.
# 
# Arguments:
#       inputList       {List}      list of elements to concatenate
#       delimiter       {str}       string to be input between items in inputList
#
# Results:
#       Returns a list of elements with the specified delimiter between "NULL" will 
#       result in no delimiter with a single space separating elements."or" will 
#       also enclose list items in parentheses for use as atomselection text.
#
#       ex. 
#       >>>  concatenateList {1 2 3 4 5} $
#       {1 $ 2 $ 3 $ 4 $ 5}

proc concatenateList {inputList delimiter} {
    if {$delimiter eq "or"} {
        set output "([lindex $inputList 0])"
    } else {
        set output [lindex $inputList 0]
    }
    
    for {set i 1} {$i < [llength $inputList]} {incr i} {
        if {$delimiter eq "or"} {
            set element "([lindex $inputList $i])"
            set output "${output} or ${element}"
        } elseif {$delimiter eq "NULL"} {
            set output "${output} [lindex $inputList $i]"
        } else {
            set output "${output} $delimiter [lindex $inputList $i]"
        }
    }
    return $output
}

# findHeadsAndTails (Previously: heads_and_tails)--
#
#       Finds the first and last bead for each lipid tail for all lipids
#       to be evaluated.
#
# Arguments: 
#       tailList        {List}      List of list of atomselections of lipid tails 
#   
# Results:
#       Returns a list of lists containing the starting beads and 
#       ending beads for a given lipid's acyl chains.
# E.G. POPC start beads would be "C1A C1B" and end beads would be "C4A C4B"

proc findHeadsAndTails { tailList } {
    for {set i 0} {$i < [llength $tailList]} {incr i} {
        set startBead []
        set endBead []
        foreach tail [lindex $tailList $i] {
            lappend startBead [lindex $tail 0]
            lappend endBead [lindex $tail end]
        }
        lappend startSelList [concatenateList $startBead "NULL"]
        lappend endSelList [concatenateList $endBead "NULL"]
    }

    return [list $startSelList $endSelList]
}

# rotateSystem (Previously: rotate_system)--
#
#       Rotates the entire system by a certian degree over frames to be 
#       evaluated
#
# Arguments: 
#       axis        {str}       axis of rotation(x, y, z)
#       degree      {int}       degree of rotation from 0 to 360
#       start       {int}       start frame
#       stop        {int}       stop frame
#
# Results:
#       results in the system being rotated by user-defined specifications

proc rotateSystem {axis degree start stop} {
    if {$stop == -1} {
        set stop [molinfo top get numframes]    
    } 
    set sel [atomselect top all]
    for {set i $start} {$i<$stop} {incr i} {
        $sel frame $i 
        $sel update 
        set com [measure center $sel weight mass]
        set matrix [transaxis $axis $degree]
        $sel moveby [vecscale -1.0 $com]
        $sel move $matrix
        $sel moveby $com
    }
    $sel delete
}

# assignLeaflet (Previously: leaflet_check)--
#       
#       Checks whether lipid tails are above/below the PO4 bead,and assigns user to 
#       1 or 2 for outer or inner leaflet.
# 
# Arguments:
#       frm                 {int}       frame to be evaluated 
#       species             {list}      list of lipid species in system
#       findHeadsAndTails   {list}      list of first and last beads in lipid tail
#       window              {int}       Window size of tilted lipids 
#       poreSort            {str}       Lipids Excluded From analysis
#
# Results:
#       
# 
#
# Necessary Revisions:
#       Needs to be revised to just take the 'top' bead in a lipid
#       rather than hard-code PO4

proc assignLeaflet {frm species findHeadsAndTails window poreSort} {
    set starts [lindex $findHeadsAndTails 0]
    set ends [lindex $findHeadsAndTails 1]

    ;# does leaflet check for different lipid species separately
    ;# because bead names may conflict between species
    for {set i 0} {$i < [llength $species]} {incr i} {
        set lipidType [lindex $species $i]
        set totalSel [atomselect top "resname $lipidType" frame $frm]
        set endNames [lindex $ends $i]

        ;# how many beads are in the given lipid species?
        set speciesBeadNum [llength [lsort -unique [$totalSel get name]]]
        
        ;# how many tails are in this lipid species?
        set numTails [llength $endNames]

        set startSel [atomselect top "resname $lipidType and name PO4" frame $frm]
        set endSel [atomselect top "resname $lipidType and name $endNames" frame $frm]
        set startZ [$startSel get z]
        set endZ [$endSel get z]
        $startSel delete
        $endSel delete
        
        set userList []
        set counter 0

        ;# iterate through each lipid in the system and calc average height of the endbeads
        for {set j 0} {$j < [llength $endZ]} {set j [expr $j+$numTails]} {
            set avgEndHeight [vecexpr [lrange $endZ $j [expr $j+$numTails-1]] mean]

            ;# subtract $avgendheight from the PO4 bead's height
            set avgHeight [expr [lindex $startZ $counter]-$avgEndHeight]

            ;# assign user value accordingly
            if {$avgHeight > $window} {
                lappend userList [lrepeat $speciesBeadNum 1.0]
            } elseif {$avgHeight < -$window} {
                lappend userList [lrepeat $speciesBeadNum 2.0]
            } else {
                lappend userList [lrepeat $speciesBeadNum 3.0]
            }
            incr counter
        }

        ;# convert list of lists into one long list
        set userVals [concatenateList $userList "NULL"]
        
        $totalSel set user $userVals
        $totalSel delete
    }

    if {$poreSort ne "NULL"} {
        ;# custom pore sorting proc for 5x29 and 7k3g
        pore_sorter_custom $frm $species $poreSort
    }
}

# printBinInfo (Previously: print_line_init)--
#
#       starts a new line in the print file that has the min/max 
#       r or x value for the bin, depending on if polar or cartesian
#
# Arguments:
#       file        {str}       File name
#       number      {int}       Bin values ranging from 0 to n bins
#       d1          {int}       Width of bins 
#       min         {int}       Minimum values for analysis
#
# Results:
#       prints line containing the min and max bin values to a file 
#
# Necessary Revisions/Problems:
#       min is needs to be removed

proc printBinInfo {file number d1 min} {
    puts -nonewline $file "[format {%0.2f} [expr $number*$d1+$min]]  [format {%0.2f} [expr ($number+1)*$d1+$min]]  "
}

# printValue (Previously: print_value)--
#
#       adds a value to be printed to a file
#
# Arguments:
#       file        {str}       File name
#       value       {}          Value to be printed
#       endLine     {int}       Determines if value should be added to same line or next line
#                               by 0 or 1
# Results:
#       prints line in file with value specified

proc printValue {file value endLine} {
    if {$endLine == 0} {
        puts -nonewline $file " $value" 
    } elseif {$endLine == 1} {
        puts $file " $value"
    } else {
        puts "Something went wrong - end_line should have value of 0 or 1 only"
        break
    }
}

# printFrame (Previously: print_frame)--
#
#       print an entire 2D array (usually 1 frame) to file
# 
# Arguments:
#       N1
#       outFiles
#       key
#       d1
#       min
#       N2
#       polar
#       selex
#
# Results:
#

proc printFrame {N1 outFiles key d1 min N2 polar selex} {

    set file [dict get $outFiles $selex $key fname]

    ;# starts new line in outfile with bin values
    for {set m 0.0} {$m < $N1} {set m [expr $m+1.0]} {
        set binstart [format {%0.2f} [expr $m*$d1+$min]]
        set binend [format {%0.2f} [expr ($m+1)*$d1+$min]]
        puts -nonewline $file "$binstart  $binend  "
        ;# prints bin values through ultimate value in one line
        for {set n 0.0} {$n < $N2} {set n [expr $n+1.0]} {
            if {[dict exists $outFiles $selex $key bin "$m,$n"]} {
                puts -nonewline $file " [dict get $outfiles $selex $key bin "$m,$n"]"
            } else {
                puts -nonewline $file " nan"
            }
        }
        ;# starts a new line
        puts $file " "
    }
}

# analyzeTails (Previously: tail_analyzer)--
#
#
# Arguments:
#       species     {list} 
#
# Results: 
#       returns a nested list: 
#       top level is by species, mid level is by tail in species
#       bottom level is by beads in tail
#
# e.g. a membrane with DO and DP lipids would be: 
# |-------------------------------------taillist------------------------------------|
#   |------------------DO-----------------| |------------------DP-----------------|
#     |-----tail0-----| |-----tail1-----|     |-----tail0-----| |-----tail1-----|
#
# { { {C1A C2A C3A C4A} {C1B C2B C3B C4B} } { {C1A D2A C3A C4A} {C1B D2B C3B C4B} } } 

proc analyzeTails { species } {
    set taillist []
    set letters "A B C D E F G H I J"
    foreach lipidtype $species {
        set tails []
        set sel [atomselect top "resname $lipidtype"]
        set res [$sel get resid]
        $sel delete
        set sel [atomselect top "resname $lipidtype and resid [lindex $res 0]"]
        set names [$sel get name]
        $sel delete
        foreach letter $letters {
            set tail []
            foreach nm $names {
                if {[string match ??${letter} $nm]} {
                    lappend tail $nm
                }
            }
            if {[llength $tail] != 0} {
                lappend tails $tail
            }
        }
        if {[llength $tails] != 0} {
            lappend taillist $tails
        }
    }

    ;# returns top/bottom beads in lipid tails for leaflet sorting
    set heads_and_tails [findHeadsAndTails $species $taillist]

    ;# one list with all the bead names for convenience
    set full_tails []
    foreach lipidtype $taillist {
        foreach tail $lipidtype {
            foreach bead $tail {
                lappend full_tails $bead
            }
        }
    }

    return [list $taillist $heads_and_tails $full_tails]
}

# numberTails (Previously: Tail_numberer)--
#
# Arguments:
#       species
#       tailList
#
# Results:
#       
;# tail_numberer changes user3 to hold a tail number.
;# This makes different tails easily separable for tilt/order analysis
proc numberTails { species tailList } {
    for {set lipidtype 0} {$lipidtype < [llength $species]} {incr lipidtype} {
        for {set tail 0} {$tail < [llength [lindex $tailList $lipidtype]]} {incr tail} {
            set sel [atomselect top "resname [lindex $species $lipidtype] and name [lindex [lindex $tailList $lipidtype] $tail]"]
            for {set frm 0} {$frm < [molinfo top get numframes]} {incr frm} {
                $sel frame $frm 
                $sel update
                $sel set user3 [expr $tail+1]
            }
            $sel delete
        }
    }
}

# assignBins (Previously: bin_assigner)--
#       makes a list of bins based on x, y values and the coordinate system
# Arguments:
#       xVals       {list}
#       yVals       {list}
#       binWidth1   {float}
#       binWidth2   {float}
#       thetaDeg    {float}
#       polar       {int}
#       frm         {int}  
#
# Results:
#       returns two lists of bins in the x or y direction

proc assignBins {xVals yVals binWidth1 binWidth2 thetaDeg polar frm} {
    
    if {$polar == 1} {
        ;# use polar (r,theta) bins

        ;#calculate r: distance from origin for all x,y pairs
        set r_vals [vecexpr [vecexpr [vecexpr $xVals sq] [vecexpr $yVals sq] add] sqrt]
        
        ;#turn into bin numbers rather than r values
        set dim1_bins [vecexpr [vecexpr $r_vals $d1 div] floor]
        
        ;#calculate theta: use atan2 to get values for al x,y pairs
        set theta_vals [vecexpr $yVals $xVals atan2 pi div 180 mult]

        ;#atan2 gives values from -180 to 180; shifting to 0 to 360
        for {set i 0} {$i<[llength $theta_vals]} {incr i} {
            if {[lindex $theta_vals $i] < 0} {
                set theta_vals [lreplace $theta_vals $i $i [expr [lindex $theta_vals $i]+360]]
            }
        }

        ;#turn into bin numbers rather than theta values
        set dim2_bins [vecexpr [vecexpr $theta_vals $thetaDeg div] floor]
        
    } elseif {$polar == 0} {
        ;# use cartesian (x,y) bins
        
        ;# shift all values so that values within unitcell are greater than 0 and less than unitcell len
        set xlen [molinfo top get a frame $frm]
        set ylen [molinfo top get b frame $frm]
        set xmin [expr -$xlen/2.0]
        set ymin [expr -$ylen/2.0]
        set xVals [vecexpr $xVals $xmin sub]
        set yVals [vecexpr $yVals $ymin sub]

        ;# any negative values or values exceeding unitcell len are lipids that flipped across PBC
        ;# and should be put back for binning purposes (but not for order params purposes!)
        for {set i 0} {$i < [llength $xVals]} {incr i} {
            if {[lindex $xVals $i] < 0} {
                lset xVals $i [expr $xlen-1.0]
            } elseif {[lindex $xVals $i] > $xlen} {
                lset xVals $i 0
            }
            if {[lindex $yVals $i] < 0} {
                lset yVals $i [expr $ylen-1.0]
            } elseif {[lindex $yVals $i] > $ylen} {
                lset yVals $i 0
            }
        }

        ;# turn into bin numbers rather than x,y values
        set dim1_bins [vecexpr [vecexpr $xVals $binWidth1 div] floor]
        set dim2_bins [vecexpr [vecexpr $yVals $binWidth2 div] floor]
    }

    return [list $dim1_bins $dim2_bins]
}


# createOutfiles (Previously: create_outfiles) --
#
# Arguments:
#       system
#       quanity
#       headNames
#       TailList
#       coordSystem
#       folderName
#
# Results:
#
;# create a dict containing all the outfile names/addresses
;# density segregates by species
;# tilt and order segregate by species and tail number
proc createOutfiles {system quantity_of_interest headnames species taillist coordsys foldername} {
    file mkdir "${foldername}/tcl_output"
    if {$quantity_of_interest eq "height_density"} {
        dict set outfiles z1z2 heights_up fname [open "${foldername}/tcl_output/${system}.zone.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles z1z2 heights_down fname [open "${foldername}/tcl_output/${system}.ztwo.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles z0 heights_zzero fname [open "${foldername}/tcl_output/${system}.zzero.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles z1z2 counts_up fname [open "${foldername}/tcl_output/${system}.zone.${headnames}.${coordsys}.totdensity.dat" w]
        dict set outfiles z1z2 counts_down fname [open "${foldername}/tcl_output/${system}.ztwo.${headnames}.${coordsys}.totdensity.dat" w]
        dict set outfiles z0 counts_zzero fname [open "${foldername}/tcl_output/${system}.zzero.${headnames}.${coordsys}.totdensity.dat" w]
        foreach lipidtype $species {
            dict set outfiles z1z2 density_up_${lipidtype} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.zone.${coordsys}.density.dat" w]
            dict set outfiles z1z2 density_down_${lipidtype} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.ztwo.${coordsys}.density.dat" w]
            dict set outfiles z0 density_zzero_${lipidtype} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.zzero.${coordsys}.density.dat" w]
        }
    } elseif {$quantity_of_interest eq "tilt_order"} {
        for {set i 0} {$i < [llength $taillist]} {incr i} {
            set lipidtype [lindex $species $i]
            for {set j 0} {$j < [llength [lindex $taillist $i]]} {incr j} {
                set tailnum "tail$j"
                set taillength [llength [lindex [lindex $taillist $i] $j]]
                dict set outfiles $taillength tilts_up_${lipidtype}_${tailnum} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.${tailnum}.zone.${coordsys}.tilt.dat" w]
                dict set outfiles $taillength tilts_down_${lipidtype}_${tailnum} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.${tailnum}.ztwo.${coordsys}.tilt.dat" w]
                dict set outfiles $taillength order_up_${lipidtype}_${tailnum} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.${tailnum}.zone.${coordsys}.order.dat" w]
                dict set outfiles $taillength order_down_${lipidtype}_${tailnum} fname [open "${foldername}/tcl_output/${system}.${lipidtype}.${tailnum}.ztwo.${coordsys}.order.dat" w]
            }
        }
    }
    
    return $outfiles
}

# measureBoxSize (Prevoiusly: measure_box_size)--
#       a proc for measuring box size
# Arguments:
#       frame
#       dimension
# Results:
#   
# Necessary Revisions/Problems: 
#       [molinfo top get a] is unreliable for certain coarse-grain simulations 
proc measureBoxSize {frame dim} {
    if {($dim ne "x") && ($dim ne "y") && ($dim ne "z")} {
        puts "dim must be x, y, or z"
        return
    }
    set sel [atomselect top all frame $frame]
    set dimvals [$sel get $dim]
    set max [lindex [lsort -real $dimvals] end]
    set min [lindex [lsort -real $dimvals] 0]
    set len [expr $max-$min]

    return $len
}

# prepareBins (Previously: bin_prep)--
#       Determines the number of bins and the step length in each dimension
# Arguments:
#       frameNumber
#       polar
#       min
#       drN1
#       N2
#
# Results:
#       
proc prepareBins {nframes polar min dr_N1 N2} {
    
    if {$polar == 1} {
    
        dict set bindims d1 $dr_N1

        #measure box size at final frame to get bin values
        set box_x [molinfo top get a frame [expr $nframes-1]]

        set box_r [expr int($box_x)/2]
        set rrange [expr $box_r-$min]
        
        #calculate number of dim1 bins from d1 and range1
        if {[expr $rrange%$dr_N1] == 0} { 
            dict set bindims N1 [expr [expr $rrange/$dr_N1]-1] 
        } else {
            dict set bindims N1 [expr $rrange/$dr_N1]
        }

        #calculate dtheta in degrees and radians
        dict set bindims dthetadeg [expr 360/[expr $N2*1.0]]
        global M_PI
        dict set bindims d2 [expr 2*$M_PI/$N2]
        dict set bindims N2 $N2
    
    } elseif {$polar == 0} {

        dict set bindims N1 $dr_N1
        ;# fix this if you ever want to implement rectangular systems
        dict set bindims N2 $dr_N1
        
        set bindims [updateDimensions $bindims 0]

        dict set bindims dthetadeg "NULL"
    }

    return $bindims
}

# calculateReferenceHeight (Previously: calc_ref_height)--
#
# Arguments:
#       configDict
#       frame
#
# Result:
#   
proc calculateReferenceHeight {config_dict frm} {
    if {[dict get $config_dict reference_point] ne "NULL"} {
        set ref_bead [atomselect top [dict get $config_dict reference_point] frame $frm]
        set ref_height [$ref_bead get z]
        $ref_bead delete
        set ref_height [vecexpr $ref_height mean]
    } else {
        set ref_height "NULL"
    }
    return $ref_height
}

# getSelectionInformation (Previously: grab_sel_info)--
#
# Arguments:
#       sel
#       referenceHeight
#
# Results:
#       
proc grab_sel_info {sel ref_height} {
    dict set sel_info xvals_list [$sel get x]
    dict set sel_info yvals_list [$sel get y]
    dict set sel_info resid_list [$sel get resid]
    dict set sel_info lipid_list [$sel get resname]
    dict set sel_info name_list [$sel get name]

    ;# warn the user if a selection is empty
    if {[llength [dict get $sel_info xvals_list]] == 0} {
        puts "selection $sel has no atoms in it"
    }

    ;# the z vals are subtracted by a reference height provided in cell_prep 
    if {$ref_height ne "NULL"} {
        dict set sel_info zvals_list [vecexpr [$sel get z] $ref_height sub]
    } else {
        dict set sel_info zvals_list [$sel get z]
    }   

    ;# user contains a 1 or 2 for outer or inner leaflet, respectively
    dict set sel_info leaflet_list [$sel get user]

    ;# user3 contains an int that describes which tail in the lipid this is
    ;# E.G. POPC will have 0 or 1 (it has two tails)
    ;# E.G. OANT will have 0, 1, 2, 3, 4, or 5 (it has 6 tails)
    set tail_list [$sel get user3]
    dict set sel_info tail_list [vecexpr $tail_list 1 sub]

    return $sel_info
}

# updateDimensions (Previously: update_dims)--
#
# Arguments: 
#       binDimension
#       frame
#
# Results:
#   
proc updateDimensions {bindims frm} {
    set x [molinfo top get a frame $frm]
    set y [molinfo top get b frame $frm]

    dict set bindims d1 [expr $x/[expr [dict get $bindims N1]*1.0]]
    dict set bindims d2 [expr $y/[expr [dict get $bindims N2]*1.0]]

    return $bindims
}

# concatenateNames (Previously: concat_names) --
#
# Arguments: 
#       headnames
#
# Results:
#   
;# concatenate all beadnames together for file naming purposes
proc concat_names { headnames } {
    if {[llength $headnames] > 1} {
        set condensed_name [lindex $headnames 0]
        for {set i 1} {$i < [llength $headnames]} {incr i} {
            set addname [lindex $headnames $i]
            set condensed_name "$condensed_name.$addname"
        }
    } elseif {[llength $headnames] == 1} {
        set condensed_name [lindex $headnames 0]
    } else {
        puts "headnames must contain a bead name"
        break
    }
    return $condensed_name
}

# createResidueDictionaries (Previously: create_res_dict)--
#
# Arguments:
#       species
#       headNames
#       lipidList
#       nameList
#       residList
#       dimensionOneBinList
#       dimensionTwoBinList
#       leafletList
#       selex
#
# Results:
#  
;# binning is controlled by the location of beads named as $headnames in cell_prep
;# creates a nested dict with keys set to the name "bin1#,bin2#,leaflet#" and
;# values corresponding to the indices of lipids in that bin
proc createResidueDictionaries { species headnames lipid_list name_list resid_list dim1_bins_list dim2_bins_list leaflet_list selex} {
    ;# initialize a nested dict with a dummy key and value
    dict set res_dict dummy "dummy"
    
    if {$selex ne "z0"} {
        for {set i 0} {$i < [llength $lipid_list]} {incr i} {
            if {([lindex $leaflet_list $i] == 3) || ([lindex $leaflet_list $i] == 4)} {
                continue
            } elseif {([lsearch $species [lindex $lipid_list $i]] != -1) && ([lsearch $headnames [lindex $name_list $i]] != -1)} {
                set bin "[lindex $dim1_bins_list $i],[lindex $dim2_bins_list $i]"
                set bin_leaf "$bin,[expr int([lindex $leaflet_list $i])]"
                if {[dict exists $res_dict $bin_leaf]} {
                    dict append res_dict $bin_leaf " $i"
                } else {
                    dict set res_dict $bin_leaf $i                
                }
            }
        }
    } else {
        ;# zzero needs to be handled separately; keys are set to "bin1#,bin2#,3"
        for {set i 0} {$i < [llength $lipid_list]} {incr i} {
            set bin "[lindex $dim1_bins_list $i],[lindex $dim2_bins_list $i]"
            set bin_leaf "$bin,3"
            if {[dict exists $res_dict $bin_leaf]} {
                dict append res_dict $bin_leaf " $i"
            } else {
                dict set res_dict $bin_leaf $i                
            }
        }
    }

    ;# delete the dummy
    dict unset res_dict dummy
    
    return $res_dict
}

# outputDensityNormalizationInformation (Previously: output_density_norm_info)
#
#       Calculates the normalization factor for 
#       density enrichment calculations
#
# Arguments:
#       start
#       frameNumber
#       step
#       species
#       system
#       headNames
#       coordSystem
#       folderName
#
# Results:
#       
;# calculates the normalization factor for density enrichment calculations
proc outputDensityNormalizationInformation {start nframes step species system headnames coordsys foldername} {
    set arealist []
    for {set frm $start} {$frm <= $nframes} {set frm [expr $frm+$step]} {
        lappend arealist [expr [molinfo top get a frame $frm]*[molinfo top get b frame $frm]]
    }
    set avgarea [vecexpr $arealist mean]
    set normfactor_outfile [open "${foldername}/tcl_output/${system}.${coordsys}.density.normfactor.dat" w]
    foreach spec $species {
        set sel [atomselect top "resname $spec"]
        set names [lsort -unique [$sel get name]]
        set Sb 0
        foreach name $names {
            if {[lsearch $headnames $name] != -1} {
                incr Sb
            }
        }
        set Nb [llength [lsort -unique [$sel get resid]]]
        $sel delete
        set normfactor [expr $avgarea/[expr $Nb*$Sb/2.0]]
        puts $normfactor_outfile "$spec $normfactor"
    }
    
    close $normfactor_outfile
}

# sortTailLength (Previously: tail_length_sorter) --
#
# Arguments:
#       species
#       acylNames
#
# Results:
#   
;# nougat tilt and order calculations segregate by lipid species and tail.
;# Several different tails in a system may be the same length as each other,
;# and so their calculations can be combined in the same loop iteration.
;# This proc determines the unique tail lengths in the system and creates a
;# list of atomselection texts that correspond.
proc sortTailLength {species acyl_names} {
    set lenlist []
    for {set i 0} {$i < [llength $acyl_names]} {incr i} {
        foreach tail [lindex $acyl_names $i] {
            lappend lenlist [llength $tail]
        }
    }
    set lengthlist [lsort -unique $lenlist]
    set sellist []
    foreach length $lengthlist {
        set resnamelist []
        set namelist []
        for {set i 0} {$i < [llength $acyl_names]} {incr i} {
            foreach tail [lindex $acyl_names $i] {
                if {[llength $tail] == $length} {
                    lappend resnamelist [lindex $species $i]
                    foreach nm $tail {
                        lappend namelist $nm 
                    }
                }
            }
        }
        lappend sellist "resname [lsort -unique $resnamelist] and name [lsort -unique $namelist]"
    }
    return [list $sellist $lengthlist]
}

# getCosineTheta (Previously: get_costheta)--
#
# Arguments:
#       start
#       end
#
# Results:
#   
;# cosine theta is the dot product of n_{1,2} and the vector [0 0 1]
;# this corresponds to the 3rd value of n_{1,2} 
proc get_costheta {start end} {   
    set r12 [vecsub $start $end]
    set n12 [vecnorm $r12]
    return [lindex $n12 2]
}

# calculateOrderParameters (Previously: order_params)--
# 
#       calculates order parameter of each lipid 
#       tail of a given length
# 
# Arguments:
#       length
#       xValues
#       yValues
#       zValues
#
# Results:
#       
proc calculateOrderParameters {length xvals yvals zvals} {
    set order_list []
    set temp_list []
    for {set i 1} {$i <= [llength $xvals]} {incr i} {
        if {[expr $i%$length] == 0} {
            ;# when this is TRUE, you've gotten cos2theta for each of the bonds in your tail
            ;# already and now you need to average them
            set avg [vecexpr $temp_list mean]
            set order [expr {$avg * 1.5 - 0.5}]
            lappend order_list [lrepeat $length $order]
            set temp_list []
        } else {
            ;# calculate cos2theta for each bond in the tail and append to a temp list
            set start [list [lindex $xvals $i] [lindex $yvals $i] [lindex $zvals $i]]
            set end [list [lindex $xvals [expr $i-1]] [lindex $yvals [expr $i-1]] [lindex $zvals [expr $i-1]]]
            set costheta [get_costheta $start $end]
            lappend temp_list [expr {$costheta * $costheta}]
        }
    }

    ;# this is a list of lists, but we want just a list
    set final_order_list [concatenateList $order_list "NULL"]
    
    return [list $final_order_list]
}

# averageTiltAndOrderParameter (Previously: tilt_order_averaging)
#
#       Uses res_dict entries to compute bin averages, then assigns them to the correct outfile
#
# Arguments:
#       residueDictionary
#       outfiles
#       leafletList
#       lipidList
#       tilts
#       orders
#       tailList
#       selex
#
# Results:
#       
proc averageTiltAndOrderParameter {res_dict outfiles leaflet_list lipid_list tilts orders tail_list selex} {
    dict set counts placeholder "dummy"
    dict for {bin indices} $res_dict {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr {[string length $bin] - 3}]]
        foreach indx $indices {
            set tailnum [expr int([lindex $tail_list $indx])]
            set species [lindex $lipid_list $indx]
            if {$leaf == 1} {
                set tilt_key "tilts_up_${species}_tail${tailnum}"
                set order_key "order_up_${species}_tail${tailnum}"
            } elseif {$leaf == 2} {
                set tilt_key "tilts_down_${species}_tail${tailnum}"
                set order_key "order_down_${species}_tail${tailnum}"
            } else {
                puts "Something has gone wrong with the binning"
                return
            }
            if {[dict exists $counts $selex $order_key bin $correct_bin]} {
                set oldcount [dict get $counts $selex $order_key bin $correct_bin]
                set newcount [expr $oldcount+1.0]
                set oldorder [dict get $outfiles $selex $order_key bin $correct_bin]
                set newordersum [expr {$oldorder * $oldcount + [lindex [lindex $orders 0] $indx]}]
                set neworder [expr $newordersum/$newcount]
                set oldtilt [dict get $outfiles $selex $tilt_key bin $correct_bin]
                set newtiltsum [vecexpr [vecexpr $oldtilt $oldcount mult] [lindex $tilts $indx] add]
                set newtilt [vecexpr $newtiltsum $newcount div]
                dict set outfiles $selex $order_key bin $correct_bin $neworder
                dict set counts $selex $order_key bin $correct_bin $newcount
                dict set outfiles $selex $tilt_key bin $correct_bin $newtilt
            } else {
                dict set outfiles $selex $order_key bin $correct_bin [lindex [lindex $orders 0] $indx]
                dict set outfiles $selex $tilt_key bin $correct_bin [lindex $tilts $indx]
                dict set counts $selex $order_key bin $correct_bin 1.0
            }
        }
    }
    dict unset counts $selex
    dict unset counts placeholder
    return $outfiles
}

# averageHeightAndDensity (Previously: height_density_averaging)
#
#       uses res_dict entries to compute bin averages, then assigns them to the correct outfile
#
# Arguments:
#       residueDictionary
#       outfiles
#       leafletList
#       lipidList
#       zValsList
#       NameList
#
# Results:
#   
proc averageHeightAndDensity {res_dict outfiles leaflet_list lipid_list zvals_list name_list} {
    dict for {bin indices} $res_dict {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr {[string length $bin] - 3}]]
        foreach indx $indices {
            set species [lindex $lipid_list $indx]
            if {$leaf == 1} {
                set field_key "z1z2"
                set dens_key "density_up_${species}"
                set height_key "heights_up"
                set counts_key "counts_up"
            } elseif {$leaf == 2} {
                set field_key "z1z2"
                set dens_key "density_down_${species}"
                set height_key "heights_down"
                set counts_key "counts_down"
            } elseif {$leaf == 3} {
                set field_key "z0"
                set dens_key "density_zzero_${species}"
                set height_key "heights_zzero"
                set counts_key "counts_zzero"
            } else {
                puts "something has gone wrong with the binning"
                return
            }
            if {[dict exists $outfiles $field_key $height_key bin $correct_bin]} {
                set oldcount [dict get $outfiles $field_key $counts_key bin $correct_bin]
                set newcount [expr $oldcount+1.0]
                set oldavg [dict get $outfiles $field_key $height_key bin $correct_bin]
                set newsum [expr {$oldavg * $oldcount + [lindex $zvals_list $indx]}]
                set newavg [expr $newsum/$newcount]
                dict set outfiles $field_key $height_key bin $correct_bin $newavg
                dict set outfiles $field_key $counts_key bin $correct_bin $newcount
                if {[dict exists $outfiles $field_key $dens_key bin $correct_bin]} {
                    set olddens [dict get $outfiles $field_key $dens_key bin $correct_bin]
                    dict set outfiles $field_key $dens_key bin $correct_bin [expr $olddens+1.0]
                } else {
                    dict set outfiles $field_key $dens_key bin $correct_bin 1.0
                }       
            } else {
                dict set outfiles $field_key $height_key bin $correct_bin [lindex $zvals_list $indx]
                dict set outfiles $field_key $counts_key bin $correct_bin 1.0
                dict set outfiles $field_key $dens_key bin $correct_bin 1.0
            }
        }
    }
    return $outfiles
}

# setBetaValues (Previously: set_beta_vals)--
#       
# Arguments:
#       inclusionSelection
#       species
#
# Results:
#       
proc setBetaValues {inclusion_sel species} {
    puts "starting to fill beta values"
    
    if {$inclusion_sel ne "NULL"} {
        set inclusion [atomselect top $inclusion_sel]
        $inclusion set beta 0
        $inclusion delete 

        set excl_sel [atomselect top "not $inclusion_sel and not resname W"]
    } else {
        set excl_sel [atomselect top "not resname W"]
    }
    
    set other_resnames [lsort -unique [$excl_sel get resname]]
    $excl_sel delete

    set counter 1

    foreach resnm $other_resnames {
        if {[string first - $resnm] != -1} {
            set resnm \'$resnm\'
        } elseif {[string first + $resnm] != -1} {
            set resnm \'$resnm\'
        }
        set res_sel [atomselect top "resname $resnm"]
        set resids [lsort -unique [$res_sel get resid]]
        
        set sel [atomselect top "resname $resnm and resid [lindex $resids 0]"]
        set num [$sel num]
        $sel delete

        set counterlist []

        for {set i 0} {$i < [llength $resids]} {incr i} {
            for {set j 0} {$j < $num} {incr j} {
                lappend counterlist $counter
            }
            incr counter 
        }

        $res_sel set beta $counterlist

        $res_sel delete
    }

    set water_sel [atomselect top "resname W"]
    set indxs [$water_sel get index]
    $water_sel set beta $indxs
    $water_sel delete

    puts "beta values complete"
    return
}

# readPolar (Previously: read_polar) --
#
# Arguments:
#       polar
#
# Results:
#       
proc readPolar {polar} {
    ;# generate string for polar or cartesian coordinates
    if {$polar == 1} {
        set coordsys "polar"
    } elseif {$polar == 0} {
        set coordsys "cart"
    } else {
        puts "polar must be 1 or 0"
        break
    }
    return $coordsys
}

# createAtomSelections (Previously: create_atomselections)--
#       
# Arguments:
#       quanity
#       configDictionary
#
# Results:
#       
proc createAtomSelections {quantity_of_interest config_dict} {
    ;#atomselections setup as dict
    if {$quantity_of_interest eq "height_density"} {
        dict set selections z1z2 [atomselect top "resname [dict get $config_dict species] and name [dict get $config_dict full_tails]"]
        dict set selections z0 [atomselect top "resname [dict get $config_dict species] and ((user 1 and within 6 of user 2) or (user 2 and within 6 of user 1))"]

    } elseif {$quantity_of_interest eq "tilt_order"} {
        set lists [sortTailLength [dict get $config_dict species] [dict get $config_dict acyl_names]]
        set sellist [lindex $lists 0]
        set lenlist [lindex $lists 1]
        foreach sel $sellist len $lenlist {
            dict set selections $len [atomselect top "$sel"] 
        }
    }

    foreach key [dict keys $selections] {
        set sel [dict get $selections $key]
        $sel uplevel 1
    }

    return $selections
}

;#********************************;#
;# Liam scripts or custom scripts ;#
;#********************************;#

;# Alignment based off vmd alignment
proc Align { stuff } {
    puts "Align start"
    set nframes [molinfo top get numframes]
    set ref [atomselect top $stuff frame 0]
    for {set frames 1} {$frames < $nframes} {incr frames} {
        set com [atomselect top $stuff frame $frames]
        set TM [measure fit $com $ref]
        $com delete
        set move_sel [atomselect top "all" frame $frames]
        $move_sel move $TM
        $move_sel delete
    }
    $ref delete
    puts "Align end"
}


proc separate_chains {molid cutoff} {
    ;# Multichain proteins default to chain X in martinize.py
    ;# separate_chains will jump from atom to atom and anywhere the 
    ;# distance is greater than $cutoff it will assume there is a new subunit.

    set sel [atomselect $molid "name BB SC1 to SC4"]
    set indxlist [$sel get index]
    set indmin [::tcl::mathfunc::min {*}$indxlist]
    set indmax [::tcl::mathfunc::max {*}$indxlist]
    $sel delete

    set list1 [list $indmin]

    for {set i $indmin} {$i < $indmax} {incr i} {
        set sel1 [atomselect $molid "index $i"]
        set sel2 [atomselect $molid "residue [expr $i+1]"]

        set loc1 [lindex [$sel1 get {x y z}] 0]
        set loc2 [lindex [$sel2 get {x y z}] 0]

        set dist [vecdist $loc1 $loc2]

            if {$dist > $cutoff} {
                lappend list1 $i
                lappend list1 [expr $i+1]
            }
        $sel1 delete
        $sel2 delete
    }

    set chars [list A B C D E]
    set k 0
    lappend list1 $indmax

    foreach {i j} $list1 {
        set sel [atomselect $molid "index $i to $j"]
        $sel set chain [lindex $chars $k]
        $sel delete
        set k [expr $k+1]
    }
}

;# THIS IS FOR 5X29!
proc set_occupancy {molid} {

    set sel [atomselect $molid "resid 0 to 39"]
    $sel set occupancy 1
    $sel delete
    set sel [atomselect $molid "resid 40 to 51"]
    $sel set occupancy 2
    $sel delete
    set sel [atomselect $molid "resid 52 to 65"]
    $sel set occupancy 3
    $sel delete

}

;# gets theta value from x and y pair
;# this is only used by Protein_Position
proc get_theta {x y} {
    set tmp  [expr {atan2($y,$x)}]
    global M_PI

    ;# atan2 returns in a range of -pi to pi
    ;# shift it to 0 to 2pi
    if {$tmp < 0} {
        set theta [expr {2*$M_PI + $tmp}]    
    } else {
        set theta $tmp
    }

    ;# change to degrees
    return [RtoD $theta]
}

;# Ouputs position of the centered protein in a membrane
;# accross both leaflets
;# only useful for analysis later - doesn't impact your polar density script
proc Protein_Position {name hnames tnames} {
    ;# in order to use this, must have your TMD chains separated and saved as occupancy 3
    set chain_names [list "A" "B" "C" "D" "E"]

    set lastframe [expr [molinfo top get numframes]-1]

    set zone_sel [atomselect top "(name $hnames and chain U) and within 6 of name BB"]
    set zone_zvals [$zone_sel get z]
    set zone_Ht [vecexpr $zone_zvals mean]
    $zone_sel delete

    set ztwo_sel [atomselect top "(name $hnames and chain L) and within 6 of name BB"]
    set ztwo_zvals [$ztwo_sel get z]
    set ztwo_Ht [vecexpr $ztwo_zvals mean]
    $ztwo_sel delete

    set zmid_sel [atomselect top "name $tnames and within 6 of name BB"]
    set zmid_zvals [$zmid_sel get z]
    set zmid_Ht [vecexpr $zmid_zvals mean]
    $zmid_sel delete

    foreach ht [list $zone_Ht $ztwo_Ht $zmid_Ht $zmid_Ht] eqtxt [list "zone" "ztwo" "zzero" "zplus"] {
    set fout [open "${name}_helcoords_${eqtxt}.dat" w]
        puts $fout  "#These are the positions of your TMD helices in polar coords"
        foreach chnm $chain_names {
                set sel [atomselect top "(chain ${chnm} and name BB and occupancy 1) and (z < [expr $ht+5] and z > [expr $ht-5])" frame $lastframe]
                set com [measure center $sel weight mass]
                $sel delete
                set x [lindex $com 0]
                set y [lindex $com 1]
                set r [expr sqrt($x*$x+$y*$y)]
                set theta [get_theta $x $y]
                puts "chain ${chnm} and $r $theta"

                puts -nonewline $fout "$r $theta "
            }
            puts $fout ""
        close $fout
    }
}

proc Center_System {wrap_sel species inclusion_sel} {
    puts "${wrap_sel}"
    puts "Center_System now running"

    setBetaValues $inclusion_sel $species
    qunwrap compound beta
    if {$inclusion_sel ne "NULL"} {
        qwrap compound beta center $inclusion_sel 
    }
    qwrap compound beta center $wrap_sel

    puts "Center_System finished!"
}

;# removes lipids from analysis (by setting user to 4)
;# current removal criteria:
;# within 5 angstroms of the pore center
;# within 30 angstroms of protein BB resid 30 (upper pore-lining residue)
proc pore_sorter_custom {frm species inc} {
    
    ;# define pore center
    set sel [atomselect top "name BB and resid 30" frame $frm]
    set com [measure center $sel]
    set x [lindex $com 0]
    set y [lindex $com 1]
    $sel delete

    ;# "same resid as" doesn't work and "same residue as" has no meaning in CG trajs,
    ;# so have to do this in two steps
    if {$inc eq "5x29"} {
        set porebeads [atomselect top "(resname $species and within 12 of (name BB and resid 30)) or (resname $species and ((x-$x)*(x-$x)+(y-$y)*(y-$y) <= 25))" frame $frm]
    } elseif {$inc eq "7k3g"} {
        set porebeads [atomselect top "resname $species and ((x-$x)*(x-$x)+(y-$y)*(y-$y) <= 25)" frame $frm]
    }
    set badresids [lsort -unique [$porebeads get resid]]
    $porebeads delete
    if {[llength $badresids] != 0} {
        set porelipids [atomselect top "resname $species and resid $badresids" frame $frm]
        $porelipids set user 4.0
        $porelipids delete
    } 
}
