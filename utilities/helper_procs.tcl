# convertRadianToDegree (Previously: RtoD)--
#
#       Converts radians to degrees
#
# Arguments:
#       rad         {float}    Angle in radians
#
# Results:
#       Returns angle in degrees

proc convertRadianToDegree {rad} {
    global M_PI
    return [expr $rad*180.0/$M_PI]
}

# fitVecsToSel (Previously: tilt_angles)--
#
#       Performs a least squares fit for each lipid tail of N bead length
#       given the x, y and z values in order. This allows you to use 1 instance 
#       of $lsqnormfactor only,rather than calculate it on the fly each time.
#
# Arguments:
#       length      {int}     Number of beads in a single lipid tail.
#       xVals       {list}    Ordered list of x coordinates. 
#       yVals       {list}    Ordered list of y coordinates. 
#       zVals       {list}    Ordered list of z coordinates.
#
# Results:
#       The result is a list of list that contains a normalized fitted vector 
#       and the number of beads in an acyl chain.   

proc fitVecsToSel {length xVals yVals zVals} {
    set tiltList []
    set lsqNormFactor [calculateLsqNormFactor $length]
    set xvec [fitTailVectors $length $xVals $lsqNormFactor]
    set yvec [fitTailVectors $length $yVals $lsqNormFactor]
    set zvec [fitTailVectors $length $zVals $lsqNormFactor]
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

    set lsqNormFactor [vecaddScalar $I [expr -1.0*$d]]

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
        if {[llength $lsqNormFactor] == 1} {
                lappend fitValues [vecsum [vecscale $coords $lsqNormFactor]]
            } else {
                lappend fitValues [vecsum [vecmul $coords $lsqNormFactor]]
            }
    }
    return $fitValues
}

# concatenateList (Previously: cat_list)--
#
#       Concatenates list elements into one long string, separated by 
#       spaces on either side of the $delimiter.
# 
# Arguments:
#       inputList  {List}  list of elements to concatenate
#       delimiter  {str}   string to be input between items in inputList
#
# Results:
#       Returns a list of elements with the specified delimiter between. "NULL" will 
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
#       tailList  {List}  nested list of tails organized by lipid type (see analyzeTails) 
#   
# Results:
#       Returns a nested lists containing the starting beads and 
#       ending beads for a given lipid's acyl chains.
#
#       Ex.
#       >>> findHeadsAndTails {{{C1A D2A C3A C4A} {C1B C2B C3B C4B}}}
#       {{C1A C1B}} {{C4A C4B}}
#       E.G. POPC start beads would be "C1A C1B" and end beads would be "C4A C4B"

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
#       axis     {str}   axis of rotation(x, y, z)
#       degree   {int}   degree of rotation from 0 to 360
#       start    {int}   start frame
#       stop     {int}   stop frame
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
#       poreSort            {str}       lipids Excluded From analysis
#
# Results:
#       
#       all lipids will contain a user value corresponding to upper or lower leaflet
#
# Necessary Revisions:
#       - Needs to be revised to just take the 'top' bead in a lipid
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
            set avgEndHeight [vecmean [lrange $endZ $j [expr $j+$numTails-1]]] 

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
#       file        {str}       file name
#       number      {int}       bin values ranging from 0 to n bins
#       d1          {int}       width of bins 
#       min         {int}       minimum values for analysis
#
# Results:
#       prints line containing the min and max bin values to a file 
#
# Necessary Revisions/Problems:
#       - min is only partially implemented

proc printBinInfo {file number d1 min} {
    puts -nonewline $file "[format {%0.2f} [expr $number*$d1+$min]]  [format {%0.2f} [expr ($number+1)*$d1+$min]]  "
}

# printValue (Previously: print_value)--
#
#       adds a value to be printed to a file
#
# Arguments:
#       file        {str}       File name
#       value       {float}     Value to be printed
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
#       N1          {float}     bin dimension in x direction/radial bin dimensions
#       outFiles    {dict}      dictionary of outputfiles and various membrane related quanities 
#       key         {str}       key for value associated with outFiles dictionary
#       d1          {float}     radial bin width
#       min         {int}       controls how much of the region around the protein is analyzed (not fully implemented)  
#       N2          {float}     bin dimension in y direction
#       polar       {int}       1 or 0 denoting if system is analyzed in polar 
#                               or cartesian coordinartes
#       selex       {str}       a dict key that holds an atomselection as its value
#
# Results:
#
#       Prints 2D arrays of information of specific frame to output files
#
#Necessary Revisions/Problems:
#       - min needs to be implemented 

proc printFrame {N1 outfiles key d1 min N2 polar selex} {
    set file [dict get $outfiles $selex $key fname]

    ;# starts new line in outfile with bin values
    for {set m 0.0} {$m < $N1} {set m [expr $m+1.0]} {
        set binstart [format {%0.2f} [expr $m*$d1+$min]]
        set binend [format {%0.2f} [expr ($m+1)*$d1+$min]]
        puts -nonewline $file "$binstart  $binend  "
        ;# prints bin values through ultimate value in one line
        for {set n 0.0} {$n < $N2} {set n [expr $n+1.0]} {
            if {[dict exists $outfiles $selex $key bin "$m,$n"]} {
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
#       Sorts tails of lipids by head group and tail type
#
# Arguments:
#       species     {list}      list of lipid species in system
#
# Results: 
#       returns a nested list: top level is by species, mid level 
#       is by tail in species bottom level is by beads in tail
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
    set findHeadsAndTails [findHeadsAndTails $taillist]

    ;# one list with all the bead names for convenience
    set full_tails []
    foreach lipidtype $taillist {
        foreach tail $lipidtype {
            foreach bead $tail {
                lappend full_tails $bead
            }
        }
    }
    return [list $taillist $findHeadsAndTails $full_tails]
}

# numberTails (Previously: Tail_numberer)--
#
#       Sorts lipid tails in user3 by an arbitrary tail number
#
# Arguments:
#       species     {list}      list of lipid species in system
#       tailList    {list}      nested list of tails organized by lipid type (see analyzeTails)
#
# Results:
#       tail list in user3 is assigned values based on number of tails in lipid

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


# vecFloor
#       Turns a list of floats into a list of ints, flooring all values in the process.
#
# Arguments:
#       inputList               {list}      a list of numbers
#
# Results:
#
#       Returns a list of (floored) ints.

proc vecFloor {inputList} {
    set outputList []
    foreach item $inputList {
        lappend outputList [expr {int($item)}]
    }
    return $outputList
}


# vecSqrt
#       Takes square root of every item in list. Homemade (slower) alt. to vecexpr.
#
# Arguments:
#       inputList               {list}      a list of numbers
#
# Results:
#
#       Returns a list of square roots.

proc vecSqrt {inputList} {
    set outputList []
    foreach item $inputList {
        lappend outputList [expr {sqrt($item)}]
    }
    return $outputList
}


# vecAtan2
#       Takes arctan2 of every item in list. Homemade (slower) alt. to vecexpr.
#
# Arguments:
#       Xlist               {list}      a list of X coordinates
#       Ylist               {list}      a list of Y coordinates
#
# Results:
#
#       Returns a list of arctan2s.

proc vecAtan2 {Xlist Ylist} {
    set outputList []
    foreach x $Xlist y $Ylist {
        lappend outputList [expr {atan2($y,$x)}]
    }
    return $outputList
}

# assignBins (Previously: bin_assigner)--
#       makes a list of bins based on x, y values and the coordinate system
# Arguments:
#       xVals       {list}      list of all x values of lipids in a particular species
#       yVals       {list}      list of all y values of lipids in a particular species
#       binWidth1   {float}     Length of x dimension for a particular frame
#       binWidth2   {float}     Length of y dimension for a particular frame
#       thetaDeg    {float}     Degree system split to get appropriate Ntheta or number of bins
#       polar       {int}       1 or 0 denoting if system is analyzed in polar or cartesian coordinartes
#       frm         {int}       current frame of system
#
# Results:
#       returns two lists of bins in the x or y direction

proc assignBins {xVals yVals binWidth1 binWidth2 thetaDeg polar frm use_vecexpr} {
    
    if {$polar == 1} {
        ;# use polar (r,theta) bins

        ;#calculate r: distance from origin for all x,y pairs
        set r_vals2 [vecadd [vecmul $xVals $xVals] [vecmul $yVals $yVals]]

        if {$use_vecexpr == "yes"} {
            set r_vals [vecexpr $r_vals2 sqrt]

            ;#calculate theta: use atan2 to get values for al x,y pairs
            set theta_vals [vecexpr $yVals $xVals atan2 pi div 180 mult]
        } else {
            set r_vals [vecSqrt $r_vals2]

            ;#calculate theta: use atan2 to get values for al x,y pairs
            global M_PI
            set theta_vals [vecscale [vecAtan2 $xVals $yVals] [expr 180/$M_PI]]
        }


        ;#atan2 gives values from -180 to 180; shifting to 0 to 360
        for {set i 0} {$i<[llength $theta_vals]} {incr i} {
            if {[lindex $theta_vals $i] < 0} {
                set theta_vals [lreplace $theta_vals $i $i [expr [lindex $theta_vals $i]+360]]
            }
        }

        ;#turn into bin numbers (floats)
        set dim1_bins_float [vecscale $r_vals [expr 1.0/$binWidth1]]
        set dim2_bins_float [vecscale $theta_vals [expr 1.0/$thetaDeg]]

        # floor all the floats to actually get bin numbers
        if {$use_vecexpr == "yes"} {
            set dim1_bins [vecexpr $dim1_bins_float floor]        
            set dim2_bins [vecexpr $dim2_bins_float floor]
        } else {
            set dim1_bins [vecFloor $dim1_bins_float]        
            set dim2_bins [vecFloor $dim2_bins_float]
            set dim1_bins [vecscale $dim1_bins 1.0]
            set dim2_bins [vecscale $dim2_bins 1.0]
        }

        
    } elseif {$polar == 0} {
        ;# use cartesian (x,y) bins
        
        ;# shift all values so that values within unitcell are greater than 0 and less than unitcell len
        set xlen [molinfo top get a frame $frm]
        set ylen [molinfo top get b frame $frm]
        set xmin [expr -$xlen/2.0]
        set ymin [expr -$ylen/2.0]
        set xVals [vecaddScalar $xVals [expr -1.0*$xmin]]
        set yVals [vecaddScalar $yVals [expr -1.0*$ymin]]

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

        ;# turn into bin numbers (floats)
        set dim1_bins_float [vecscale $xVals [expr 1.0/$binWidth1]]
        set dim2_bins_float [vecscale $yVals [expr 1.0/$binWidth2]]

        ;# floor all floats to actually get bin numbers
        if {$use_vecexpr == "yes"} {
            set dim1_bins [vecexpr $dim1_bins_float floor]
            set dim2_bins [vecexpr $dim2_bins_float floor]
        } else {
            set dim1_bins [vecFloor $dim1_bins_float]
            set dim1_bins [vecscale $dim1_bins 1.0]
            set dim2_bins [vecFloor $dim2_bins_float]
            set dim2_bins [vecscale $dim2_bins 1.0]
        }
    }

    return [list $dim1_bins $dim2_bins]
}

# createOutfiles (Previously: create_outfiles) --
#
# Arguments:
#       quantity        {str}       quanities being evaluated, either height_density or tilt_order
#       species         {list}      species of lipids in system
#       tailList        {list}      nested list of tail names organized by lipid species
#       folderName      {str}       name of folder
#
# Results:
#       returns a dict containing all the outfile names/addresses
#       seperates density by species and tilt and order segregate 
#       by species and tail number

proc createOutfiles {quantity species tailList folderName} {
    file mkdir "${folderName}/tcl_output"
    
    if {$quantity eq "height"} {
        set height "${folderName}/tcl_output/height"
        file mkdir $height
        dict set outfiles z1z2 heights_up fname [open "${height}/zone.dat" w]
        dict set outfiles z1z2 heights_down fname [open "${height}/ztwo.dat" w]
        dict set outfiles z0 heights_zzero fname [open "${height}/zzero.dat" w]
    } elseif {$quantity eq "density"} {
        set density "${folderName}/tcl_output/density"
        file mkdir $density
        file mkdir "${density}/combined"
        dict set outfiles z1z2 counts_up fname [open "${density}/combined/zone.dat" w]
        dict set outfiles z1z2 counts_down fname [open "${density}/combined/ztwo.dat" w]
        dict set outfiles z0 counts_zzero fname [open "${density}/combined/zzero.dat" w]
        foreach lipidtype $species {
            file mkdir "${density}/${lipidtype}"
            dict set outfiles z1z2 density_up_${lipidtype} fname [open "${density}/${lipidtype}/zone.dat" w]
            dict set outfiles z1z2 density_down_${lipidtype} fname [open "${density}/${lipidtype}/ztwo.dat" w]
            dict set outfiles z0 density_zzero_${lipidtype} fname [open "${density}/${lipidtype}/zzero.dat" w]
        }
    } elseif {$quantity eq "tilt_order"} {
        set order "${folderName}/tcl_output/order"
        set tilt "${folderName}/tcl_output/tilt"
        file mkdir $order
        file mkdir $tilt
        for {set i 0} {$i < [llength $tailList]} {incr i} {
            set lipidtype [lindex $species $i]
            for {set j 0} {$j < [llength [lindex $tailList $i]]} {incr j} {
                set tailnum "tail$j"
                set order_folder "${order}/${lipidtype}/${tailnum}"
                set tilt_folder "${tilt}/${lipidtype}/${tailnum}"
                file mkdir $order_folder
                file mkdir $tilt_folder
                set taillength [llength [lindex [lindex $tailList $i] $j]]
                dict set outfiles $taillength tilts_up_${lipidtype}_${tailnum} fname [open "${tilt_folder}/zone.dat" w]
                dict set outfiles $taillength tilts_down_${lipidtype}_${tailnum} fname [open "${tilt_folder}/ztwo.dat" w]
                dict set outfiles $taillength order_up_${lipidtype}_${tailnum} fname [open "${order_folder}/zone.dat" w]
                dict set outfiles $taillength order_down_${lipidtype}_${tailnum} fname [open "${order_folder}/ztwo.dat" w]
            }
        }
    }
    
    return $outfiles
}

# measureBoxSize (Prevoiusly: measure_box_size)--
#
#       measures length of a single boxedge
#
# Arguments:
#       frame        {int}       current frame being analyzed 
#       coordinate   {str}       coordinate, either "x","y" or "z" 
# Results:
#   
#       returns the length of the box in a single direction 
#
# Necessary Revisions/Problems: 
#       [molinfo top get a] is unreliable for certain coarse-grain simulations
# 
proc measureBoxSize {frame coordinate} {
    if {($coordinate ne "x") && ($coordinate ne "y") && ($coordinate ne "z")} {
        puts "dim must be x, y, or z"
        return
    }
    set sel [atomselect top all frame $frame]
    set dimvals [$sel get $coordinate]
    set max [lindex [lsort -real $dimvals] end]
    set min [lindex [lsort -real $dimvals] 0]
    set len [expr $max-$min]

    return $len
}

# prepareBins (Previously: bin_prep)--
#
#       determines the number of bins and the step length in each dimension
#
# Arguments:
#       frameNumber     {int}       current frame being analyzed 
#       polar           {int}       1 or 0 denoting if system is analyzed in polar or cartesian coordinartes
#       min             {int}       controls how much of the region around the protein is analyzed (not fully implemented)
#       drN1            {int}       bin dimension in x direction/radial bin dimensions
#       N2              {int}       bin dimension in y direction
#
# Results:
#
#       returns a dictionary containing bin and step length information
#
# Necessary Revisions/Problems:
#       - Assumes a square system for cartesian plotting
#       - min is only partially implemented

proc prepareBins {frameNumber polar min drN1 N2} {
    
    if {$polar == 1} {
    
        dict set bindims d1 $drN1

        #measure box size at final frame to get bin values
        set box_x [molinfo top get a frame [expr $frameNumber]]

        set box_r [expr int($box_x)/2]
        set rrange [expr $box_r-$min]
        
        #calculate number of dim1 bins from d1 and range1
        if {[expr $rrange%$drN1] == 0} { 
            dict set bindims N1 [expr [expr $rrange/$drN1]-1] 
        } else {
            dict set bindims N1 [expr $rrange/$drN1]
        }

        #calculate dtheta in degrees and radians
        dict set bindims dthetadeg [expr 360/[expr $N2*1.0]]
        global M_PI
        dict set bindims d2 [expr 2*$M_PI/$N2]
        dict set bindims N2 $N2

        polarAreaWarning [dict get $bindims d1] [dict get $bindims N1] $min [dict get $bindims d2]
    
    } elseif {$polar == 0} {

        dict set bindims N1 $drN1
        ;# fix this if you ever want to implement rectangular systems
        dict set bindims N2 $drN1
        
        set bindims [updateDimensions $bindims 0]

        dict set bindims dthetadeg "NULL" 
    
    }

    return $bindims
}

# polarAreaWarning --
#
#       Prints a warning if the bin area gets too small
#
# Arguments:
#       d1          {float}         length of r bin
#       N1          {int}           number of r bins
#       min         {float}         starting r value for bin 0
#       d2          {float}         length of theta bin
#
# Result:
#
#       Printed warning; no return.

proc polarAreaWarning {d1 N1 min d2} {
    set baseArea [expr $d1*$d2]
    set i 0
    set area 0
    while {($i < $N1) && ($area < 66.67)} {
        set distToCenter [expr [expr $min + $i * $d1] + [expr $d1 / 2.0]]
        set area [expr $baseArea * $distToCenter]
        incr i 
    }
    if {$i==0} {
        return
    } elseif {$i == $N1} {
        puts "WARNING: All bins are smaller than .67 nm^2"
        puts "Consider resizing your bins, or take results with a grain of salt."
    } else {
        puts "WARNING: Bins closer than edge of radial ring [expr $i-2] are smaller than .67 nm^2"
        puts "Consider resizing your bins, or take results with a grain of salt."
    }
}

# calculateReferenceHeight (Previously: calc_ref_height)--
#
#       Calculates reference height of inclusion 
#
# Arguments:
#       configDict      {dict}      dictionary of values created from the config file 
#       frame           {int}       current frame
#
# Result:
#
#       returns the reference height of the inclusion for a specific frame
  
proc calculateReferenceHeight {configDict frm} {
    if {[dict get $configDict reference_point] ne "NULL"} {
        set ref_sel [atomselect top [dict get $configDict reference_point] frame $frm]
        set ref_height [lindex [measure center $ref_sel] 2]
        $ref_sel delete
    } else {
        set ref_height "NULL"
    }
    return $ref_height
}

# getSelInfo (Previously: grab_sel_info)--
#
#       gets various pieces of information of lipid atomselection
#
# Arguments:
#       sel             {atomsel}       atomselection for lipids 
#       refHeight       {float}         reference height to subtract from lipids
#
# Results:
#
#       returns dictionary containing various pieces of information 
#       for lipids in system

proc getSelInfo {sel refHeight} {
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
    if {$refHeight ne "NULL"} {
        dict set sel_info zvals_list [vecaddScalar [$sel get z] [expr -1.0*$refHeight]]
    } else {
        dict set sel_info zvals_list [$sel get z]
    }   

    ;# user contains a 1 or 2 for outer or inner leaflet, respectively
    dict set sel_info leaflet_list [$sel get user]

    ;# user3 contains an int that describes which tail in the lipid this is
    ;# E.G. POPC will have 0 or 1 (it has two tails)
    ;# E.G. OANT will have 0, 1, 2, 3, 4, or 5 (it has 6 tails)
    set tail_list [$sel get user3]
    dict set sel_info tail_list [vecaddScalar $tail_list -1]

    return $sel_info
}

# updateDimensions (Previously: update_dims)--
#
#       Updates the bin dimensions based on frame
#
# Arguments: 
#       bindims         {dict}  previous bin dimensions
#       frame           {int}   current frame
#
# Results:
#   
#       returns updated bin dimensions for the current frame

proc updateDimensions {bindims frame} {
    set x [molinfo top get a frame $frame]
    set y [molinfo top get b frame $frame]

    set d1 [expr $x/[expr [dict get $bindims N1]*1.0]]
    set d2 [expr $y/[expr [dict get $bindims N2]*1.0]]

    dict set bindims d1 $d1 
    dict set bindims d2 $d2 

    if {[expr $d1*$d2] < 6.7} {
        puts "WARNING: bin size is less than .67 nm^2"
        puts "consider resizing your bins to be bigger."
    }

    return $bindims
}

# concatenateNames (Previously: concat_names) --
#
#       concatenate all beadnames together for file naming purposes
#
# Arguments: 
#       headNames       {str}       names of beads that define neutral surface
#
# Results:
#   
#       returns string of bead names

proc concatenateNames { headNames } {
    if {[llength $headNames] > 1} {
        set condensed_name [lindex $headNames 0]
        for {set i 1} {$i < [llength $headNames]} {incr i} {
            set addname [lindex $headNames $i]
            set condensed_name "$condensed_name.$addname"
        }
    } elseif {[llength $headNames] == 1} {
        set condensed_name [lindex $headNames 0]
    } else {
        puts "headnames must contain a bead name"
        break
    }
    return $condensed_name
}

# createResidueDictionaries (Previously: create_res_dict)--
#
#       creates dictionary containing various values needed downstream
#
# Arguments:
#       species         {list}      species of lipids in system
#       headNames       {str}       names of beads that define neutral surface
#       lipidList       {list}      all lipid resnames in system
#       nameList        {list}      all lipid bead names in system
#       dimOneBinList   {list}      list of bins in the x direction/radial bin's
#       dimTwoBinList   {list}      list ob bins in y direction
#       leafletList     {list}      list of user values for lipids denoting upper, lower, 
#                                   or middle leaflet
#       selex           {atomsel}   from the dictionary created by the createAtomSelections proc
#
# Results:
#  
#       returns dictionary with various values corresponding to lipids in 
#       a particular bin

proc createResidueDictionaries { species headNames lipidList nameList dimOneBinList dimTwoBinList leafletList selex} {
    ;# initialize a nested dict with a dummy key and value
    dict set res_dict dummy "dummy"
    
    if {$selex ne "z0"} {
        for {set i 0} {$i < [llength $lipidList]} {incr i} {
            if {([lindex $leafletList $i] == 3) || ([lindex $leafletList $i] == 4)} {
                continue
            } elseif {([lsearch $species [lindex $lipidList $i]] != -1) && ([lsearch $headNames [lindex $nameList $i]] != -1)} {
                set bin "[lindex $dimOneBinList $i],[lindex $dimTwoBinList $i]"
                set bin_leaf "$bin,[expr int([lindex $leafletList $i])]"
                if {[dict exists $res_dict $bin_leaf]} {
                    dict append res_dict $bin_leaf " $i"
                } else {
                    dict set res_dict $bin_leaf $i                
                }
            }
        }
    } else {
        ;# zzero needs to be handled separately; keys are set to "bin1#,bin2#,3"
        for {set i 0} {$i < [llength $lipidList]} {incr i} {
            set bin "[lindex $dimOneBinList $i],[lindex $dimTwoBinList $i]"
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


proc calc_bin_info {start end step N1 N2 coordSystem d1 d2} {
    set arealist []
    set d1list []
    set d2list []
    for {set frm $start} {$frm <= $end} {set frm [expr $frm+$step]} {
        set L1 [molinfo top get a frame $frm]
        set L2 [molinfo top get b frame $frm]
        lappend arealist [expr $L1*$L2]
        if {$coordSystem == "cart"} {
            lappend d1list [expr $L1/$N1*1.0]
            lappend d2list [expr $L2/$N2*1.0]
        }
    }
    set avgarea [vecmean $arealist]
    if {$coordSystem == "cart"} {
        set avgd1 [vecmean $d1list]
        set avgd2 [vecmean $d2list]
    } elseif {$coordSystem == "polar"} {
        set avgd1 $d1 
        set avgd2 $d2
    } else {
        puts "Something went wrong with calc_bin_info"
    }

    return [list $avgarea $avgd1 $avgd2]
}


proc outputNougatLog {start end step species tailList system headNames coordSystem folderName N1 N2 d1 d2} {
    set logFile [open "${folderName}/tcl_output/nougat.log" w]

    ;# calculate average area, d1, and d2
    set binInfo [calc_bin_info $start $end $step $N1 $N2 $coordSystem $d1 $d2]
    set avgArea [lindex $binInfo 0]
    set avgd1 [lindex $binInfo 1]
    set avgd2 [lindex $binInfo 2]

    ;# output system info
    puts $logFile "#NAME AND COORDINATE SYSTEM USED"
    puts $logFile "$system"
    puts $logFile "$coordSystem"
    puts $logFile ""

    ;# output species names and bead names
    puts $logFile "#SYSTEM CONTENTS"
    puts $logFile "$species"
    puts $logFile ""

    ;# output density normalization info
    puts $logFile "#NUMBER OF TAILS"
    for {set i 0} {$i < [llength $tailList]} {incr i} {
        set lipidtype [lindex $species $i]
        set ntails [llength [lindex $tailList $i]]
        puts $logFile "${lipidtype}:${ntails}"
    }
    puts $logFile ""

    ;# output density normalization info
    puts $logFile "#DENSITY NORMALIZATION"
    set density_norm_factor [outputDensityNormInfo $species $avgArea]
    foreach spec_norm_pair $density_norm_factor {
        puts $logFile "$spec_norm_pair"
    }
    puts $logFile ""

    ;# output bin number
    puts $logFile "#BIN INFO"
    puts $logFile "$N1 $N2"
    puts $logFile "$avgd1 $avgd2"
    puts $logFile ""

    close $logFile
}


# outputDensityNormInfo (Previously: output_density_norm_info)
#
#       Calculates the normalization factor for 
#       density enrichment calculations
#
# Arguments:
#       species         {list}      species of lipids in system
#       avgArea         {flt}       average area of box across portion of traj under analysis
#
# Results:
#       
#       calculates the normalization factor for density enrichment calculations

proc outputDensityNormInfo {species avgArea} {
    set normlist []
    foreach spec $species {
        set sel [atomselect top "resname $spec"]
        set beads_per_lipid [llength [lsort -unique [$sel get name]]]
        set number_of_lipids [llength [lsort -unique [$sel get resid]]]
        set lipids_in_one_leaflet [expr {int([expr $number_of_lipids / 2.0])}]
        $sel delete
        set number_of_beads [expr $number_of_lipids*$beads_per_lipid]
        set number_of_beads_per_leaflet [expr $number_of_beads / 2.0]
        set bulk_density [expr $number_of_beads_per_leaflet / $avgArea]
        set normfactor [expr 1 / $bulk_density]
        lappend normlist "${spec}:${normfactor}"
    }
    return $normlist
}


# sortTailLength (Previously: tail_length_sorter) --
#
#       Determines tails with unique tail lengths and 
#       makes an atomselection for them
#
# Arguments:
#       species         {list}      List of lipid types in system
#       acylNames       {list}      List of list of all bead names of all tails in the system 
#
# Results:
#   
#       Returns a list with an atomselection and the number of unique 
#       tails (tails with different lengths)

proc sortTailLength {species acylNames} {
    set lenlist []
    for {set i 0} {$i < [llength $acylNames]} {incr i} {
        foreach tail [lindex $acylNames $i] {
            lappend lenlist [llength $tail]
        }
    }
    set lengthlist [lsort -unique $lenlist]
    set sellist []
    foreach length $lengthlist {
        set resnamelist []
        set namelist []
        for {set i 0} {$i < [llength $acylNames]} {incr i} {
            foreach tail [lindex $acylNames $i] {
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
#       start       {vector}    Position of the start of the vector
#       end         {vector}    Position of the end of the vector
#
# Results:
#       returns the cosine theta of the vector, 
#       this is the dot product of n_{1,2} and the vector [0 0 1]
#       this corresponds to the 3rd value of n_{1,2} 
proc getCosineTheta {start end} {
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
#       length      {int}       length of lipids being evaluated
#       xValues     {list}      x coordinate of tail beads of all lipids being evaluated 
#       yValues     {list}      y coordinate of tail beads of all lipids being evaluated
#       zValues     {list}      z coordinate of tail beads of all lipids being evaluated
#
# Results:
#
#       Returns a nested list of the average order parameter of a lipid tail
#       ex.
#       {{0.4 0.4 0.4 0.4} {0.9 0.9 0.9 0.9} ...}  

proc calculateOrderParameters {length xValues yValues zValues} {
    set order_list []
    set temp_list []
    for {set i 1} {$i <= [llength $xValues]} {incr i} {
        if {[expr $i%$length] == 0} {
            ;# when this is TRUE, you've gotten cos2theta for each of the bonds in your tail
            ;# already and now you need to average them
            set avg [vecmean $temp_list]
            set order [expr {$avg * 1.5 - 0.5}]
            lappend order_list [lrepeat $length $order]
            set temp_list []
        } else {
            ;# calculate cos2theta for each bond in the tail and append to a temp list
            set start [list [lindex $xValues $i] [lindex $yValues $i] [lindex $zValues $i]]
            set end [list [lindex $xValues [expr $i-1]] [lindex $yValues [expr $i-1]] [lindex $zValues [expr $i-1]]]
            set costheta [getCosineTheta $start $end]
            lappend temp_list [expr {$costheta * $costheta}]
        }
    }

    ;# this is a list of lists, but we want just a list
    set final_order_list [concatenateList $order_list "NULL"]
    
    return [list $final_order_list]
}

# averageTiltAndOrderParameter (Previously: tilt_order_averaging)
#
#       Uses residueDictionary entries to compute bin averages, then assigns them to the correct outfile
#
# Arguments:
#       residueDictionary       {dict}      contains various attributes of the system (check:createResidueDictionaries)
#       outfiles                {dict}      contains all output files
#       lipidList               {list}      all lipid resnames in system 
#       tilts                   {list}      tilts of lipids 
#       orders                  {list}      order parameter of lipid tails
#       tailList                {list}      nested list of tails organized by lipid type (see analyzeTails)
#       selex                   {atomsel}   from the dictionary created by the createAtomSelections proc 
#
# Results:
#
#       Returns dictionary with average tilt and order parameter 
#       of lipids for a specific bin 

proc averageTiltAndOrderParameter {residueDictionary outfiles lipidList tilts orders tailList selex} {
    dict set counts placeholder "dummy"
    dict for {bin indices} $residueDictionary {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr {[string length $bin] - 3}]]
        foreach indx $indices {
            set tailnum [expr int([lindex $tailList $indx])]
            set species [lindex $lipidList $indx]
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
                set newtiltsum [vecadd [vecscale $oldtilt $oldcount] [lindex $tilts $indx]]
                set newtilt [vecscale $newtiltsum [expr 1.0/$newcount]]
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


# vecaddScalar
#
#       Adds a scalar to every item in list. Homemade (slower) alternative to vecexpr.
#
# Arguments:
#       inputList               {list}      a list of numbers
#       scalar                  {float}     the scalar you wish to add to each item of inputList
#
# Results:
#
#       Returns list that has been transformed by the scalar value provided.

proc vecaddScalar {inputList scalar} {
    set outputList []
    foreach item $inputList {
        lappend outputList [expr $item + $scalar]
    }
    return $outputList
}


# averageHeight (Previously: height_density_averaging)
#
#       uses residueDictionary entries to compute bin averages, then assigns them to the correct outfile
#
# Arguments:
#       residueDictionary       {dict}      contains various attributes of the system (check:createResidueDictionaries) 
#       outfiles                {dict}      contains all output files 
#       lipidList               {list}      all lipid resnames in system 
#       zValsList               {list}      height of lipids from inclusion if present
#
# Results:
#
#       Returns dictionary with average height, density and counts 
#       of lipids for a specific bin 

proc averageHeight {residueDictonary outfiles lipidList zValsList} {
    dict for {bin indices} $residueDictonary {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr {[string length $bin] - 3}]]
        foreach indx $indices {
            set species [lindex $lipidList $indx]
            if {$leaf == 1} {
                set field_key "z1z2"
                set height_key "heights_up"
                set counts_key "counts_up"
            } elseif {$leaf == 2} {
                set field_key "z1z2"
                set height_key "heights_down"
                set counts_key "counts_down"
            } elseif {$leaf == 3} {
                set field_key "z0"
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
                set newsum [expr {$oldavg * $oldcount + [lindex $zValsList $indx]}]
                set newavg [expr $newsum/$newcount]
                dict set outfiles $field_key $height_key bin $correct_bin $newavg
                dict set outfiles $field_key $counts_key bin $correct_bin $newcount   
            } else {
                dict set outfiles $field_key $height_key bin $correct_bin [lindex $zValsList $indx]
                dict set outfiles $field_key $counts_key bin $correct_bin 1.0
            }
        }
    }
    return $outfiles
}

# setBetaValues (Previously: set_beta_vals)--
#
#       assigns beta values to various components of system
#       
# Arguments:
#       selectInclusion     {str}       selection for the inclusion
#
# Results:
#
#       sets all beta values in vmd based on group

proc setBetaValues {selectInclusion} {
    puts "starting to fill beta values"
    
    if {$selectInclusion ne "NULL"} {
        set inclusion [atomselect top $selectInclusion]
        $inclusion set beta 0
        $inclusion delete 

        set excl_sel [atomselect top "not $selectInclusion and not resname W"]
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
#       Converts bianary value (1 or 0) to "polar" or "cart"
#
# Arguments:
#       polar    {int}      1 or 0 denoting if system is analyzed in polar or cartesian coordinartes
#
# Results:
#       
#       returns string denoting if system is analyzed in polar or cartesian    

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
#       creates dictionaries of atomselections  
#       
# Arguments:
#       quantity            {str}       quanities being evaluated, either height or tilt_order       
#       configDictionary    {dict}      dictionary of values created from the config file 
#
# Results:
#       
#       returns a dictionary with updated atomselections for height&density 
#       or tilt&order processing 
#       
# Necessary Revisions/Problems:
#       Min is not fully implemented

proc createAtomSelections {quantity configDictionary} {
    ;#atomselections setup as dict
    if {$quantity eq "height"} {
        dict set selections z1z2 [atomselect top "resname [dict get $configDictionary species] and name [dict get $configDictionary full_tails]"]
        dict set selections z0 [atomselect top "resname [dict get $configDictionary species] and ((user 1 and within 6 of user 2) or (user 2 and within 6 of user 1))"]
    } elseif {$quantity eq "tilt_order"} {
        set lists [sortTailLength [dict get $configDictionary species] [dict get $configDictionary acyl_names]]
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



# Align --
#
#       Alignment based off vmd alignment 
#       
# Arguments:
#       align_seltext       {str}       Selection text for aligning    
#	tilt_flag	    {bool}	Rotational fit around z axis only
# 					used to prevent spin/swivel motion of a membrane protein,
# 					without changing membrane orientation (no tilt)
#					Values: 1 (align all axes), 0 (don't align tilt)
#					Default: 1 (align all)
#	M_PI		    {float}	The value of pi.
#	molid		    {string} 	The molid to align. Defuault: top
# Results:
#       
#       The original trajector
#       
# Necessary Revisions/Problems:
#	Why doesn't it know what M_PI is unless we pass it as an argument?
#      
# Based on a script by Jérôme Hénin <jerome.henin@cnrs.fr>

proc Align { align_seltext {tilt_flag 1} {molid top}} {
    puts "Align start"
    global M_PI
    set nframes [molinfo $molid get numframes]
    set system [atomselect $molid "all"]
    set ref [atomselect $molid $align_seltext frame 0]
    set align_by [atomselect $molid "index [$ref list]"]
    for {set the_frame 1} {$the_frame < $nframes} {incr the_frame} {
        $align_by frame $the_frame
        $system frame $the_frame
        set TM [measure fit $align_by $ref]
        if {$tilt_flag==0} {
		    # Sine and cosine of z rotation are in first column of 4x4 matrix
		    set m00 [lindex $TM 0 0]
		    set m01 [lindex $TM 0 1]
		    # Negative sign to reverse rotation
		    # Use atan2 for safety
		    set alpha [expr {-180.0 / $M_PI * atan2($m01, $m00)}]
		    set TM [transaxis z $alpha deg]
        }
        $system move $TM
	    
    }
    $ref delete
    $system delete
    $align_by delete
    puts "Align end"
}


;#********************************;#
;# Liam scripts or custom scripts ;#
;#********************************;#




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
    return [convertRadianToDegree $theta]
}

;# Ouputs position of the centered protein in a membrane
;# accross both leaflets
;# only useful for analysis later - doesn't impact your polar density script
proc Protein_Position {name hnames chainNames folderName} {
    ;# in order to use this, must have your TMD chains separated and saved as occupancy 1

    set lastframe [expr [molinfo top get numframes]-1]

    set zone_sel [atomselect top "(name $hnames and user 1) and within 6 of name BB"]
    set zone_zvals [$zone_sel get z]
    set zone_Ht [vecmean $zone_zvals]
    $zone_sel delete

    set ztwo_sel [atomselect top "(name $hnames and user 2) and within 6 of name BB"]
    set ztwo_zvals [$ztwo_sel get z]
    set ztwo_Ht [vecmean $ztwo_zvals]
    $ztwo_sel delete

    set zmid_sel [atomselect top "((user 1 and within 6 of user 2) or (user 2 and within 6 of user 1)) and within 6 of name BB"]
    set zmid_zvals [$zmid_sel get z]
    set zmid_Ht [vecmean $zmid_zvals]
    $zmid_sel delete

    foreach ht [list $zone_Ht $ztwo_Ht $zmid_Ht $zmid_Ht] eqtxt [list "zone" "ztwo" "zzero" "zplus"] {
        puts "$eqtxt"
        set fout [open "${folderName}/tcl_output/${name}_helcoords_${eqtxt}.dat" w]
        puts $fout  "#These are the positions of your TMD helices in polar coords"
        foreach chnm $chainNames {
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

proc run_qwrap {wrap_sel inclusion_sel} {
    puts "${wrap_sel}"
    puts "Qwrap now running"

    setBetaValues $inclusion_sel
    qunwrap compound beta
    if {$inclusion_sel ne "NULL"} {
        qwrap compound beta center $inclusion_sel 
    }
    if {$inclusion_sel ne $wrap_sel} {
        qwrap compound beta center $wrap_sel
    }

    puts "run_qwrap finished!"
}

;# removes lipids from analysis (by setting user to 4)
;# current removal criteria:
;# within 5 angstroms of the pore center
;# within 30 angstroms of protein BB resid 30 (upper pore-lining residue)
proc pore_sorter_custom {frm species inc} {
    ;# "same resid as" doesn't work and "same residue as" has no meaning in CG trajs,
    ;# so have to do this in two steps
    if {$inc eq "5x29"} {
        set porebeads [atomselect top "(resname $species and within 12 of (name BB and resid 30)) or (resname $species and ((x*x+y*y) <= 25))" frame $frm]
    } elseif {$inc eq "7k3g"} {
        set porebeads [atomselect top "resname $species and (x*x+y*y <= 25)" frame $frm]
    }
    set badresids [lsort -unique [$porebeads get resid]]
    $porebeads delete
    if {[llength $badresids] != 0} {
        set porelipids [atomselect top "resname $species and resid $badresids" frame $frm]
        $porelipids set user 4.0
        $porelipids delete
    } 
}
