
# Converts radians to degrees
proc RtoD {r} {
    global M_PI
    return [expr $r*180.0/$M_PI]
}

# Returns a least squares fit for each lipid tail with N_beads = $length
# in the system. This allows you to use 1 instance of $lsqnormfactor only,
# rather than calculate it on the fly each time. 
# $lsqnormfactor is a list of values (i-(N-1)/2)
proc tilt_angles {length xvals yvals zvals} {
    set tilt_list []
    set lsqnormfactor [calc_lsq_normfactor $length]
    set xvec [fit_all_tails $length $xvals $lsqnormfactor]
    set yvec [fit_all_tails $length $yvals $lsqnormfactor]
    set zvec [fit_all_tails $length $zvals $lsqnormfactor]
    for {set i 0} {$i < [llength $xvec]} {incr i} {
        set vector "[lindex $xvec $i] [lindex $yvec $i] [lindex $zvec $i]"
        set norm [vecnorm $vector]

        # it is desirable to have a list that is the same length as $xvals
        # for ease of indexing in the binning process; hence lrepeat
        lappend tilt_list [lrepeat $length $norm]
    }

    # this is now a list of lists, but we want it all in one level
    set final_tilt_list [cat_list $tilt_list "NULL"]
    
    return $final_tilt_list
}

# returns a list of (i-(N-1)/2)
proc calc_lsq_normfactor { length } {
    set diff [expr $length-1]
    set d [expr 0.5*$diff] ;#normalization factor
    set I []

    #Make a list from 0 to length
    for {set k 0} {$k < $length} {incr k} {
        lappend I $k
    }

    set lsqnormfactor [vecexpr $I $d sub]

    return $lsqnormfactor
}

# Fit the points x to x = ai + b, i=0...N-1, and return the value 
# a = sum[ (i-(N-1)/2) * x_i] ; reference: Bevington
proc fit_all_tails {tail_length list_of_tail_coords lsqnormfactor} {
    set fit_values []
    set N_lipids [expr [llength $list_of_tail_coords] / $tail_length]
    set start_idx 0

    # iterate through list of all tail coords, separating them into one 
    # lipid tail at a time with lrange
    for {set i 0} {$i < $N_lipids} {incr i} {
        set start_idx [expr $i * $tail_length]
        set end_idx [expr $start_idx + $tail_length - 1]
        set coords [lrange $list_of_tail_coords $start_idx $end_idx]

        # Perform least squares fitting:
        # Multiply x_i and the differences (i-(N-1)/2). stack=1vec
        # Sum over the vector. stack=1scalar
        # Append to $fit_values
        lappend fit_values [vecexpr [vecexpr $coords $lsqnormfactor mult] sum]
    }
    return $fit_values
}

# Concatenates list items into one long string, separated by 
# spaces on either side of the $delimiter.
# "NULL" will result in no delimiter with a single space separating elements.
# "or" will also enclose list items in parentheses for use as 
# atomselection text. 
proc cat_list {inputlist delimiter} {
    if {$delimiter eq "or"} {
        set output "([lindex $inputlist 0])"
    } else {
        set output [lindex $inputlist 0]
    }
    
    for {set i 1} {$i < [llength $inputlist]} {incr i} {
        if {$delimiter eq "or"} {
            set element "([lindex $inputlist $i])"
            set output "${output} or ${element}"
        } elseif {$delimiter eq "NULL"} {
            set output "${output} [lindex $inputlist $i]"
        } else {
            set output "${output} $delimiter [lindex $inputlist $i]"
        }
    }
    return $output
}

# Returns a list of lists containing the starting beads and ending beads for a 
# given lipid's acyl chains.
# E.G. POPC start beads would be "C1A C1B" and end beads would be "C4A C4B"
proc heads_and_tails {species taillist} {
    for {set i 0} {$i < [llength $taillist]} {incr i} {
        set startbead []
        set endbead []
        foreach tail [lindex $taillist $i] {
            lappend startbead [lindex $tail 0]
            lappend endbead [lindex $tail end]
        }
        lappend startsellist [cat_list $startbead "NULL"]
        lappend endsellist [cat_list $endbead "NULL"]
    }

    return [list $startsellist $endsellist]
}

proc new_leaflet_check {frm species taillist window} {
    set heads_and_tails [heads_and_tails $species $taillist]
    set starts [lindex $heads_and_tails 0]
    set ends [lindex $heads_and_tails 1]

    for {set i 0} {$i < [llength $species]} {incr i} {
        set lipidtype [lindex $species $i]
        set total_sel [atomselect top "resname $lipidtype" frame $frm]
        set species_bead_num [llength [lsort -unique [$total_sel get name]]]
        set startnames "PO4"
        set endnames [lindex $ends $i]
        set numbeads [llength $endnames]
        set start_sel [atomselect top "resname $lipidtype and name $startnames" frame $frm]
        set end_sel [atomselect top "resname $lipidtype and name $endnames" frame $frm]
        set start_z [$start_sel get z]
        set end_z [$end_sel get z]
        $start_sel delete
        $end_sel delete
        set userlist []
        set counter 0
        for {set j 0} {$j < [llength $end_z]} {set j [expr $j+$numbeads]} {
            set avgendheight [vecexpr [lrange $end_z $j [expr $j+$numbeads-1]] mean]
            set avgheight [expr [lindex $start_z $counter] - $avgendheight]
            if {$avgheight > $window} {
                lappend userlist [lrepeat $species_bead_num 1.0]
            } elseif {$avgheight < -$window} {
                lappend userlist [lrepeat $species_bead_num 2.0]
            } else {
                lappend userlist [lrepeat $species_bead_num 3.0]
            }
            incr counter
        }
        set user_vals [cat_list $userlist "NULL"]
        $total_sel set user $user_vals
        $total_sel delete
    }
    set sel [atomselect top "name BB and resid 30" frame $frm]
    set com [measure center $sel]
    set x [lindex $com 0]
    set y [lindex $com 1]
    $sel delete
    set porelipids [atomselect top "(resname $species and same residue as within 12 of (name BB and resid 30)) or (resname $species and same residue as ((x-$x)*(x-$x)+(y-$y)*(y-$y) <= 25))" frame $frm]
    set badresids [lsort -unique [$porelipids get resid]]
    $porelipids delete
    if {[llength $badresids] != 0} {
        set porelipids [atomselect top "resname $species and resid $badresids" frame $frm]
        $porelipids set user 4.0
        $porelipids delete
    } 
}

# leaflet_check checks whether a lipid has flipped from one leaflet to another.
# It solely relies on whether the end-bead of each tail is above/below the first bead of each tail.
proc leaflet_check {frm species taillist} {

    set counter 0
    foreach lipidtype $species {
        set sel [atomselect top "resname $lipidtype and name [lindex [lindex [lindex $taillist $counter] 0] 0]"]
        set resids [$sel get resid]
        $sel delete
        set topsum [lrepeat [llength $resids] 0.0]
        set bottomsum [lrepeat [llength $resids] 0.0]
        foreach tail [lindex $taillist $counter] {
            set topsel [atomselect top "resname $lipidtype and name [lindex $tail 0]" frame 0]
            set bottomsel [atomselect top "resname $lipidtype and name [lindex $tail end]" frame 0]
            set topsum [vecexpr $topsum [$topsel get z] add]
            set bottomsum [vecexpr $bottomsum [$bottomsel get z] add] 
            set resids [$topsel get resid]
            set chains [$topsel get chain]
            $topsel delete
            $bottomsel delete
        } 
        set topavg [vecexpr $topsum [llength [lindex $taillist $counter]] div]
        set bottomavg [vecexpr $bottomsum [llength [lindex $taillist $counter]] div]
        set diff_top_to_bottom [vecexpr $topsum $bottomsum sub]
        for {set i 0} {$i < [llength $diff_top_to_bottom]} {incr i} {
            if {[expr abs([lindex $diff_top_to_bottom $i])] > 100} {
                set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
                $sel set chain "Z"
                $sel delete
            } else {
                if {[lindex $diff_top_to_bottom $i] > 0 && ([lindex $chains $i] == "L" || [lindex $chains $i] == "Z")} {
                    set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
                    $sel set chain "U"
                    $sel delete
                } elseif {[lindex $diff_top_to_bottom $i] < 0 && ([lindex $chains $i] == "U" || [lindex $chains $i] == "Z")} {
                    set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
                    $sel set chain "L"
                    $sel delete
                }
            }
        }
        incr counter
    }

    set sel [atomselect top "name BB and resid 30" frame $frm]
    set com [measure center $sel]
    set x [lindex $com 0]
    set y [lindex $com 1]
    $sel delete
    set porelipids [atomselect top "(resname $species and same residue as within 12 of (name BB and resid 30)) or (resname $species and same residue as ((x-$x)*(x-$x)+(y-$y)*(y-$y) <= 25))" frame $frm]
    set badresids [$porelipids get resid]
    $porelipids delete
    if {[llength $badresids] != 0} {
        set porelipids [atomselect top "resname $species and resid $badresids"]
        $porelipids set chain "Z"
        $porelipids delete
    }
    
    set upper [atomselect top "chain U" frame $frm]
    $upper set user 1
    $upper delete
    set lower [atomselect top "chain L" frame $frm]
    $lower set user 2
    $lower delete
    set bad_chains [atomselect top "chain Z" frame $frm]
    $bad_chains set user 3
    $bad_chains delete
}

proc get_theta {x y} {
    set tmp  [expr {atan2($y,$x)}]
    global M_PI
    if {$tmp < 0} {
        set theta [expr 2*$M_PI + $tmp]    
    } else {
        set theta $tmp
    }
    return [RtoD $theta]
}

;# Ouputs position of the centered protein in a membrane
;# accross both leaflets
;# only useful for analysis later - doesn't impact your polar density script
proc Protein_Position {name hnames tnames} {
    ;# in order to use this, must have your TMD chains separated and saved as occupancy 3
    set chain_names [list "A" "B" "C" "D" "E"]

    set lastframe [expr [molinfo top get numframes] -1]

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
                set sel [atomselect top "(chain ${chnm} and name BB and occupancy 1) and (z < [expr $ht + 5] and z > [expr $ht - 5])" frame $lastframe]
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

proc Center_System {inpt} {
    puts "${inpt}"
    puts "Center_System now running"
    ;# confirms your box is either square or paraelleogram-ish
    ;# will preform qwrap or pbc wrap depending

    set pbc_angles [molinfo top get {alpha beta gamma}]
    
    set sel [atomselect top "$inpt" frame 0]
    set com [measure center $sel weight mass]
    
    set counter_i 0
    ;# continues to try and recenter's box until ~ 0'ed out
    while {[expr abs([lindex $com 0])] > 1.0 &&  [expr abs([lindex $com 1])] > 1.0} {
        
        if {$counter_i > 5} {
            puts "Script was unable to converge system to (0,0,0)"
            puts "Please check your system visually, there may be"
            puts "unintended artifacts"
            $sel delete
            return
        }
        
        if {([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 1]!=90.0) && ([lindex $pbc_angles 2]!=90.0)} {
            puts "qwrap may not be optimal for your system...\n"
            puts "Running pbc wrap. To verify proper centering"
            puts "pbc wrap will be run multiple times" ; after 100
            foreach i {0 1 2 3} {
                pbc wrap -centersel "$inpt" -all
            }
        } else {
            #qunwrap compound none
            qwrap sel all center "$inpt" ;#center entire system at ~0,0,0
        }
        set com [measure center $sel weight mass]
        incr counter_i
    }
    $sel delete
    puts "Center_System finished!"
}

;# starts a new line in the print file that has the min/max r or x value for the bin, depending on if polar or cartesian
proc print_line_init {file number d1 min} {
    puts -nonewline $file "[format {%0.2f} [expr $number * $d1 + $min]]  [format {%0.2f} [expr ($number+1) * $d1 + $min]]  "
}

;# adds a value to the print file
proc print_value {file value end_line} {
    if {$end_line == 0} {
        puts -nonewline $file " $value" 
    } elseif {$end_line == 1} {
        puts $file " $value"
    } else {
        puts "Something went wrong - end_line should have value of 0 or 1 only"
        break
    }
}


;#print an entire array (usually 1 frame) to file
proc print_frame {N1 outfiles key d1 min N2 polar selex} {

    set file [dict get $outfiles $selex $key fname]

    if {$polar == 1} {
        set N2 [expr $N2 -1.0]
    }

    ;# starts new line in outfile with bin values
    for {set m 0.0} {$m <= $N1} {set m [expr $m + 1.0]} {
        print_line_init $file $m $d1 $min
        ;# prints bin values through ultimate value in one line
        for {set n 0.0} {$n <= $N2} {set n [expr $n + 1.0]} {
            if {[dict exists $outfiles $selex $key bin "$m,$n"]} {
                print_value $file [dict get $outfiles $selex $key bin "$m,$n"] 0
            } else {
                print_value $file "nan" 0
            }
        }
        ;# starts a new line
        print_value $file " " 1
    }
}

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

;# THIS IS FOR 5X29!
proc set_occupancy {molid} {

    set sel [atomselect $molid "name BB SC1 to SC4"]
    set residuelist [$sel get residue]
    set resmax [::tcl::mathfunc::max {*}$residuelist]
    $sel delete

    set list1 0

    for {set i 0} {$i < $resmax} {incr i} {
        set sel1 [atomselect $molid "residue $i"]
        set sel2 [atomselect $molid "residue [expr $i+1]"]

        set loc1 [lindex [$sel1 get {x y z}] 0]
        set loc2 [lindex [$sel2 get {x y z}] 0]

        set dist [vecdist $loc1 $loc2]

            if {$dist > 15} {
                lappend list1 $i
                lappend list1 [expr $i+1]
            }
        $sel1 delete
        $sel2 delete
    }

    set chars [list A B C D E]
    set k 0
    lappend list1 $resmax

    foreach {i j} $list1 {
        set sel [atomselect $molid "residue $i to $j"]
        $sel set chain [lindex $chars $k]
        $sel delete
        set k [expr $k+1]
    }

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


proc leaflet_sorter {species taillist analysis_start} {
    puts "Starting leaflet sorting"
    set counter 0
    foreach lipidtype $species {
        set sel [atomselect top "resname $lipidtype and name [lindex [lindex [lindex $taillist $counter] 0] 0]"]
        set resids [$sel get resid]
        $sel delete
        set topsum [lrepeat [llength $resids] 0.0]
        set bottomsum [lrepeat [llength $resids] 0.0]
        foreach tail [lindex $taillist $counter] {
            set topsel [atomselect top "resname $lipidtype and name [lindex $tail 0]" frame 0]
            set bottomsel [atomselect top "resname $lipidtype and name [lindex $tail end]" frame 0]
            set topsum [vecexpr $topsum [$topsel get z] add]
            set bottomsum [vecexpr $bottomsum [$bottomsel get z] add] 
            set resids [$topsel get resid]
            $topsel delete
            $bottomsel delete
        } 
        set topavg [vecexpr $topsum [llength [lindex $taillist $counter]] div]
        set bottomavg [vecexpr $bottomsum [llength [lindex $taillist $counter]] div]
        set diff_top_to_bottom [vecexpr $topsum $bottomsum sub]
        set uplist []
        set downlist []
        for {set i 0} {$i < [llength $diff_top_to_bottom]} {incr i} {
            if {[lindex $diff_top_to_bottom $i] > 0} {
                lappend uplist [lindex $resids $i]
            } elseif {[lindex $diff_top_to_bottom $i] < 0} {
                lappend downlist [lindex $resids $i]
            }
        }
        set upsel [atomselect top "resid $uplist"]
        set downsel [atomselect top "resid $downlist"]
        $upsel set user 1
        $upsel set chain "U"
        $downsel set user 2
        $downsel set chain "L"
        $upsel delete
        $downsel delete
        incr counter
    }

    for {set frm 0} {$frm <= $analysis_start} {incr frm} {
        leaflet_check $frm $species $taillist
    }
    
    puts "Leaflet sorting complete!"
}

# tail_analyzer returns $taillist, a nested list: 
# top level is by species
# mid level is by tail in species
# bottom level is by beads in tail
# e.g. a membrane with DO and DP lipids would be: 

# |-------------------------------------taillist------------------------------------|
#   |------------------DO-----------------| |------------------DP-----------------|
#     |-----tail0-----| |-----tail1-----|     |-----tail0-----| |-----tail1-----|
#
# { { {C1A C2A C3A C4A} {C1B C2B C3B C4B} } { {C1A D2A C3A C4A} {C1B D2B C3B C4B} } } 

proc tail_analyzer { species } {
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

    return $taillist
}

proc lipid_analyzer {species acyl_names} {
    set numlist []
    foreach lipid $species {
        set sel [atomselect top "resname $lipid"]
        set numbeads [llength [lsort -unique [$sel get name]]]
        $sel delete
        lappend numlist $numbeads
    }

    set starts_and_ends [heads_and_tails $species $acyl_names]
    set starts [lindex $starts_and_ends 0]
    set ends [lindex $starts_and_ends 1]
    
    ;# figure out how many tails are in each lipid
    set lenlist []
    for {set i 0} {$i < [llength $species]} {incr i} {
        lappend lenlist [llength [lindex $starts $i]]
    }
    set uniquelengths [lsort -unique $lenlist]

    foreach length $uniquelengths {
        set startsel []
        set endsel []
        set indices [lsearch -all $lenlist $length]
        foreach index $indices {
            set startseltext "resname [lindex $species $index] and name [lindex $starts $i]"
            set endseltext "resname [lindex $species $index] and name [lindex $ends $i]"
            lappend startsel $startseltext
            lappend endsel $endseltext
        }
        set uniquelengthstartseltext [cat_list $startsel "or"]
        set uniquelengthendseltext [cat_list $endsel "or"]
        dict set lipids $length start $uniquelengthstartseltext
        dict set lipids $length end $uniquelengthendseltext
    }

    return $lipids
}

;# tail_numberer changes user3 to hold a tail number.
;# This makes different tails easily separable for tilt/order analysis
proc tail_numberer { species taillist } {
    for {set lipidtype 0} {$lipidtype < [llength $species]} {incr lipidtype} {
        for {set tail 0} {$tail < [llength [lindex $taillist $lipidtype]]} {incr tail} {
            set sel [atomselect top "resname [lindex $species $lipidtype] and name [lindex [lindex $taillist $lipidtype] $tail]"]
            for {set frm 0} {$frm < [molinfo top get numframes]} {incr frm} {
                $sel frame $frm 
                $sel update
                $sel set user3 [expr $tail+1]
            }
            $sel delete
        }
    }
}

proc bin_assigner {x_vals y_vals d1 d2 dthetadeg polar} {
    
    if {$polar == 1} {
        ;# use polar (r,theta) bins

        ;#calculate r: distance from origin for all x,y pairs
        set r_vals [vecexpr [vecexpr [vecexpr $x_vals sq] [vecexpr $y_vals sq] add] sqrt]
        
        ;#turn into bin numbers rather than r values
        set dim1_bins [vecexpr [vecexpr $r_vals $d1 div] floor]
        
        ;#calculate theta: use atan2 to get values for al x,y pairs
        set theta_vals [vecexpr $y_vals $x_vals atan2 pi div 180 mult]

        ;#atan2 gives values from -180 to 180; shifting to 0 to 360
        for {set i 0} {$i<[llength $theta_vals]} {incr i} {
            if {[lindex $theta_vals $i] < 0} {
                set theta_vals [lreplace $theta_vals $i $i [expr [lindex $theta_vals $i]+360]]
            }
        }

        ;#turn into bin numbers rather than theta values
        set dim2_bins [vecexpr [vecexpr $theta_vals $dthetadeg div] floor]
        
    } elseif {$polar == 0} {
        ;# use cartesian (x,y) bins
        
        ;# shift all values so that they are temporarily positive
        ;# no negative bin numbers allowed!
        set xmin [vecexpr $x_vals min]
        set ymin [vecexpr $y_vals min]
        set x_vals [vecexpr $x_vals $xmin sub]
        set y_vals [vecexpr $y_vals $ymin sub]

        ;# turn into bin numbers rather than x,y values
        set dim1_bins [vecexpr [vecexpr $x_vals $d1 div] floor]
        set dim2_bins [vecexpr [vecexpr $y_vals $d2 div] floor]
    }

    return [list $dim1_bins $dim2_bins]
}

proc create_outfiles {system quantity_of_interest headnames species taillist coordsys} {

    if {$quantity_of_interest eq "height_density"} {
        dict set outfiles z1z2 heights_up fname [open "${system}.zone.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles z1z2 heights_down fname [open "${system}.ztwo.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles z0 heights_zzero fname [open "${system}.zzero.${headnames}.${coordsys}.height.dat" w]
        foreach lipidtype $species {
            dict set outfiles z1z2 density_up_${lipidtype} fname [open "${system}.${lipidtype}.zone.${coordsys}.density.dat" w]
            dict set outfiles z1z2 density_down_${lipidtype} fname [open "${system}.${lipidtype}.ztwo.${coordsys}.density.dat" w]
            dict set outfiles z0 density_zzero_${lipidtype} fname [open "${system}.${lipidtype}.zzero.${coordsys}.density.dat" w]
        }
    } elseif {$quantity_of_interest eq "tilt_order"} {
        for {set i 0} {$i < [llength $taillist]} {incr i} {
            set lipidtype [lindex $species $i]
            for {set j 0} {$j < [llength [lindex $taillist $i]]} {incr j} {
                set tailnum "tail$j"
                set taillength [llength [lindex [lindex $taillist $i] $j]]
                dict set outfiles $taillength tilts_up_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.zone.${coordsys}.tilt.dat" w]
                dict set outfiles $taillength tilts_down_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.ztwo.${coordsys}.tilt.dat" w]
                dict set outfiles $taillength order_up_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.zone.${coordsys}.order.dat" w]
                dict set outfiles $taillength order_down_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.ztwo.${coordsys}.order.dat" w]
            }
        }
    }
    
    return $outfiles
}

proc bin_prep {nframes polar min d1 N2} {
    #measure box size at final frame to get bin values
    set box_x [molinfo top get a frame [expr $nframes - 1]]

    if {$polar == 1} {
        set box_r [expr int($box_x) / 2]
        set range1 [expr $box_r - $min]
    } elseif {$polar == 0} {
        set range1 [expr int([vecexpr $box_x floor])]
    }
    
    #calculate number of dim1 bins from d1 and range1
    if {[expr $range1 % $d1] == 0} { 
        set N1 [expr [expr $range1 / $d1] - 1] 
    } else {
        set N1 [expr $range1 / $d1]
    }

    #calculate dim2 values, based on whether polar or cartesian
    if {$polar == 1} {
        set dthetadeg [expr 360/$N2]
        global M_PI
        set d2 [expr 2 * $M_PI / $N2]
    } elseif {$polar == 0} {
        set d2 $d1
        set box_y [molinfo top get b frame [expr $nframes - 1]]
        set range2 [expr int([vecexpr $box_y floor])]
        if {[expr $range2 % $d2] == 0} { 
            set N2 [expr [expr $range2 / $d2] - 1] 
        } else {
            set N2 [expr $range2 / $d2]
        }
        set dthetadeg "NA"
    }
    return [list $d1 $d2 $N1 $N2 $dthetadeg]
}

proc concat_names { headnames } {
    ;# concatenate all beadnames together for file naming purposes
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

proc create_res_dict { species headnames lipid_list name_list resid_list dim1_bins_list dim2_bins_list leaflet_list selex} {
    dict set res_dict dummy "dummy"
    if {$selex ne "z0"} {
        for {set i 0} {$i < [llength $lipid_list]} {incr i} {
            if {([lsearch $species [lindex $lipid_list $i]] != -1) && ([lsearch $headnames [lindex $name_list $i]] != -1)} {
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
    dict unset res_dict dummy
    return $res_dict
}

proc output_density_norm_info {start nframes step species system headnames coordsys} {
    set arealist []
    for {set frm $start} {$frm < $nframes} {set frm [expr $frm + $step]} {
        lappend arealist [expr [molinfo top get a frame $frm] * [molinfo top get b frame $frm]]
    }
    set avgarea [vecexpr $arealist mean]
    set normfactor_outfile [open "${system}.${coordsys}.density.normfactor.dat" w]
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
        set normfactor [expr $avgarea / [expr $Nb * $Sb / 2.0]]
        puts $normfactor_outfile "$spec $normfactor"
    }
    
    close $normfactor_outfile
}

proc tail_length_sorter {species acyl_names} {
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

proc get_costheta {i xvals yvals zvals} {
    set start [list [lindex $xvals $i] [lindex $yvals $i] [lindex $zvals $i]]
    set end [list [lindex $xvals [expr $i-1]] [lindex $yvals [expr $i-1]] [lindex $zvals [expr $i-1]]]
    set r12 [vecsub $start $end]
    set n12 [vecnorm $r12]
    return [lindex $n12 2]
}

proc order_params {length xvals yvals zvals} {
    set order_list []
    set temp_list []
    for {set i 1} {$i <= [llength $xvals]} {incr i} {
        if {[expr $i % $length] == 0} {
            set avg [vecexpr $temp_list mean]
            set order [expr $avg*1.5 - 0.5]
            lappend order_list $order
            set temp_list []
        } else {
            set costheta [get_costheta $i $xvals $yvals $zvals]
            lappend temp_list [expr $costheta * $costheta]
        }
    }
    set final_order_list []
    foreach item $order_list {
        set final_order_list [concat $final_order_list [lrepeat $length $item]]
    }
    return [list $final_order_list]
}

proc do_tilt_order_binning {res_dict outfiles leaflet_list lipid_list tilts orders tail_list selex} {
    dict set counts placeholder "dummy"
    dict for {bin indices} $res_dict {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr [string length $bin] - 3]]

        foreach indx $indices {
            set tailnum [expr int([lindex $tail_list $indx])]
            set species [lindex $lipid_list $indx]
            if {$leaf == 1} {
                set tilt_key "tilts_up_${species}_tail${tailnum}"
                set order_key "order_up_${species}_tail${tailnum}"
            } else {
                set tilt_key "tilts_down_${species}_tail${tailnum}"
                set order_key "order_down_${species}_tail${tailnum}"
            }
            if {[dict exists $counts $selex $order_key bin $correct_bin]} {
                set oldcount [dict get $counts $selex $order_key bin $correct_bin]
                set newcount [expr $oldcount + 1.0]
                set oldorder [dict get $outfiles $selex $order_key bin $correct_bin]
                set newordersum [expr $oldorder * $oldcount + [lindex [lindex $orders 0] $indx]]
                set neworder [expr $newordersum / $newcount]
                set oldtilt [dict get $outfiles $selex $tilt_key bin $correct_bin]
                set newtiltsum [vecexpr [vecexpr $oldtilt $oldcount mult] [lindex [lindex $tilts 0] $indx] add]
                set newtilt [vecexpr $newtiltsum $newcount div]
                dict set outfiles $selex $order_key bin $correct_bin $neworder
                dict set counts $selex $order_key bin $correct_bin $newcount
                dict set outfiles $selex $tilt_key bin $correct_bin $newtilt
            } else {
                dict set outfiles $selex $order_key bin $correct_bin [lindex [lindex $orders 0] $indx]
                dict set outfiles $selex $tilt_key bin $correct_bin [lindex [lindex $tilts 0] $indx]
                dict set counts $selex $order_key bin $correct_bin 1.0
            }
        }
    }
    dict unset counts $selex
    return $outfiles
}

proc do_height_density_binning {res_dict outfiles leaflet_list lipid_list zvals_list name_list} {
    dict for {bin indices} $res_dict {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr [string length $bin] - 3]]

        foreach indx $indices {
            set species [lindex $lipid_list $indx]
            if {$leaf == 1} {
                set field_key "z1z2"
                set dens_key "density_up_${species}"
                set height_key "heights_up"
            } elseif {$leaf == 2} {
                set field_key "z1z2"
                set dens_key "density_down_${species}"
                set height_key "heights_down"
            } elseif {$leaf == 3} {
                set field_key "z0"
                set dens_key "density_zzero_${species}"
                set height_key "heights_zzero"
            } else {
                puts "something has gone wrong with the binning"
                return
            }
            if {[dict exists $outfiles $field_key $height_key bin $correct_bin]} {
                set oldcount [dict get $outfiles $field_key $dens_key bin $correct_bin]
                set newcount [expr $oldcount + 1.0]
                set oldavg [dict get $outfiles $field_key $height_key bin $correct_bin]
                set newsum [expr $oldavg * $oldcount + [lindex $zvals_list $indx]]
                set newavg [expr $newsum / $newcount]
                dict set outfiles $field_key $height_key bin $correct_bin $newavg
                dict set outfiles $field_key $dens_key bin $correct_bin $newcount
            } else {
                dict set outfiles $field_key $height_key bin $correct_bin [lindex $zvals_list $indx]
                dict set outfiles $field_key $dens_key bin $correct_bin 1.0
            }
        }
    }
    return $outfiles
}