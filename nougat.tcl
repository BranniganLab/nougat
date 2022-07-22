
package require pbctools

# EDIT THE PATHS HERE
# TELL nougat WHERE TO FIND YOUR VERSIONS OF qwrap AND vecexpr
#set QWRAP "~/qwrap-master"
#set VEC "~/utilities/vecexpr"

set UTILS "~/PolarHeightBinning/utilities"

load ${UTILS}/qwrap.so
load ${UTILS}/vecexpr.so

proc cell_prep {system end} {

    ;# set lastframe based on $end input
    if {$end == -1} {
        set lastframe [expr [molinfo top get numframes] -1]
    } else {
        set lastframe $end
    }
    
    ;#**********************************************************
    ;#          MAKE EDITS BELOW BEFORE STARTING
    ;#********************************************************** 

    ;# provide atomselection-style text that defines what is in your inclusion 
    set inclusion_sel "name BB SC1 to SC4"

    ;# provide atomselection-style text that defines anything that isn't your inclusion_sel 
    ;# or membrane
    ;# E.G. solvent, ions, other molecules that aren't membrane lipids
    set excluded_sel "resname W ION"

    ;# figures out which lipids are in the system
    ;# no edits required
    set lipidsel [atomselect top "not $inclusion_sel and not $excluded_sel"]
    set species [lsort -unique [$lipidsel get resname]]
    $lipidsel delete

    ;# provide atomselection-style text that defines what bead(s) should be centered and wrapped around
    ;# usually, this would be name BB for proteins
    ;# for 5x29 we had absolute position restraints and a small box z dimension, so I'm using the membrane itself here
    set wrap_sel "resname $species"

    ;# provide atomselection-style text that defines what beads to align around if you want to prevent xy rotation from interfering with results
    ;# if your inclusion tumbles in the membrane (like a nanoparticle), comment out the align command below
    set align_sel "name BB"

    ;# provide atomselection-style text that defines the reference point that should correspond with height 0 in your plots
    ;# E.G. for 5x29 we decided resid 15 would be the 'zero-point' and all heights would be provided with reference to 
    ;# the average position of resid 15
    set reference_point "name BB and resid 15"

    ;# provide the beadnames that you consider to form the surface of your membrane
    ;# we chose the top tail beads because they are what form the 'hydrophobic surface'
    ;# in our opinion
    set headnames "C1A C1B"

    ;# center, wrap, and align the system
    ;# if your inclusion 'tumbles' in the membrane (like a nanoparticle) comment out Align!
    Center_System "$wrap_sel"
    Align "$align_sel"

    ;# custom proc to set my TMD helices to occupancy 1
    ;# this allows Protein_Position to work
    ;# comment this out or customize it for your inclusion
    set_occupancy top 

    ;# figures out which beads are in your lipids
    ;# then assigns lipids to user value 1 or 2, depending on leaflet
    ;# no edits required
    set acyl_names [tail_analyzer $species]
    set tail_ends [list_ends $acyl_names]
    leaflet_sorter $species $acyl_names $lastframe

    ;# this will only work if your TMD helices are set to occupancy 1
    ;# otherwise, comment it out
    ;# all it does is put a dot on the polar heatmap where a TMD helix should be, so not essential at all
    ;#Protein_Position $system $headnames $acyl_names ;# FIX ME

    ;#**********************************************************
    ;#          MAKE EDITS ABOVE BEFORE STARTING
    ;#********************************************************** 

    set return_list [] 
    lappend return_list $species 
    lappend return_list $headnames 
    lappend return_list $acyl_names
    lappend return_list $tail_ends
    lappend return_list $reference_point
    return $return_list  
}


proc list_ends {acyl_names} {
    set tail_list []
    foreach lipid $acyl_names {
        set lipidlist []
        foreach tail $lipid {
            lappend lipidlist [lindex $tail end]
        }
        lappend tail_list $lipidlist
    }
    return $tail_list
}

proc RtoD {r} {
    global M_PI
    return [expr $r*180.0/$M_PI]
}

#  returns a least squares fit for each lipid tail in the system at once
#  Fit the points x to x = ai + b, i=0...N-1, and return the value of a 
# a = 12/( (N(N^2 - 1)) ) sum[ (i-(N-1)/2) * xi]
# reference: Bevington
proc lsq_vecexpr { tail_length list_of_tail_coords } {
  set i_list []
  set counter 0
  for {set i 0} {$i < [llength $list_of_tail_coords]} {incr i} {
    lappend i_list $counter
    set counter [expr $counter + 1]
    if {[expr $counter%$tail_length]==0} {
      set counter 0
    }
  }
  set d [expr {0.5*($tail_length-1)}]
  vecexpr $i_list $d sub >multiplier
  set vector_components [vecexpr $list_of_tail_coords $multiplier mult]
  set vectors []
  for {set i 0} {$i < [expr [llength $list_of_tail_coords] / $tail_length]} {incr i} {
    set startnum [expr $i*$tail_length]
    set endnum [expr [expr [expr $i+1] * $tail_length] -1]
    set quantity_to_sum [lrange $vector_components $startnum $endnum]
    lappend vectors [vecexpr $quantity_to_sum sum]
  }

  return $vectors
}

#written in 2022 to be fast and dirty
# just relies on whether the bottom of each tail is above/below the top of each tail
# to sort leaflets
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
            if {[expr abs([lindex $diff_top_to_bottom $i])] > 40} {
                set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
                $sel set chain "U"
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

    #remove pore lipids
    #rewrite to make this useable for you

    #set sel [atomselect top "name BB and resid 30" frame $frm]
    #set com [measure center $sel]
    #set x [lindex $com 0]
    #set y [lindex $com 1]
    #$sel delete
    #set porelipids [atomselect top "(resname $species and same resid as within 9 of resid 30) or (resname $species and same resid as ((x-$x)*(x-$x)+(y-$y)*(y-$y) <= 16))" frame $frm]
    #$porelipids set chain "Z"
    #$porelipids delete

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
            qwrap sel "not $inpt" center "$inpt" ;#center entire system at ~0,0,0
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
proc print_frame {N1 file d1 min N2 data polar} {
    
    #needed in order to make the proc understand this is array data
    array set data_copy $data

    if {$polar == 1} {
        set N2 [expr $N2-1]
    }

    ;# starts new line in outfile with bin values
    for {set m 0} {$m <= $N1} {incr m} {
        print_line_init $file $m $d1 $min
        ;# adds bin values through penultimate value in one line
        for {set n 0} {$n < $N2} {incr n} {
            print_value $file $data_copy($m,$n) 0
        }
        ;# adds final value and starts new line in outfile
        print_value $file $data_copy($m,$N2) 1
    }
}


;# sets all values of array to zero
proc initialize_array {N1 N2 init_value} {

    for {set m 0} {$m <= $N1} {incr m} {
        for {set n 0} {$n <= $N2} {incr n} {
            set data($m,$n) $init_value
        }
    }

    return [array get data]
}


;# calculate bin density
proc calculate_density {data min d1 d2 N1 N2 totfrms polar} {
    
    array set data_copy $data

    if {$polar == 1} {
        for {set m 0} {$m <= $N1} {incr m} {
            set rval [expr $m * $d1 + $min + [expr 0.5 * $d1]]
            set area [expr $rval * $d1 * $d2]
            for {set n 0} {$n <= $N2} {incr n} {
                set data_copy($m,$n) [expr $data_copy($m,$n) / [expr $totfrms * $area]]
            }
        }
    } elseif {$polar == 0} {
        set area [expr $d1 * $d2]
        set divisor [expr $area * $totfrms]
        for {set m 0} {$m <= $N1} {incr m} {
            for {set n 0} {$n <= $N2} {incr n} {
                set data_copy($m,$n) [expr $data_copy($m,$n) / $divisor]
            }
        }
    }

    return [array get data_copy]
}


;# normalize bin density to get enrichment
proc normalize_density {data N1 N2 normfactor numbeads} {
    
    array set data_copy $data

    set normfactor [expr $normfactor * $numbeads]

    for {set m 0} {$m <= $N1} {incr m} {
        for {set n 0} {$n <= $N2} {incr n} {
            set data_copy($m,$n) [expr $data_copy($m,$n) / $normfactor]
        }
    }

    return [array get data_copy]
}

;# Alignment based off vmd alignment
proc Align { stuff } {
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
}

proc calc_vec_magn {xvals yvals zvals} {
    set x2 [vecexpr $xvals sq]
    set y2 [vecexpr $yvals sq]
    set z2 [vecexpr $zvals sq]
    return [vecexpr [vecexpr [vecexpr $x2 $y2 add] $z2 add] sqrt]
}

proc calc_vector_btw_beads {b1x b1y b1z b2x b2y b2z} {
    set xvals [vecexpr b2x b1x sub]
    set yvals [vecexpr b2y b1y sub]
    set zvals [vecexpr b2z b1z sub]
    return [list $xvals $yvals $zvals]
}

proc calc_avg_order {angle_list} {
    set threecos2minusone [vecexpr [vecexpr [vecexpr [vecexpr $angle_list cos] sq] 3 mult] 1 sub]
    return [expr [vecexpr $threecos2minusone mean] / 2.0]
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

proc tail_analyzer { species nframes } {
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
    
    ;# change user3 to match tail number
    ;# makes tails separable for tilt/order analysis
    for {set lipidtype 0} {$lipidtype < [llength $species]} {incr lipidtype} {
        for {set tail 0} {$tail < [llength [lindex $taillist $lipidtype]]} {incr tail} {
            set sel [atomselect top "resname [lindex $species $lipidtype] and name [lindex [lindex $taillist $lipidtype] $tail]"]
            for {set frm 0} {$frm < $nframes} {incr frm} {
                $sel frame $frm 
                $sel update
                $sel set user3 [expr $tail+1]
            }
            $sel delete
        }
    }

    return $taillist
}

proc tilt_angles {tail_length tail_one tail_two} {
    set vector_list []
    set t1xvals [lsq_vecexpr $tail_length [$tail_one get x]]
    set t1yvals [lsq_vecexpr $tail_length [$tail_one get y]]
    set t1zvals [lsq_vecexpr $tail_length [$tail_one get z]]
    set t2xvals [lsq_vecexpr $tail_length [$tail_two get x]]
    set t2yvals [lsq_vecexpr $tail_length [$tail_two get y]]
    set t2zvals [lsq_vecexpr $tail_length [$tail_two get z]]
    for {set i 0} {$i < [llength $t1xvals]} {incr i} {
        set vector1 "[lindex $t1xvals $i] [lindex $t1yvals $i] [lindex $t1zvals $i]"
        set norm1 [vecnorm $vector1]
        set vector2 "[lindex $t2xvals $i] [lindex $t2yvals $i] [lindex $t2zvals $i]"
        set norm2 [vecnorm $vector2]
        lappend vector_list $norm1
        lappend vector_list $norm2
    }
    return $vector_list
}

proc prep_thickness_lists {sellist ref_height} {
    set thickness_list []
    set counter 0
    foreach bead [lrange $sellist 4 end] {
        set thickness_list "$thickness_list [vecexpr [$bead get z] $ref_height sub]"
        incr counter
    }
    if {$counter > 1} {
        set mins_list [vecexpr $thickness_list $counter min_ew]
        ;# there is no max_ew in vecexpr, so multiply by -1, do min_ew, then multiply by -1 again
        set maxs_list [vecexpr [vecexpr [vecexpr $thickness_list -1.0 mult] $counter min_ew] -1.0 mult] 
    } elseif {$counter == 1} {
        set mins_list $thickness_list
        set maxs_list $thickness_list
    }
    return [list $mins_list $maxs_list]
}

proc bin_generator {x_vals y_vals d1 d2 dthetadeg polar} {
    if {$polar == 1} {
        ;#get theta values for all x,y pairs
        set theta_vals [vecexpr $y_vals $x_vals atan2 pi div 180 mult]

        ;#atan2 gives values from -180 to 180; shifting to 0 to 360
        for {set i 0} {$i<[llength $theta_vals]} {incr i} {
            if {[lindex $theta_vals $i] < 0} {
                set theta_vals [lreplace $theta_vals $i $i [expr [lindex $theta_vals $i]+360]]
            }
        }

        ;#turn into bin numbers rather than theta values
        vecexpr [vecexpr $theta_vals $dthetadeg div] floor &dim2_bins
        
        ;#calculate distance from origin for all x,y pairs
        set r_vals [vecexpr [vecexpr [vecexpr $x_vals sq] [vecexpr $y_vals sq] add] sqrt]
        
        ;#turn into bin numbers rather than r values
        vecexpr [vecexpr $r_vals $d1 div] floor &dim1_bins
    } elseif {$polar == 0} {
        set xmin [vecexpr $x_vals min]
        set ymin [vecexpr $y_vals min]
        set x_vals [vecexpr $x_vals $xmin sub]
        set y_vals [vecexpr $y_vals $ymin sub]
        vecexpr [vecexpr $x_vals $d1 div] floor &dim1_bins
        vecexpr [vecexpr $y_vals $d2 div] floor &dim2_bins
    }
    return [list $dim1_bins $dim2_bins]
}

proc create_outfiles {quantity_of_interest system headnames species taillist coordsys} {
    if {$quantity_of_interest eq "height"} {
        dict set outfiles heights_up [open "${system}.zone.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles heights_down [open "${system}.ztwo.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles heights_zplus [open "${system}.zplus.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles heights_zzero [open "${system}.zzero.${headnames}.${coordsys}.height.dat" w]
    } elseif {$quantity_of_interest eq "density"} {
        foreach lipidtype $species {
            dict set outfiles density_up_$lipidtype [open "${system}.${lipidtype}.zone.${headnames}.${coordsys}.density.dat" w]
            dict set outfiles density_down_$lipidtype [open "${system}.${lipidtype}.ztwo.${headnames}.${coordsys}.density.dat" w]
        }
    } elseif {$quantity_of_interest eq "tilt"} {
        for {set i 0} {$i < [llength $taillist]} {incr i} {
            set lipidtype [lindex $species $i]
            for {set j 0} {$j < [llength [lindex $taillist $i]]} {incr j} {
                set tailnum "tail$j"
                foreach bead [lindex [lindex $taillist $i] $j] {
                    dict set outfiles tilts_up_$lipidtype_$tailnum_$bead [open "${system}.${lipidtype}.${tailnum}.zone.${bead}.${coordsys}.tilt.dat" w]
                    dict set outfiles tilts_down_$lipidtype_$tailnum_$bead [open "${system}.${lipidtype}.${tailnum}.ztwo.${bead}.${coordsys}.tilt.dat" w]
                    dict set outfiles order_up_$lipidtype_$tailnum_$bead [open "${system}.${lipidtype}.${tailnum}.zone.${bead}.${coordsys}.order.dat" w]
                    dict set outfiles order_down_$lipidtype_$tailnum_$bead [open "${system}.${lipidtype}.${tailnum}.ztwo.${bead}.${coordsys}.order.dat" w]
                }
            }
        } 
    }
    return $outfiles
}

proc bin_prep {nframes polar min d1} {
    #measure box size to get bin values
    set box_x [molinfo top get a frame [expr $nframes - 1]]
    set min 0
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

proc create_atomselections {quantity_of_interest system beadname species acyl_names coordsys} {
    if {$quantity_of_interest == "height_density"} {
        dict set selections z1z2 [atomselect top "resname $species and name $beadname"]
        dict set selections z0 [atomselect top "(user 1 and within 6 of user 2) or (user 2 and within 6 of user 1)"]
    } elseif {$quantity_of_interest == "tilt_order_thickness"} {
        foreach lipidtype $species beadlist $acyl_names {
            foreach tail $beadlist {
                foreach bead $tail {
                    dict set selections $species.$bead [atomselect top "resname $species and name $bead"]
                }
            }
        }
    }
}

;########################################################################################
;# polarHeight Functions

proc start_nougat {system d1 N2 start end step polar} {

    set important_variables [cell_prep $system $start]
    set min 0 ;# change this value if you want to exclude an inner radius in polar coords
    set headnames [lindex $important_variables 1]
    set reference_point [lindex $important_variables 4]

    #need to calculate heights relative to some point on the inclusion:
    set ref_bead [atomselect top "$reference_point"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [vecexpr $ref_height mean]

    ;# set nframes based on $end input
    set maxframes [molinfo top get numframes]
    if {$end == -1} {
        set nframes $maxframes
    } elseif {($end <= $maxframes) || ($end > $start)} {
        set nframes $end
    } else {
        puts "you specified a frame number that is outside the allowable range"
        break
    }

    ;# determine number and size of bins
    set bindims [bin_prep $nframes $polar $min $d1]

    lappend important_variables $start 
    lappend important_variables $nframes 
    lappend important_variables $step 
    lappend important_variables $ref_height
    lappend important_variables $min

    run_nougat $system $headnames $important_variables $bindims $polar "height"
    run_nougat $system $headnames $important_variables $bindims $polar "density"
    run_nougat $system $headnames $important_variables $bindims $polar "tilt"
}

proc run_nougat {system beadname important_variables bindims polar quantity_of_interest} {  

    set boxarea []

    set species [lindex $important_variables 0]
    set headnames [lindex $important_variables 1]
    set acyl_names [lindex $important_variables 2]
    set tail_ends [lindex $important_variables 3]
    set start [lindex $important_variables 5]
    set nframes [lindex $important_variables 6]
    set step [lindex $important_variables 7]
    set ref_height [lindex $important_variables 8]
    set min [lindex $important_variables 9]
    
    set d1 [lindex $bindims 0]
    set d2 [lindex $bindims 1]
    set N1 [lindex $$bindims 2]
    set N2 [lindex $$bindims 3]
    set dthetadeg [lindex $bindims 4]

    if {$polar == 1} {
        set coordsys "polar"
    } elseif {$polar == 0} {
        set coordsys "cart"
    } else {
        puts "polar must be 1 or 0"
        break
    }

    if {[llength $beadname] > 1} {
        set condensed_name [lindex $beadname 0]
        for {set i 1} {$i < [llength $beadname]} {incr i} {
            set addname [lindex $beadname $i]
            set condensed_name "$condensed_name.$addname"
        }
    } elseif {[llength $beadname] == 1} {
        set condensed_name $beadname
    } else {
        puts "beadname must contain a bead name"
        break
    }

    #outfiles setup
    set outfiles [create_outfiles $quantity_of_interest $system $condensed_name $species $acyl_names $coordsys]

    puts "Setup complete. Starting analysis now."	

    set selections [create_atomselections $quantity_of_interest $system $beadname $species $acyl_names $coordsys]

    set heads [atomselect top "name $beadname"]
    
    if {$separate_beads == 0} {
        set interface [atomselect top "(user 1 and within 6 of user 2) or (user 2 and within 6 of user 1)"]
        set blist [list $heads $interface]
        set t1tiltsel [atomselect top "resname $species and name GL1 $tail_one"]
        set t2tiltsel [atomselect top "resname $species and name GL2 $tail_two"]
        set sellist [list $heads $interface $t1tiltsel $t2tiltsel]

        for {set beadpair 0} {$beadpair < [llength $tail_list]} {incr beadpair} {
            set newsel [atomselect top "resname $species and name [lindex $tail_list $beadpair]"]
            lappend sellist $newsel
        }
    } elseif {$separate_beads == 1} {
        set sellist [list $heads]
        set blist [list $heads]
    }  

    array set density_up [initialize_array $N1 $N2 0.0]
    array set density_down [initialize_array $N1 $N2 0.0]
    array set density_zzero [initialize_array $N1 $N2 0.0]
    
    ;#start frame looping here
    for {set frm $start} {$frm < $nframes} {incr frm $step} {
        
        if {$step == 1} {
            leaflet_check $frm $species "PO4" $tailnames
        } elseif {$step > 1} {
            if {$frm != $start} {
                for {set frame [expr $step - 1]} {$frame >= 0} {set frame [expr $frame - 1]} {
                    leaflet_check [expr $frm - $frame] $species "PO4" $tailnames
                }
            }
        }

        foreach selex $sellist {
            $selex frame $frm 
            $selex update
        }
        
        puts "$system $frm"

        set box_height [molinfo top get c]
        
        set box_area_per_frame [expr [molinfo top get a] * [molinfo top get b]]
        lappend boxarea $box_area_per_frame

        set meas_z_zero 0

        set taillength [expr [llength $tail_one] + 1]

        if {$separate_beads == 0} {
            set tilts [tilt_angles $taillength $t1tiltsel $t2tiltsel]
        }

        ;# identify lowest/highest beads for thickness calculations later
        set thickness_lists [prep_thickness_lists $sellist $ref_height]
        set mins_list [lindex $thickness_lists 0]
        set maxs_list [lindex $thickness_lists 1]

        ;# slow implementation of order params
        ;#set order_list [measure_order $sellist $ref_height]

        foreach bead $blist {

            set x_vals [$bead get x]
            set y_vals [$bead get y]
            set z_vals [vecexpr [$bead get z] $ref_height sub]
            set leaflet [$bead get user]

            set bins [bin_generator $x_vals $y_vals $d1 $d2 $dthetadeg $polar]
            set dim1_bins [lindex $bins 0]
            set dim2_bins [lindex $bins 1]

            ;#initialize arrays to zeros
            if {$meas_z_zero == 0} {
                array set totals_up [initialize_array $N1 $N2 0.0]
                array set counts_up [initialize_array $N1 $N2 0.0]
                array set totals_down [initialize_array $N1 $N2 0.0]
                array set counts_down [initialize_array $N1 $N2 0.0]
                array set totals_zplus [initialize_array $N1 $N2 0.0]
                array set thickness_up [initialize_array $N1 $N2 0.0]
                array set thickness_down [initialize_array $N1 $N2 0.0]
                if {$separate_beads == 0} {
                    array set tilts_up [initialize_array $N1 $N2 {0.0 0.0 0.0}]
                    array set tilts_down [initialize_array $N1 $N2 {0.0 0.0 0.0}]
                    ;#array set chain_order_up [initialize_array $N1 $N2 0.0]
                    ;#array set chain_order_down [initialize_array $N1 $N2 0.0]
                }
            } elseif {$meas_z_zero == 1} {
                array set totals_zzero [initialize_array $N1 $N2 0.0]
                array set counts_zzero [initialize_array $N1 $N2 0.0]
            }

            ;#fill in total/count arrays with z sum and count sum
            for {set i 0} {$i < [llength $x_vals]} {incr i} {
                set m [lindex $dim1_bins $i]
                set n [lindex $dim2_bins $i]
                if {$m <= $N1 && $n <= $N2} {
                    if {$meas_z_zero == 0} {
                        if {[lindex $leaflet $i] == 1} {
                            set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                            set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                            set density_up($m,$n) [expr {$density_up($m,$n) + 1}]
                            set thickness_up($m,$n) [expr {$thickness_up($m,$n) + [lindex $mins_list $i]}]
                            if {$separate_beads == 0} {
                                set tilts_up($m,$n) [vecexpr $tilts_up($m,$n) [lindex $tilts $i] add]
                                ;#set chain_order_up($m,$n) [expr {$chain_order_up($m,$n) + [lindex $order $i]}]
                            }
                        } elseif {[lindex $leaflet $i] == 2} {
                            set totals_down($m,$n) [expr {$totals_down($m,$n) + [lindex $z_vals $i]}]
                            set counts_down($m,$n) [expr {$counts_down($m,$n) + 1}]
                            set density_down($m,$n) [expr {$density_down($m,$n) + 1}]
                            set thickness_down($m,$n) [expr {$thickness_down($m,$n) + [lindex $maxs_list $i]}]
                            if {$separate_beads == 0} {
                                set tilts_down($m,$n) [vecexpr $tilts_down($m,$n) [lindex $tilts $i] add]
                                ;#set chain_order_down($m,$n) [expr {$chain_order_down($m,$n) + [lindex $order $i]}]
                            }
                        }
                    } elseif {$meas_z_zero == 1} {
                        set totals_zzero($m,$n) [expr {$totals_zzero($m,$n) + [lindex $z_vals $i]}]
                        set counts_zzero($m,$n) [expr {$counts_zzero($m,$n) + 1}]
                    }
                }
            }    

            ;#turn the z sum into a z avg
            for {set m 0} {$m <= $N1} {incr m} {
                for {set n 0} {$n <= $N2} {incr n} {
                    if {$meas_z_zero == 0} {
                        if {$counts_up($m,$n) != 0.0} {
                            set totals_up($m,$n) [expr $totals_up($m,$n) / $counts_up($m,$n)]
                            set thickness_up($m,$n) [expr $thickness_up($m,$n) / $counts_up($m,$n)]
                            set thickness_up($m,$n) [expr $totals_up($m,$n) - $thickness_up($m,$n)]
                            if {$separate_beads == 0} {
                                set tilts_up($m,$n) [vecexpr $tilts_up($m,$n) $counts_up($m,$n) div]
                                set tilts_up($m,$n) [vecnorm $tilts_up($m,$n)]
                                ;# ADD ORDER HERE
                            }
                        } else {
                            set totals_up($m,$n) "nan"
                            set thickness_up($m,$n) "nan"
                            if {$separate_beads == 0} {
                                set tilts_up($m,$n) "nan nan nan"
                                ;# ADD ORDER HERE
                            }
                        }
                        if {$counts_down($m,$n) != 0.0} {
                            set totals_down($m,$n) [expr $totals_down($m,$n) / $counts_down($m,$n)]
                            set thickness_down($m,$n) [expr $thickness_down($m,$n) / $counts_down($m,$n)]
                            set thickness_down($m,$n) [expr $thickness_down($m,$n) - $totals_down($m,$n)]
                            if {$separate_beads == 0} {
                                set tilts_down($m,$n) [vecexpr $tilts_down($m,$n) $counts_down($m,$n) div]
                                set tilts_down($m,$n) [vecnorm $tilts_down($m,$n)]
                                ;# ADD ORDER HERE
                            }
                        } else {
                            set totals_down($m,$n) "nan"
                            set thickness_down($m,$n) "nan"
                            if {$separate_beads == 0} {
                                set tilts_down($m,$n) "nan nan nan"
                                ;# ADD ORDER HERE
                            }
                        }
                        if {$counts_up($m,$n) != 0.0 && $counts_down($m,$n) != 0.0} {
                            set totals_zplus($m,$n) [expr [expr $totals_up($m,$n) + $counts_down($m,$n)] / 2.0]
                        } else {
                            set totals_zplus($m,$n) "nan"
                        }
                        


                    } elseif {$meas_z_zero == 1} {
                        if {$counts_zzero($m,$n) != 0} {
                            set totals_zzero($m,$n) [expr $totals_zzero($m,$n) / $counts_zzero($m,$n)]
                        } else {
                            set totals_zzero($m,$n) "nan"
                        }
                    }
                }
            }

            ;#output heights to files
            if { $meas_z_zero == 0 } {
                print_frame $N1 $heights_up $d1 $min $N2 [array get totals_up] $polar 
                print_frame $N1 $heights_down $d1 $min $N2 [array get totals_down] $polar 
                print_frame $N1 $heights_zplus $d1 $min $N2 [array get totals_zplus] $polar 
                print_frame $N1 $tilt_up $d1 $min $N2 [array get tilts_up] $polar 
                print_frame $N1 $tilt_down $d1 $min $N2 [array get tilts_down] $polar 
                print_frame $N1 $thick_up $d1 $min $N2 [array get thickness_up] $polar 
                print_frame $N1 $thick_down $d1 $min $N2 [array get thickness_down] $polar 
            } elseif {$meas_z_zero == 1} {
                print_frame $N1 $heights_zzero $d1 $min $N2 [array get totals_up] $polar 
            }

            if {$separate_beads == 0} {
                set meas_z_zero 1
            }
        }
    }

    ;# calculate density
    set delta_frame [expr ($nframes - $start) / $step]
    array set density_up [calculate_density [array get density_up] $min $d1 $d2 $N1 $N2 $delta_frame $polar]
    array set density_down [calculate_density [array get density_down] $min $d1 $d2 $N1 $N2 $delta_frame $polar]
    array set density_zzero [calculate_density [array get density_zzero] $min $d1 $d2 $N1 $N2 $delta_frame $polar]

    ;# prepare to normalize density
    set avg_area [vecexpr $boxarea mean]
    
    ;# x_B * N_L = N_B
    set lipidnum []
    foreach leaflet "1 2" {
        set lipidsel [atomselect top "resname $species and user $leaflet" frame [expr $nframes -1]]
        set lipidres [lsort -unique [$lipidsel get resid]]
        lappend lipidnum [llength $lipidres]
        $lipidsel delete
    }
    
    set up_normfactor [expr [lindex $lipidnum 0] / $avg_area]
    set down_normfactor [expr [lindex $lipidnum 1] / $avg_area]
    set combined_normfactor [expr [vecexpr $lipidnum sum] / [expr 2 * $avg_area]]

    ;# calculate normalized density
    array set density_up [normalize_density [array get density_up] $N1 $N2 $up_normfactor [llength $headnames]]
    array set density_down [normalize_density [array get density_down] $N1 $N2 $down_normfactor [llength $headnames]]
    array set density_zzero [normalize_density [array get density_zzero] $N1 $N2 $combined_normfactor [expr [llength $tailnames] * 2.0]]  

    ;# output densities to files
    print_frame $N1 $dens_up $d1 $min $N2 [array get density_up] $polar 
    print_frame $N1 $dens_down $d1 $min $N2 [array get density_down] $polar 
    print_frame $N1 $dens_zzero $d1 $min $N2 [array get density_zzero] $polar 

    ;#clean up
    foreach channel [file channels "file*"] {
        close $channel
    }
    
    foreach selection [atomselect list] {
        $selection delete
    }
}

proc run_field_mult {list_of_systems polar} {
    foreach item $list_of_systems {
        set gro "/u1/home/js2746/Bending/PC/${item}/${item}.gro"
        set xtc "/u1/home/js2746/Bending/PC/${item}/${item}.xtc"
        #set gro "/u1/home/js2746/Bending/Jam_test/nougattest/${item}/insane.gro"
        #set xtc "/u1/home/js2746/Bending/Jam_test/nougattest/${item}/md_reduced.xtc"
        
        #set gro "/home/jesse/Bending/sims/PG/${item}.gro"
        #set xtc "/home/jesse/Bending/sims/PG/${item}.xtc"
        mol new $gro
        mol addfile $xtc waitfor all
        puts $gro
        puts $xtc
        animate delete beg 0 end 0 skip 0 top
        if {$polar == 1} {
            nougatByField $item 12 30 200 -1 1 1
        } elseif {$polar == 0} {
            nougatByField $item 12 30 200 -1 1 0
        }
        mol delete top
    }
}

proc run_prep {list_of_systems} {
    foreach item $list_of_systems {
        set gro "/u1/home/js2746/Bending/PC/${item}/${item}.gro"
        set xtc "/u1/home/js2746/Bending/PC/${item}/${item}.xtc"
        mol new $gro
        mol addfile $xtc waitfor all
        puts $gro
        puts $xtc
        animate delete beg 0 end 0 skip 0 top
        cell_prep $item
        mol delete top
    }
}
