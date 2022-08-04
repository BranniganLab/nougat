
package require pbctools

# EDIT THE PATHS HERE
# TELL nougat WHERE TO FIND YOUR VERSIONS OF qwrap AND vecexpr
#set QWRAP "~/qwrap-master"
#set VEC "~/utilities/vecexpr"

set UTILS "~/PolarHeightBinning/utilities"

load ${UTILS}/qwrap.so
load ${UTILS}/vecexpr.so
#load ~/qwrap/qwrap.so 
#load ~/vecexpr/vecexpr.so 

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
    set excluded_sel "resname W ION 'CL-' 'NA+'"

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
    set full_tails []
    foreach lipidtype $acyl_names {
        foreach tail $lipidtype {
            foreach bead $tail {
                lappend full_tails $bead
            }
        }
    }
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
    lappend return_list $full_tails
    lappend return_list $reference_point
    return $return_list  
}

proc RtoD {r} {
    global M_PI
    return [expr $r*180.0/$M_PI]
}

#  returns a least squares fit for each lipid tail in the system at once
#  Fit the points x to x = ai + b, i=0...N-1, and return the value of a 
# a = 12/( (N(N^2 - 1)) ) sum[ (i-(N-1)/2) * xi]
# reference: Bevington
proc lsq_vec { X } {
  #initialize
  set len [llength $X]
  set diff [expr $len-1]
  set d [expr 0.5*$diff] ;#normalization factor
  set I {}

  #Make a list from 0 to len
  for {set k 0} {$k < $len} {incr k} {
    lappend I $k
  }
  
  #Get element-wise diff between I(integers 0 to len) and d(0.5*(len-1)). stack=1vector
  #push X to stack. stack=2vecs
  #Multiply X and the differences. stack=1vec
  #sum over the vector. stack=1scalar
  #store the scalar in tot
  vecexpr $I $d sub <X mult sum >tot
  
  return $tot
}

proc fit_all_tails { tail_length list_of_tail_coords } {
  set vectors []
  set N_lipids [expr [llength $list_of_tail_coords] / $tail_length]
  set start_idx 0
  for {set i 0} {$i < $N_lipids} {incr i} {
    set next_idx [expr $start_idx + $tail_length]
    set end_idx [expr $next_idx - 1]
    set coords [lrange $list_of_tail_coords $start_idx $end_idx]
    lappend vectors [lsq_vec $coords]
    set start_idx $next_idx
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
proc print_frame {N1 outfiles key d1 min N2 polar} {

    set file [dict get $outfiles $key fname]

    if {$polar == 1} {
        set N2 [expr $N2-1]
    }

    ;# starts new line in outfile with bin values
    for {set m 0} {$m <= $N1} {incr m} {
        print_line_init $file $m $d1 $min
        ;# adds bin values through penultimate value in one line
        for {set n 0} {$n < $N2} {incr n} {
            if {[dict exists $outfiles $key bin "$m,$n"]} {
                print_value $file [dict get $outfiles $key bin "$m,$n"] 0
            } else {
                print_value $file "nan" 0
            }
        }
        ;# adds final value and starts new line in outfile
        if {[dict exists $outfiles $key bin "$m,$N2"]} {
            print_value $file [dict get $outfiles $key bin "$m,$N2"] 1
        } else {
            print_value $file "nan" 1
        }
    }
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
    
    ;# change user3 to match tail number
    ;# makes tails separable for tilt/order analysis
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

    return $taillist
}

proc tilt_angles {length xvals yvals zvals} {
    set tilt_list []
    set xvec [fit_all_tails $length $xvals]
    set yvec [fit_all_tails $length $yvals]
    set zvec [fit_all_tails $length $zvals]
    for {set i 0} {$i < [llength $xvec]} {incr i} {
        set vector "[lindex $xvec $i] [lindex $yvec $i] [lindex $zvec $i]"
        set norm [vecnorm $vector]
        lappend tilt_list $norm
    }
    set final_tilt_list []
    foreach item $tilt_list {
        set final_tilt_list [concat $final_tilt_list [lrepeat $length $item]]
    }
    return [list $final_tilt_list]
}

proc bin_assigner {x_vals y_vals d1 d2 dthetadeg polar} {
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

proc create_outfiles {system quantity_of_interest headnames species taillist coordsys} {

    if {$quantity_of_interest eq "height_density"} {
        dict set outfiles heights_up fname [open "${system}.zone.${headnames}.${coordsys}.height.dat" w]
        dict set outfiles heights_down fname [open "${system}.ztwo.${headnames}.${coordsys}.height.dat" w]
        foreach lipidtype $species {
            dict set outfiles density_up_${lipidtype} fname [open "${system}.${lipidtype}.zone.${coordsys}.density.dat" w]
            dict set outfiles density_down_${lipidtype} fname [open "${system}.${lipidtype}.ztwo.${coordsys}.density.dat" w]
        }
    } elseif {$quantity_of_interest eq "tilt_order"} {
        for {set i 0} {$i < [llength $taillist]} {incr i} {
            set lipidtype [lindex $species $i]
            for {set j 0} {$j < [llength [lindex $taillist $i]]} {incr j} {
                set tailnum "tail$j"
                dict set outfiles tilts_up_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.zone.${coordsys}.tilt.dat" w]
                dict set outfiles tilts_down_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.ztwo.${coordsys}.tilt.dat" w]
                dict set outfiles order_up_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.zone.${coordsys}.order.dat" w]
                dict set outfiles order_down_${lipidtype}_${tailnum} fname [open "${system}.${lipidtype}.${tailnum}.ztwo.${coordsys}.order.dat" w]
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

proc concat_names { headnames } {
    ;# concatenate all beadnames together for file naming purposes
    if {[llength $headnames] > 1} {
        set condensed_name [lindex $headnames 0]
        for {set i 1} {$i < [llength $headnames]} {incr i} {
            set addname [lindex $headnames $i]
            set condensed_name "$condensed_name.$addname"
        }
    } elseif {[llength $headnames] == 1} {
        set condensed_name $beadname
    } else {
        puts "headnames must contain a bead name"
        break
    }
    return $condensed_name
}

proc create_res_dict { species headnames lipid_list name_list resid_list dim1_bins_list dim2_bins_list leaflet_list} {
    dict set res_dict dummy "dummy"
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
    dict unset res_dict dummy
    return $res_dict
}

proc do_height_density_binning {res_dict outfiles leaflet_list lipid_list zvals_list} {
    dict for {bin indices} $res_dict {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr [string length $bin] - 3]]
        set newlist []
        foreach indx $indices {
            if {$leaf == 1} {
                set dens_key "density_up_[lindex $lipid_list $indx]"
                set height_key "heights_up"
            } else {
                set dens_key "density_down_[lindex $lipid_list $indx]"
                set height_key "heights_down"
            }
            lappend newlist [lindex $zvals_list $indx]
        }
        if {[llength $newlist] == 1} {
            dict set outfiles $height_key bin $correct_bin [lindex $newlist 0]
            dict set outfiles $dens_key bin $correct_bin 1
        } elseif {[llength $newlist] > 1} {
            set avg [vecexpr $newlist mean]
            dict set outfiles $height_key bin $correct_bin $avg 
            dict set outfiles $dens_key bin $correct_bin [llength $newlist]
        }
    }
    return $outfiles
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

proc order_params {length xvals yvals zvals} {
    set order_list []
    set temp_list []
    for {set i 1} {$i <= [llength $xvals]} {incr i} {
        if {[expr $i % $length] == 0} {
            set avg [vecexpr $temp_list mean]
            set step2 [expr $avg*3.0 - 1]
            set order [expr $step2 / 2.0]
            lappend order_list $order
            set temp_list []
            continue
        } else {
            set start [list [lindex $xvals $i] [lindex $yvals $i] [lindex $zvals $i]]
            set end [list [lindex $xvals [expr $i-1]] [lindex $yvals [expr $i-1]] [lindex $zvals [expr $i-1]]]
            set magn_a [vecdist $start $end]
            set theta [expr [expr [lindex $zvals $i] - [lindex $zvals $i-1]] / $magn_a]
            if {($theta > 1) || ($theta < -1)} {
                puts "Something is wrong with your order params"
                puts "This is out of the range allowed for arccos"
                return
            }
            lappend temp_list [expr $theta * $theta]
        }
    }
    set final_order_list []
    foreach item $order_list {
        set final_order_list [concat $final_order_list [lrepeat $length $item]]
    }
    return [list $final_order_list]
}

proc do_tilt_order_binning {res_dict outfiles leaflet_list lipid_list tilts orders tail_list} {
    dict for {bin indices} $res_dict {
        set leaf [string range $bin end end]
        set correct_bin [string range $bin 0 [expr [string length $bin] - 3]]
        set tiltlist []
        set orderlist []
        foreach indx $indices {
            if {$leaf == 1} {
                set tilt_key "tilts_up_[lindex $lipid_list $indx]_tail[expr int([expr [lindex $tail_list $indx] - 1])]"
                set order_key "order_up_[lindex $lipid_list $indx]_tail[expr int([expr [lindex $tail_list $indx] - 1])]"
            } else {
                set tilt_key "tilts_down_[lindex $lipid_list $indx]_tail[expr int([expr [lindex $tail_list $indx] - 1])]"
                set order_key "order_down_[lindex $lipid_list $indx]_tail[expr int([expr [lindex $tail_list $indx] - 1])]"
            }
            lappend tiltlist [lindex [lindex $tilts 0] $indx]
            lappend orderlist [lindex [lindex $orders 0] $indx]
        }
        if {[llength $tiltlist] == 1} {
            dict set outfiles $tilt_key bin $correct_bin [lindex $tiltlist 0]
            dict set outfiles $order_key bin $correct_bin [lindex $orderlist 0]
        } elseif {[llength $tiltlist] > 1} {
            set tiltsum [list 0 0 0]
            for {set i 0} {$i < [llength $tiltlist]} {incr i} {
                set tiltsum [vecexpr [lindex $tiltlist $i] $tiltsum add]
            }
            set tiltavg [vecexpr $tiltsum [llength $tiltlist] div]
            dict set outfiles $tilt_key bin $correct_bin $tiltavg
            set orderavg [vecexpr $orderlist mean] 
            dict set outfiles $order_key bin $correct_bin $orderavg
        }
    }
    return $outfiles
}

;########################################################################################
;# polarHeight Functions

proc start_nougat {system d1 N2 start end step polar} {

    ;# running cell_prep will do some important initial configuration based on user input. 
    ;# check the extensive documentation at the top of this file for instructions.
    set important_variables [cell_prep $system $start]
    
    ;# unpack user-provided info from cell_prep
    set species [lindex $important_variables 0]
    set headnames [lindex $important_variables 1]
    set acyl_names [lindex $important_variables 2]
    set full_tails [lindex $important_variables 3]
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

    ;# change this value if you want to exclude an inner radius in polar coords
    set min 0 

    ;# determine number and size of bins
    set bindims [bin_prep $nframes $polar $min $d1]

    ;# add all these new values to important_values for easy transfer
    lappend important_variables $start 
    lappend important_variables $nframes 
    lappend important_variables $step 
    lappend important_variables $ref_height
    lappend important_variables $min

    ;# run nougat twice, once to compute height and density and once to compute
    ;# lipid tail vectors and order parameters
    run_nougat $system $important_variables $bindims $polar "height_density" 
    run_nougat $system $important_variables $bindims $polar "tilt_order" 

}

proc run_nougat {system important_variables bindims polar quantity_of_interest} {  

    ;# unpack important_variables
    set species [lindex $important_variables 0]
    set headnames [lindex $important_variables 1]
    set acyl_names [lindex $important_variables 2]
    set full_tails [lindex $important_variables 3]
    set start [lindex $important_variables 5]
    set nframes [lindex $important_variables 6]
    set step [lindex $important_variables 7]
    set ref_height [lindex $important_variables 8]
    set min [lindex $important_variables 9]
    
    ;# unpack bindims
    set d1 [lindex $bindims 0]
    set d2 [lindex $bindims 1]
    set N1 [lindex $bindims 2]
    set N2 [lindex $bindims 3]
    set dthetadeg [lindex $bindims 4]

    set boxarea []

    ;# generate string for polar or cartesian coordinates
    if {$polar == 1} {
        set coordsys "polar"
    } elseif {$polar == 0} {
        set coordsys "cart"
    } else {
        puts "polar must be 1 or 0"
        break
    }

    ;# outfiles setup as dict
    set outfiles [create_outfiles $system $quantity_of_interest [concat_names $headnames] $species $acyl_names $coordsys]

    ;#atomselections setup as dict
    if {$quantity_of_interest eq "height_density"} {
        dict set selections z1z2 [atomselect top "resname $species and name $full_tails"]
        #dict set selections z0 [atomselect top "(user 1 and within 6 of user 2) or (user 2 and within 6 of user 1)"]
    } elseif {$quantity_of_interest eq "tilt_order"} {
        set lists [tail_length_sorter $species $acyl_names]
        set sellist [lindex $lists 0]
        set lenlist [lindex $lists 1]
        foreach sel $sellist len $lenlist {
            dict set selections $len [atomselect top "$sel"] 
        }
    }


    puts "Setup complete. Starting frame analysis now."   

    ;# start frame looping here
    for {set frm $start} {$frm < $nframes} {incr frm $step} {
        
        ;# update leaflets in case lipids have flip-flopped
        if {$step == 1} {
            leaflet_check $frm $species $acyl_names
        } elseif {$step > 1} {
            if {$frm != $start} {
                for {set frame [expr $step - 1]} {$frame >= 0} {set frame [expr $frame - 1]} {
                    leaflet_check [expr $frm - $frame] $species $acyl_names
                }
            }
        }

        puts "$system $quantity_of_interest $frm"

        ;# calculate avg box area for density normalization later
        set box_height [molinfo top get c]
        set box_area_per_frame [expr [molinfo top get a] * [molinfo top get b]]
        lappend boxarea $box_area_per_frame
        
        ;# height_density only has one selection, so this will execute once.
        ;# tilt_order has different selections, one for each tail length present
        ;# in the system, so this will execute as many times as there are
        ;# unique tail lengths.
        foreach selex [dict keys $selections] {
            
            ;# $selex is a dict key that holds an atomselection as its value
            set sel [dict get $selections $selex]
            
            $sel frame $frm 
            $sel update

            set xvals_list [$sel get x]
            set yvals_list [$sel get y]
            set resid_list [$sel get resid]
            set lipid_list [$sel get resname]
            set name_list [$sel get name]
            
            ;# the z vals are subtracted by a reference height provided in cell_prep 
            set zvals_list [vecexpr [$sel get z] $ref_height sub]
            
            ;# user contains a 1 or 2 for outer or inner leaflet, respectively
            set leaflet_list [$sel get user]

            ;# user3 contains an int that describes which tail in the lipid this is
            ;# E.G. POPC will have 0 or 1 (it has two tails)
            ;# E.G. OANT will have 0, 1, 2, 3, 4, or 5 (it has 6 tails)
            set tail_list [$sel get user3]

            ;# calculate which bins each bead belongs in along both axes
            ;# and return as two lists of same length as the lists above
            set bins [bin_assigner $xvals_list $yvals_list $d1 $d2 $dthetadeg $polar]
            set dim1_bins_list [lindex $bins 0]
            set dim2_bins_list [lindex $bins 1]

            ;# Binning is controlled by the bead designated in $headnames.
            ;# Creates a dict that contains the bin and leaflet information linked to
            ;# a resid and index number. Facilitates easy binning later. 
            set res_dict [create_res_dict $species $headnames $lipid_list $name_list $resid_list $dim1_bins_list $dim2_bins_list $leaflet_list]
            
            ;# Make necessary calculations (in case of tilts/orders) and then bin them
            if {$quantity_of_interest eq "height_density"} {
                set outfiles [do_height_density_binning $res_dict $outfiles $leaflet_list $lipid_list $zvals_list]
            } elseif {$quantity_of_interest eq "tilt_order"} {
                set tilts [tilt_angles [dict keys $selections] $xvals_list $yvals_list $zvals_list]
                set orders [order_params [dict keys $selections] $xvals_list $yvals_list $zvals_list]
                set outfiles [do_tilt_order_binning $res_dict $outfiles $leaflet_list $lipid_list $tilts $orders $tail_list]
            }

            ;# Now that all information has been binned, print it to files
            foreach key [dict keys $outfiles] {
                print_frame $N1 $outfiles $key $d1 $min $N2 $polar

                ;# precautionary cleanup before next step
                set outfiles [dict unset outfiles $key bin]
            } 

            ;# precautionary cleanup before next step
            dict remove res_dict
        }
    }

    ;# close all outfiles
    foreach channel [file channels "file*"] {
        close $channel
    }

    ;# delete all atomselections in scope
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