
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
    ;# Only works for systems with 100% one type of lipid
    ;# will be fixed in June update
    ;# no edits required
    set acyl_names [tail_analyzer $species]
    set tail_one [lindex $acyl_names 0]
    set tail_two [lindex $acyl_names 1]
    set t1T [lindex $tail_one end]
    set t2T [lindex $tail_two end]
    set tailnames "$t1T $t2T"
    set tail_list []
    if {[llength $tail_one] == [llength $tail_two]} {
        for {set i 1} {$i < [llength $tail_one]} {incr i} {
                set t1bead [lindex $tail_one $i]
                set t2bead [lindex $tail_two $i]
                set names "$t1bead $t2bead"
                lappend tail_list $names
        }
    } else {
        puts "Tail lengths are different. Teach me what to do."
    }

    ;# Assigns lipids to user value 1 or 2 depending on leaflet
    ;# no edits required
    leaflet_sorter $species $tailnames $lastframe 

    ;# this will only work if your TMD helices are set to occupancy 1
    ;# otherwise, comment it out
    ;# all it does is put a dot on the polar heatmap where a TMD helix should be, so not essential at all
    Protein_Position $system $headnames $tailnames 

    ;#**********************************************************
    ;#          MAKE EDITS ABOVE BEFORE STARTING
    ;#********************************************************** 

    set return_list [] 
    lappend return_list $species 
    lappend return_list $headnames 
    lappend return_list $tailnames 
    lappend return_list $tail_one 
    lappend return_list $tail_two 
    lappend return_list $tail_list 
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
# just relies on whether the tails are above/below the heads
# to sort leaflets
proc leaflet_check {frm species headname tailname} {
    set headsel [atomselect top "resname $species and name $headname" frame $frm]
    set head_ht [$headsel get z]
    set resids [$headsel get resid]
    set chains [$headsel get chain]
    $headsel delete
    set tail_hts [lrepeat [llength $head_ht] 0]
    foreach bead $tailname {
        set tailsel [atomselect top "resname $species and name $bead" frame $frm]
        set tail_ht [$tailsel get z]
        set tail_hts [vecexpr $tail_hts $tail_ht add]
        $tailsel delete
    }
    set avg_tail_ht [vecexpr $tail_hts [llength $tailname] div]

    set test [vecexpr $head_ht $avg_tail_ht sub]
    for {set i 0} {$i < [llength $test]} {incr i} {
        if {[expr abs([lindex $test $i])] > 40} {
            set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
            $sel set chain "U"
            $sel delete
        } else {
            if {[lindex $test $i] > 0 && ([lindex $chains $i] == "L" || [lindex $chains $i] == "Z")} {
                set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
                $sel set chain "U"
                $sel delete
            } elseif {[lindex $test $i] < 0 && ([lindex $chains $i] == "U" || [lindex $chains $i] == "Z")} {
                set sel [atomselect top "resname $species and resid [lindex $resids $i]" frame $frm]
                $sel set chain "L"
                $sel delete
            }
        }
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


proc leaflet_sorter {species tailnames sample_frame} {
    puts "Starting leaflet sorting"
    set nframes [molinfo top get numframes]
    set sel [atomselect top "name PO4" frame 0]
    set resids [$sel get resid]
    set indexs [$sel get index]
    $sel delete
    foreach resd $resids indx $indexs {
        set lipid [atomselect top "resname $species and resid $resd" frame 0]
        set po4 [atomselect top "index $indx" frame 0]
        set lipid_com [measure center $lipid weight mass]
        set lipid_com_z [lindex $lipid_com 2]
        set po4_z [$po4 get z]
        if {$po4_z > $lipid_com_z} {
            $lipid set chain U
            $lipid set user 1
        } else {
            $lipid set chain L
            $lipid set user 2
        }
        $lipid delete
        $po4 delete
    }

    for {set frm 0} {$frm <= $sample_frame} {incr frm} {
        leaflet_check $frm $species "PO4" $tailnames
    }
    
    puts "Leaflet sorting complete!"
}

proc tail_analyzer { species } {

    set sel [atomselect top "resname $species"]
    set res [$sel get resid]
    $sel delete
    set sel [atomselect top "resname $species and resid [lindex $res 0]"]
    set names [$sel get name]
    $sel delete
 
    set tail_one []
    set tail_two []

    foreach nm $names {
        if {[string match ??A $nm]} {
            lappend tail_one $nm
        } elseif {[string match ??B $nm]} {
            lappend tail_two $nm
        }
    }  

    return [list $tail_one $tail_two]
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

;# NEED TO UPDATE

;########################################################################################
;# polarHeight Functions

proc nougatByField {system d1 N2 start end step polar separate_beads} {

    set important_variables [cell_prep $system $start]
    set species [lindex $important_variables 0]
    set headnames [lindex $important_variables 1]
    set tailnames [lindex $important_variables 2]
    set tail_one [lindex $important_variables 3]
    set tail_two [lindex $important_variables 4]
    set tail_list [lindex $important_variables 5]
    set reference_point [lindex $important_variables 6]

    #need to calculate heights relative to some point on the inclusion:
    set ref_bead [atomselect top "$reference_point"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [vecexpr $ref_height mean]

    ;# set nframes based on $end input
    if {$end == -1} {
        set nframes [molinfo top get numframes]
    } else {
        set nframes $end
    }

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
        set N1 [expr [expr $range1 / $d1] - 1.0] 
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
            set N2 [expr [expr $range2 / $d2] - 1.0] 
        } else {
            set N2 [expr $range2 / $d2]
        }
    }

    if {$polar == 1} {
        set coordsys "polar"
    } elseif {$polar == 0} {
        set coordsys "cart"
    }

    lappend important_variables $d1 
    lappend important_variables $d2 
    lappend important_variables $N1
    lappend important_variables $N2
    lappend important_variables $start 
    lappend important_variables $nframes 
    lappend important_variables $step 
    lappend important_variables $ref_height
    lappend important_variables $min 

    run_nougat $system $headnames $coordsys $important_variables $polar 0
    if {$separate_beads == 1} {
        foreach beadpair $tail_list {
            run_nougat $system $beadpair $coordsys $important_variables $polar 1
        }
    }
}

proc run_nougat {system beadname coordsys important_variables polar separate_beads} {  

    set boxarea []

    set species [lindex $important_variables 0]
    set headnames [lindex $important_variables 1]
    set tailnames [lindex $important_variables 2]
    set tail_one [lindex $important_variables 3]
    set tail_two [lindex $important_variables 4]
    set tail_list [lindex $important_variables 5]
    set reference_point [lindex $important_variables 6]
    set d1 [lindex $important_variables 7]
    set d2 [lindex $important_variables 8]
    set N1 [lindex $important_variables 9]
    set N2 [lindex $important_variables 10]
    set start [lindex $important_variables 11]
    set nframes [lindex $important_variables 12]
    set step [lindex $important_variables 13]
    set ref_height [lindex $important_variables 14]
    set min [lindex $important_variables 15]

    set name1 [lindex $beadname 0]
    set name2 [lindex $beadname 1]
    set condensed_name "${name1}.${name2}"

    #outfiles setup
    set heights_up [open "${system}.zone.${condensed_name}.${coordsys}.height.dat" w]
    set heights_down [open "${system}.ztwo.${condensed_name}.${coordsys}.height.dat" w]
    set heights_zplus [open "${system}.zplus.${condensed_name}.${coordsys}.height.dat" w]
    set dens_up [open "${system}.zone.${condensed_name}.${coordsys}.density.dat" w]
    set dens_down [open "${system}.ztwo.${condensed_name}.${coordsys}.density.dat" w]
    if {$separate_beads == 0} {
        set heights_zzero [open "${system}.zzero.${condensed_name}.${coordsys}.height.dat" w]
        set dens_zzero [open "${system}.zzero.${condensed_name}.${coordsys}.density.dat" w]
        set tilt_up [open "${system}.zone.${condensed_name}.${coordsys}.tilt.dat" w]
        set tilt_down [open "${system}.ztwo.${condensed_name}.${coordsys}.tilt.dat" w]
        ;#set thick_up [open "${system}.zone.${condensed_name}.${coordsys}.thickness.dat" w]
        ;#set thick_down [open "${system}.ztwo.${condensed_name}.${coordsys}.thickness.dat" w]
        ;#set order_up [open "${system}.ztwo.${condensed_name}.${coordsys}.order.dat" w]
        ;#set order_down [open "${system}.ztwo.${condensed_name}.${coordsys}.order.dat" w]
    }

    puts "Setup complete. Starting analysis now."	


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
        
        leaflet_check $frm $species "PO4" $tailnames

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

        foreach bead $blist {

            set x_vals [$bead get x]
            set y_vals [$bead get y]
            set z_vals [vecexpr [$bead get z] $ref_height sub]
            set leaflet [$bead get user]
            set resids [$bead get resid]
            set indexs [$bead get index]

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

            ;#initialize arrays to zeros
            if {$meas_z_zero == 0} {
                array set totals_up [initialize_array $N1 $N2 0.0]
                array set counts_up [initialize_array $N1 $N2 0.0]
                array set totals_down [initialize_array $N1 $N2 0.0]
                array set counts_down [initialize_array $N1 $N2 0.0]
                array set totals_zplus [initialize_array $N1 $N2 0.0]
                if {$separate_beads == 0} {
                    array set tilts_up [initialize_array $N1 $N2 {0.0 0.0 0.0}]
                    array set tilts_down [initialize_array $N1 $N2 {0.0 0.0 0.0}]
                    ;#array set thickness_up [initialize_array $N1 $N2 0.0]
                    ;#array set thickness_down [initialize_array $N1 $N2 0.0]
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
                            if {$separate_beads == 0} {
                                set tilts_up($m,$n) [vecexpr $tilts_up($m,$n) [lindex $tilts $i] add]
                                ;#set thickness_up($m,$n) [expr {$thickness_up($m,$n) + [lindex $thickness $i]}]
                                ;#set chain_order_up($m,$n) [expr {$chain_order_up($m,$n) + [lindex $order $i]}]
                            }
                        } elseif {[lindex $leaflet $i] == 2} {
                            set totals_down($m,$n) [expr {$totals_down($m,$n) + [lindex $z_vals $i]}]
                            set counts_down($m,$n) [expr {$counts_down($m,$n) + 1}]
                            set density_down($m,$n) [expr {$density_down($m,$n) + 1}]
                            if {$separate_beads == 0} {
                                set tilts_down($m,$n) [vecexpr $tilts_down($m,$n) [lindex $tilts $i] add]
                                ;#set thickness_down($m,$n) [expr {$thickness_down($m,$n) + [lindex $thickness $i]}]
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
                        if {$counts_up($m,$n) != 0} {
                            set totals_up($m,$n) [expr $totals_up($m,$n) / [expr 1.0 * $counts_up($m,$n)]]
                            if {$separate_beads == 0} {
                                set tilts_up($m,$n) [vecexpr $tilts_up($m,$n) $counts_up($m,$n) div]
                                set tilts_up($m,$n) [vecnorm $tilts_up($m,$n)]
                                ;# ADD THICKNESS AND ORDER HERE
                            }
                        } else {
                            set totals_up($m,$n) "nan"
                            if {$separate_beads == 0} {
                                set tilts_up($m,$n) "nan nan nan"
                                ;# ADD THICKNESS AND ORDER HERE
                            }
                        }
                        if {$counts_down($m,$n) != 0} {
                            set totals_down($m,$n) [expr $totals_down($m,$n) / [expr 1.0 * $counts_down($m,$n)]]
                            if {$separate_beads == 0} {
                                set tilts_down($m,$n) [vecexpr $tilts_down($m,$n) $counts_down($m,$n) div]
                                set tilts_down($m,$n) [vecnorm $tilts_down($m,$n)]
                                ;# ADD THICKNESS AND ORDER HERE
                            }
                        } else {
                            set totals_down($m,$n) "nan"
                            if {$separate_beads == 0} {
                                set tilts_down($m,$n) "nan nan nan"
                                ;# ADD THICKNESS AND ORDER HERE
                            }
                        }
                        if {$counts_up($m,$n) != 0 && $counts_down($m,$n) != 0} {
                            set totals_zplus($m,$n) [expr [expr $totals_up($m,$n) + $counts_down($m,$n)] / 2.0]
                        } else {
                            set totals_zplus($m,$n) "nan"
                        }
                        


                    } elseif {$meas_z_zero == 1} {
                        if {$counts_zzero($m,$n) != 0} {
                            set totals_zzero($m,$n) [expr $totals_zzero($m,$n) / [expr 1.0 * $counts_zzero($m,$n)]]
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
    close $heights_up
    close $heights_down
    close $heights_zplus
    close $heights_zzero
    close $dens_up
    close $dens_down
    close $dens_zzero
    close $tilt_up 
    close $tilt_down 
    ;#close $thick_up 
    ;#close $thick_down 
    ;#close $order_up 
    ;#close $order_down 
    
    foreach selection [atomselect list] {
        $selection delete
    }
}

proc nougatByBead {system} {
    
    set important_variables [cell_prep $system $start]
    set species [lindex $important_variables 0]
    set headnames [lindex $important_variables 1]
    set tailnames [lindex $important_variables 2]
    set tail_one [lindex $important_variables 3]
    set tail_two [lindex $important_variables 4]
    set tail_list [lindex $important_variables 5]
    set reference_point [lindex $important_variables 6]

    #need to calculate heights relative to some point on the protein:
    set ref_bead [atomselect top "$reference_point"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [vecexpr $ref_height mean]

    puts "Helper scripts complete. Starting analysis now."  

    foreach beadpair $tail_list {
        set name1 [lindex $beadpair 0]
        set name2 [lindex $beadpair 1]
        set condensed_name "${name1}.${name2}"

        #outfiles setup
        set heights_up [open "${system}.${condensed_name}.zone.height.dat" w]
        set heights_down [open "${system}.${condensed_name}.ztwo.height.dat" w]
        set heights_zplus [open "${system}.${condensed_name}.zplus.height.dat" w]
        set dens_up [open "${system}.${condensed_name}.zone.density.dat" w]
        set dens_down [open "${system}.${condensed_name}.ztwo.density.dat" w]
        set dens_zplus [open "${system}.${condensed_name}.zplus.density.dat" w]

        set bead [atomselect top "name $beadpair"]

        array set density_up [initialize_array $Nr $Ntheta 0]
        array set density_down [initialize_array $Nr $Ntheta 0]
        array set density_zplus [initialize_array $Nr $Ntheta 0]

        ;#start frame looping here
        for {set frm $sample_frame} {$frm < $nframes} {incr frm $dt} {
            puts "$system $beadpair $frm"

            leaflet_check $frm $species "PO4" $tailnames
            ;#rmv_outliers $frm $species "PO4" 15

            set box_height [molinfo top get c]
            
            set box_area_per_frame [expr [molinfo top get a] * [molinfo top get a]]
            lappend boxarea $box_area_per_frame

            $bead frame $frm 
            $bead update

            set x_vals [$bead get x] 
            set y_vals [$bead get y]
            set z_vals [vecexpr [$bead get z] $ref_height sub]
            set chains [$bead get user]
            set resids [$bead get resid]
            set indexs [$bead get index]

            ;#get theta values for all x,y pairs
            set theta_vals [vecexpr $y_vals $x_vals atan2 pi div 180 mult]  

            ;#atan2 gives values from -180 to 180; shifting to 0 to 360                            
            for {set i 0} {$i<[llength $theta_vals]} {incr i} {                                         
                if {[lindex $theta_vals $i] < 0} {                                                      
                    set theta_vals [lreplace $theta_vals $i $i [expr [lindex $theta_vals $i]+360]]
                }
            }

            ;#turn into bin numbers rather than theta values
            vecexpr [vecexpr $theta_vals $dthetadeg div] floor &theta_bins
            
            ;#calculate distance from origin for all x,y pairs
            set r_vals [vecexpr [vecexpr [vecexpr $x_vals sq] [vecexpr $y_vals sq] add] sqrt]
            
            ;#turn into bin numbers rather than r values
            vecexpr [vecexpr $r_vals $dr div] floor &r_bins

            ;#initialize arrays to zeros
            array set totals_up [initialize_array $Nr $Ntheta 0]
            array set counts_up [initialize_array $Nr $Ntheta 0]
            array set totals_down [initialize_array $Nr $Ntheta 0]
            array set counts_down [initialize_array $Nr $Ntheta 0]
            array set totals_zplus [initialize_array $Nr $Ntheta 0]
            array set counts_zplus [initialize_array $Nr $Ntheta 0]

            ;#fill in arrays with z sum and count sum
            for {set i 0} {$i < [llength $r_vals]} {incr i} {
                set m [lindex $r_bins $i]
                set n [lindex $theta_bins $i]
                if {[lindex $z_vals $i] > [expr $protein_top + 5]} {
                    set [lindex $z_vals $i] [expr [lindex $z_vals $i] - $box_height]
                }
                if {$m <= $Nr} {
                    if {[lindex $chains $i] == 1} {
                        set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                        set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                        set density_up($m,$n) [expr {$density_up($m,$n) + 1}]
                    } elseif {[lindex $chains $i] == 2} {
                        set totals_down($m,$n) [expr {$totals_down($m,$n) + [lindex $z_vals $i]}]
                        set counts_down($m,$n) [expr {$counts_down($m,$n) + 1}]
                        set density_down($m,$n) [expr {$density_down($m,$n) + 1}]
                    }
                }
            }    

            ;#turn the z sum into a z avg
            for {set m 0} {$m <= $Nr} {incr m} {
                for {set n 0} {$n <= $Ntheta} {incr n} {
                    if {$counts_up($m,$n) != 0} {
                        set totals_up($m,$n) [expr $totals_up($m,$n) / $counts_up($m,$n)]
                    } else {
                        set totals_up($m,$n) "nan"
                    }
                    if {$counts_down($m,$n) != 0} {
                        set totals_down($m,$n) [expr $totals_down($m,$n) / $counts_down($m,$n)]
                        if {$counts_up($m,$n) != 0} {
                            set totals_zplus($m,$n) [expr [expr $totals_up($m,$n) + $totals_down($m,$n)]/2.0]
                            set counts_zplus($m,$n) [expr $counts_up($m,$n) + $counts_down($m,$n)]
                        } else {
                            set totals_zplus($m,$n) "nan"
                        }
                    } else {
                        set totals_down($m,$n) "nan"
                        set totals_zplus($m,$n) "nan"
                    }
                }
            }

            ;#output heights to files
            print_frame $Nr $heights_up $dr $Rmin $Ntheta [array get totals_up] $polar 
            print_frame $Nr $heights_down $dr $Rmin $Ntheta [array get totals_down] $polar 
            print_frame $Nr $heights_zplus $dr $Rmin $Ntheta [array get totals_zplus] $polar 
        }
        ;# calculate density
        set delta_frame [expr ($nframes - $sample_frame) / $dt]
        array set density_up [calculate_density [array get density_up] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]
        array set density_down [calculate_density [array get density_down] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]
        array set density_zplus [calculate_density [array get density_zplus] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]

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
        array set density_up [normalize_density [array get density_up] $Nr $Ntheta $up_normfactor [llength $headnames]]
        array set density_down [normalize_density [array get density_down] $Nr $Ntheta $down_normfactor [llength $headnames]]
        array set density_zplus [normalize_density [array get density_zplus] $Nr $Ntheta $combined_normfactor [expr [llength $headnames] * 2.0]]

        ;# output densities to files
        print_frame $Nr $dens_up $dr $Rmin $Ntheta [array get density_up] $polar 
        print_frame $Nr $dens_down $dr $Rmin $Ntheta [array get density_down] $polar 
        print_frame $Nr $dens_zplus $dr $Rmin $Ntheta [array get density_zplus] $polar 

        ;#clean up
        close $heights_up
        close $heights_down
        close $heights_zplus
        close $dens_up
        close $dens_down
        close $dens_zplus
        $bead delete
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

proc run_bead_mult {list_of_systems} {
    foreach item $list_of_systems {
        set gro "/u1/home/js2746/Bending/PC/${item}/${item}.gro"
        set xtc "/u1/home/js2746/Bending/PC/${item}/${item}.xtc"
        #set gro "/home/jesse/Bending/sims/PG/${item}.gro"
        #set xtc "/home/jesse/Bending/sims/PG/${item}.xtc"
        mol new $gro
        mol addfile $xtc waitfor all
        puts $gro
        puts $xtc
        animate delete beg 0 end 0 skip 0 top
        polarHeightByBead $item
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
