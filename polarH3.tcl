
package require pbctools

#locations on JS home computer
#set UTILS "~/Bending/scripts/PolarHeightBinning/utilities" 
#set QWRAP "~/qwrap"
#set VEC "~/vecexpr"

#locations on Belenus
set UTILS "~/PolarHeightBinning/utilities"
set QWRAP "~/qwrap-master"
set VEC "~/utilities/vecexpr"

source $UTILS/assign_tilts.tcl
source $UTILS/leaflet_sorter_scripts.tcl

load ${QWRAP}/qwrap.so
load ${VEC}/vecexpr.so

proc RtoD {r} {
    global M_PI
    return [expr $r*180.0/$M_PI]
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

proc Protein_Position {name nframes hnames tnames} {
    ;# in order to use this, must have your TMD chains separated and saved as occupancy 3
    set chain_names [list "A" "B" "C" "D" "E"]

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
	set fout [open "${name}/${name}_helcoords_${eqtxt}.dat" w]
        puts $fout  "#These are the positions of your TMD helices in polar coords"
        foreach chnm $chain_names {
                set sel [atomselect top "(chain ${chnm} and name BB and occupancy 3) and (z < [expr $ht + 5] and z > [expr $ht - 5])" frame $nframes]
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
            qwrap centersel "$inpt" ;#center entire system at ~0,0,0
        }
        set com [measure center $sel weight mass]
        incr counter_i
    }
    $sel delete
    puts "Center_System finished!"
}

;# starts a new line in the print file that has the min/max r value for the bin
proc print_line_init {file number dr Rmin} {
    puts -nonewline $file "[format {%0.2f} [expr $number * $dr + $Rmin]]  [format {%0.2f} [expr ($number+1) * $dr + $Rmin]]  "
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
proc print_frame {Nr file dr Rmin Ntheta data} {
    
    #needed in order to make the proc understand this is array data
    array set data_copy $data

    ;# starts new line in outfile with radial bin values
    for {set m 0} {$m <= $Nr} {incr m} {
        print_line_init $file $m $dr $Rmin
        ;# adds bin values through penultimate value in one line
        for {set n 0} {$n < [expr $Ntheta - 1]} {incr n} {
            print_value $file $data_copy($m,$n) 0
        }
        ;# adds final value and starts new line in outfile
        set final_ndx [expr $Ntheta-1]
        print_value $file $data_copy($m,$final_ndx) 1
    }
}


;# sets all values of array to zero
proc initialize_array {Nr Ntheta init_value} {

    for {set m 0} {$m <= $Nr} {incr m} {
        for {set n 0} {$n <= $Ntheta} {incr n} {
            set data($m,$n) $init_value
        }
    }

    return [array get data]
}


;# calculate bin density
proc calculate_density {data Rmin dr dtheta Nr Ntheta totfrms} {
    
    array set data_copy $data
    for {set m 0} {$m <= $Nr} {incr m} {
        set rval [expr $m * $dr + $Rmin + [expr 0.5 * $dr]]
        set area [expr $rval * $dr * $dtheta]
        for {set n 0} {$n <= $Ntheta} {incr n} {
            set data_copy($m,$n) [expr $data_copy($m,$n) / [expr $totfrms * $area]]
        }
    }

    return [array get data_copy]
}


;# normalize bin density to get enrichment
proc normalize_density {data Nr Ntheta normfactor numbeads} {
    
    array set data_copy $data

    set normfactor [expr $normfactor * $numbeads]

    for {set m 0} {$m <= $Nr} {incr m} {
        for {set n 0} {$n <= $Ntheta} {incr n} {
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
    $sel set occupancy 3
    $sel delete
    set sel [atomselect $molid "resid 40 to 51"]
    $sel set occupancy 2
    $sel delete
    set sel [atomselect $molid "resid 52 to 65"]
    $sel set occupancy 1
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
proc cell_prep {system} {
    
    set dr 6
    set Ntheta 30
    set sample_frame 200    ;# what frame would you like to start analysis with?
    set dt 1                ;# need to fix this if you want to use it
    
    set nframes [molinfo top get numframes]

    #measure box size to get Rmax value
    set box_x [molinfo top get a frame [expr $nframes - 1]]
    set box_r [expr int($box_x) / 2]
    set Rmax [expr [expr $box_r / $dr] - [expr $dr / 2]]
    set Rmin 0
    set Rrange [expr $Rmax - $Rmin]


    set headnames "C1A C1B" ;# which beads define the surface of your membrane?
    set boxarea []

    ;#figure out which lipids are in the system
    ;#if there are other mols in your system, add them as exceptions here!
    set sel [atomselect top "not name BB SC1 to SC4 and not resname W ION"]
    set species [lsort -unique [$sel get resname]]
    $sel delete
    
    #calculate dtheta from number of theta bins
    set dthetadeg [expr 360/$Ntheta]
    global M_PI
    set dtheta [expr 2 * $M_PI / $Ntheta]
    
    #calculate number of r bins from dr and Rrange
    if {[expr $Rrange % $dr] == 0} { 
        set Nr [expr $Rrange / $dr - 1] 
    } else {
        set Nr [expr $Rrange / $dr]
    }
    
    #will need to make this less breakable
    #only works for two tails of equal length
    #only works for lipids of all same tail type
    set acyl_names [tail_analyzer $species]
    set tail_one [lindex $acyl_names 0]
    set tail_two [lindex $acyl_names 1]
    set t1T [lindex $tail_one end]
    set t2T [lindex $tail_two end]
    set tailnames "$t1T $t2T"

    #Helper scripts
    set_occupancy top ;#formats 5x29 to have separable chains and occupancies
    Center_System "resname $species"
    Align "name BB"
    leaflet_sorter $species $tailnames $nframes    ;#assigns lipids to chain U or L depending on leaflet based on 1st frame locations
    Protein_Position $system $nframes $headnames $tailnames ;#outputs a file that contains the location of the TMD helix of each monomer
}

;########################################################################################
;# polarHeight Functions

proc polarHeightByField {system} {

    set dr 6
    set Ntheta 30
    set sample_frame 200    ;# what frame would you like to start analysis with?
    set dt 1                ;# need to fix this if you want to use it
    
    set nframes [molinfo top get numframes]

    #measure box size to get Rmax value
    set box_x [molinfo top get a frame [expr $nframes - 1]]
    set box_r [expr int($box_x) / 2]
    set Rmax [expr [expr $box_r / $dr] - [expr $dr / 2]]
    set Rmin 0
    set Rrange [expr $Rmax - $Rmin]


    set headnames "C1A C1B" ;# which beads define the surface of your membrane?
    set boxarea []

    ;#figure out which lipids are in the system
    ;#if there are other mols in your system, add them as exceptions here!
    set sel [atomselect top "not name BB SC1 to SC4 and not resname W ION"]
    set species [lsort -unique [$sel get resname]]
    $sel delete
    
    #calculate dtheta from number of theta bins
    set dthetadeg [expr 360/$Ntheta]
    global M_PI
    set dtheta [expr 2 * $M_PI / $Ntheta]
    
    #calculate number of r bins from dr and Rrange
    if {[expr $Rrange % $dr] == 0} { 
        set Nr [expr [expr $Rrange / $dr] - 1.0] 
    } else {
        set Nr [expr $Rrange / $dr]
    }
    
    #will need to make this less breakable
    #only works for systems with 100% one type of lipid
    set acyl_names [tail_analyzer $species]
    set tail_one [lindex $acyl_names 0]
    set tail_two [lindex $acyl_names 1]
    set t1T [lindex $tail_one end]
    set t2T [lindex $tail_two end]
    set tailnames "$t1T $t2T"

    #Helper scripts
    set_occupancy top ;#formats 5x29 to have separable chains and occupancies
    Center_System "resname $species"
    Align "name BB"
    leaflet_sorter $species $tailnames $sample_frame    ;#assigns lipids to chain U or L depending on leaflet based on 1st frame locations
    Protein_Position $system $nframes $headnames $tailnames ;#outputs a file that contains the location of the TMD helix of each monomer

    #need to calculate heights relative to some point on the protein:
    #For 5x29 we chose the juncture between TMD and protein cap
    #because this corresponds to height zero in our elastic simulations
    set ref_bead [atomselect top "name BB and resid 15"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [vecexpr $ref_height mean]

    #outfiles setup
    set heights_up [open "${system}/${system}.zone.height.dat" w]
    set heights_down [open "${system}/${system}.ztwo.height.dat" w]
    set heights_zplus [open "${system}/${system}.zplus.height.dat" w]
    set heights_zzero [open "${system}/${system}.zzero.height.dat" w]
    set dens_up [open "${system}/${system}.zone.density.dat" w]
    set dens_down [open "${system}/${system}.ztwo.density.dat" w]
    set dens_zplus [open "${system}/${system}.zplus.density.dat" w]
    set dens_zzero [open "${system}/${system}.zzero.density.dat" w]
    set tilt_up [open "${system}/${system}.zone.tilt.dat" w]
    set tilt_down [open "${system}/${system}.ztwo.tilt.dat" w]

    puts "Helper scripts complete. Starting analysis now."	

    #position 0 is the hydrophobic interface bead; position end is the interleaflet interface bead (nominally)
    #position 0 is used for z1, z2, and zplus; position end is used for z_zero
    set heads [atomselect top "name $headnames"]
    set tails [atomselect top "(user 1 and within 6 of user 2) or (user 2 and within 6 of user 1)"]
    set t1tiltsel [atomselect top "resname $species and name GL1 $tail_one"]
    set t2tiltsel [atomselect top "resname $species and name GL2 $tail_two"]
        
    array set density_up [initialize_array $Nr $Ntheta 0.0]
    array set density_down [initialize_array $Nr $Ntheta 0.0]
    array set density_zplus [initialize_array $Nr $Ntheta 0.0]
    array set density_zzero [initialize_array $Nr $Ntheta 0.0]
    
    ;#start frame looping here
    for {set frm $sample_frame} {$frm < $nframes} {incr frm $dt} {
        set sellist [list $t1tiltsel $t2tiltsel]
        foreach selex $sellist {
            $selex frame $frm 
            $selex update
        }
        
        puts "$system $frm"

        leaflet_check $frm $species "PO4" $tailnames

        set box_height [molinfo top get c]
        
        set box_area_per_frame [expr [molinfo top get a] * [molinfo top get b]]
        lappend boxarea $box_area_per_frame

        set meas_z_zero 0

        set taillength [expr [llength $tail_one] + 1]
        set tilts [tilt_angles $taillength $t1tiltsel $t2tiltsel]

        set blist [list $heads $tails]

        foreach bead $blist {
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
            array set totals_up [initialize_array $Nr $Ntheta 0.0]
            array set counts_up [initialize_array $Nr $Ntheta 0.0]
            array set totals_down [initialize_array $Nr $Ntheta 0.0]
            array set counts_down [initialize_array $Nr $Ntheta 0.0]
            array set totals_zplus [initialize_array $Nr $Ntheta 0.0]
            array set counts_zplus [initialize_array $Nr $Ntheta 0.0]
            array set tilts_up [initialize_array $Nr $Ntheta  {0.0 0.0 0.0}]
            array set tilts_down [initialize_array $Nr $Ntheta {0.0 0.0 0.0}]

            ;#fill in total/count arrays with z sum and count sum
            for {set i 0} {$i < [llength $r_vals]} {incr i} {
                set m [lindex $r_bins $i]
                set n [lindex $theta_bins $i]
                if {$m <= $Nr} {
                    if {$meas_z_zero == 0} {
                        if {[lindex $chains $i] == 1} {
                            set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                            set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                            set density_up($m,$n) [expr {$density_up($m,$n) + 1}]
                            set tilts_up($m,$n) [vecexpr $tilts_up($m,$n) [lindex $tilts $i] add]
                        } elseif {[lindex $chains $i] == 2} {
                            set totals_down($m,$n) [expr {$totals_down($m,$n) + [lindex $z_vals $i]}]
                            set counts_down($m,$n) [expr {$counts_down($m,$n) + 1}]
                            set density_down($m,$n) [expr {$density_down($m,$n) + 1}]
                            set tilts_down($m,$n) [vecexpr $tilts_down($m,$n) [lindex $tilts $i] add]
                        }
                    } elseif {$meas_z_zero == 1} {
                        ;# reuse the up arrays for zzero
                        set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                        set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                        set density_zzero($m,$n) [expr {$density_zzero($m,$n) + 1}]
                    }
                }
            }    

            ;#turn the z sum into a z avg
            for {set m 0} {$m <= $Nr} {incr m} {
                for {set n 0} {$n <= $Ntheta} {incr n} {
                    if {$counts_up($m,$n) != 0} {
                        set totals_up($m,$n) [expr $totals_up($m,$n) / $counts_up($m,$n)]
                        if {$meas_z_zero == 0} {
                            set tilts_up($m,$n) [vecexpr $tilts_up($m,$n) $counts_up($m,$n) div]
                            set tilts_up($m,$n) [vecnorm $tilts_up($m,$n)]
                        } else {
                            set tilts_up($m,$n) "nan nan nan"
                        }
                    } else {
                        set totals_up($m,$n) "nan"
                        set tilts_up($m,$n) "nan nan nan"
                    }
                    if {$meas_z_zero == 0} {
                        if {$counts_down($m,$n) != 0} {
                            set totals_down($m,$n) [expr $totals_down($m,$n) / $counts_down($m,$n)]
                            set tilts_down($m,$n) [vecexpr $tilts_down($m,$n) $counts_down($m,$n) div]
                            set tilts_down($m,$n) [vecnorm $tilts_down($m,$n)]
                            if {$counts_up($m,$n) != 0} {
                                set totals_zplus($m,$n) [expr [expr $totals_up($m,$n) + $totals_down($m,$n)]/2.0]
                                set counts_zplus($m,$n) [expr $counts_up($m,$n) + $counts_down($m,$n)]
                            } else {
                                set totals_zplus($m,$n) "nan"
                            }
                        } else {
                            set totals_down($m,$n) "nan"
                            set totals_zplus($m,$n) "nan"
                            set tilts_down($m,$n) "nan nan nan"
                        }
                    } elseif {$meas_z_zero == 1} {
                        set totals_down($m,$n) "nan"
                        set totals_zplus($m,$n) "nan"
                        set tilts_up($m,$n) "nan nan nan"
                        set tilts_down($m,$n) "nan nan nan"
                    }
                }
            }

            ;#output heights to files
            if { $meas_z_zero == 0 } {
                print_frame $Nr $heights_up $dr $Rmin $Ntheta [array get totals_up]
                print_frame $Nr $heights_down $dr $Rmin $Ntheta [array get totals_down]
                print_frame $Nr $heights_zplus $dr $Rmin $Ntheta [array get totals_zplus]
                print_frame $Nr $tilt_up $dr $Rmin $Ntheta [array get tilts_up]
                print_frame $Nr $tilt_down $dr $Rmin $Ntheta [array get tilts_down]
            } elseif {$meas_z_zero == 1} {
                print_frame $Nr $heights_zzero $dr $Rmin $Ntheta [array get totals_up]
            }

            set meas_z_zero 1
        }
    }

    ;# calculate density
    set delta_frame [expr ($nframes - $sample_frame) / $dt]
    array set density_up [calculate_density [array get density_up] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]
    array set density_down [calculate_density [array get density_down] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]
    array set density_zplus [calculate_density [array get density_zplus] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]
    array set density_zzero [calculate_density [array get density_zzero] $Rmin $dr $dtheta $Nr $Ntheta $delta_frame]

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
    array set density_zzero [normalize_density [array get density_zzero] $Nr $Ntheta $combined_normfactor [expr [llength $tailnames] * 2.0]]  

    ;# output densities to files
    print_frame $Nr $dens_up $dr $Rmin $Ntheta [array get density_up]
    print_frame $Nr $dens_down $dr $Rmin $Ntheta [array get density_down]
    print_frame $Nr $dens_zplus $dr $Rmin $Ntheta [array get density_zplus]
    print_frame $Nr $dens_zzero $dr $Rmin $Ntheta [array get density_zzero]

    ;#clean up
    close $heights_up
    close $heights_down
    close $heights_zplus
    close $heights_zzero
    close $dens_up
    close $dens_down
    close $dens_zplus
    close $dens_zzero
    close $tilt_up 
    close $tilt_down 
    $heads delete
    $tails delete
    $t1tiltsel delete
    $t2tiltsel delete
}

proc polarHeightByBead {system} {
    
    set dr 6
    set Ntheta 30
    set sample_frame 200    ;# what frame would you like to start analysis with?
    set dt 1                ;# need to fix this if you want to use it
    
    set nframes [molinfo top get numframes]

    #measure box size to get Rmax value
    set box_x [molinfo top get a frame [expr $nframes - 1]]
    set box_r [expr int($box_x) / 2]
    set Rmax [expr [expr $box_r / $dr] - [expr $dr / 2]]
    set Rmin 0
    set Rrange [expr $Rmax - $Rmin]


    set headnames "C1A C1B" ;# which beads define the surface of your membrane?
    set boxarea []

    ;#figure out which lipids are in the system
    ;#if there are other mols in your system, add them as exceptions here!
    set sel [atomselect top "not name BB SC1 to SC4 and not resname W ION"]
    set species [lsort -unique [$sel get resname]]
    $sel delete
    
    #calculate dtheta from number of theta bins
    set dthetadeg [expr 360/$Ntheta]
    global M_PI
    set dtheta [expr 2 * $M_PI / $Ntheta]
    
    #calculate number of r bins from dr and Rrange
    if {[expr $Rrange % $dr] == 0} { 
        set Nr [expr $Rrange / $dr - 1] 
    } else {
        set Nr [expr $Rrange / $dr]
    }
    
    #will need to make this less breakable
    #only works for two tails of equal length
    #only works for systems with 100% one type of lipid
    set tail_list []
    set acyl_names [tail_analyzer $species]
    set tail_one [lindex $acyl_names 0]
    set tail_two [lindex $acyl_names 1]
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
    set t1T [lindex $tail_one end]
    set t2T [lindex $tail_two end]
    set tailnames "$t1T $t2T"

    #Helper scripts
    set_occupancy top ;#formats 5x29 to have separable chains and occupancies
    Center_System "resname $species"
    Align "occupancy 1 to 3 and name BB"
    leaflet_sorter $species $tailnames $sample_frame    ;#assigns lipids to chain U or L depending on leaflet based on 1st frame locations
    Protein_Position $system $nframes $headnames $tailnames ;#outputs a file that contains the location of the TMD helix of each monomer

    #need to calculate heights relative to some point on the protein
    #for 5x29 we chose the juncture between TMD and protein cap
    #because this corresponds to height zero in our elastic simulations
    set ref_bead [atomselect top "name BB and resid 15"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [vecexpr $ref_height mean]

    #find top of protein for pbc crossing event
    set sel [atomselect top "name BB"]
    set prot_z [$sel get z]
    $sel delete
    set protein_top [::tcl::mathfunc::max {*}$prot_z]
    set protein_top [expr $protein_top - $ref_height]

    puts "Helper scripts complete. Starting analysis now."  

    foreach beadpair $tail_list {
        set name1 [lindex $beadpair 0]
        set name2 [lindex $beadpair 1]
        set condensed_name "${name1}.${name2}"

        #outfiles setup
        set heights_up [open "${system}/${system}.${condensed_name}.zone.height.dat" w]
        set heights_down [open "${system}/${system}.${condensed_name}.ztwo.height.dat" w]
        set heights_zplus [open "${system}/${system}.${condensed_name}.zplus.height.dat" w]
        set dens_up [open "${system}/${system}.${condensed_name}.zone.density.dat" w]
        set dens_down [open "${system}/${system}.${condensed_name}.ztwo.density.dat" w]
        set dens_zplus [open "${system}/${system}.${condensed_name}.zplus.density.dat" w]

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
            print_frame $Nr $heights_up $dr $Rmin $Ntheta [array get totals_up]
            print_frame $Nr $heights_down $dr $Rmin $Ntheta [array get totals_down]
            print_frame $Nr $heights_zplus $dr $Rmin $Ntheta [array get totals_zplus]
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
        print_frame $Nr $dens_up $dr $Rmin $Ntheta [array get density_up]
        print_frame $Nr $dens_down $dr $Rmin $Ntheta [array get density_down]
        print_frame $Nr $dens_zplus $dr $Rmin $Ntheta [array get density_zplus]

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

proc run_field_mult {list_of_systems} {
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
        polarHeightByField $item
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