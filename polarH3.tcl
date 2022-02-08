
package require pbctools

#locations on JS home computer
#set UTILS "~/Bending/scripts/PolarHeightBinning/utilities" 
#set QWRAP "~/qwrap"
#set VEC "~/vecexpr"

#locations on Belenus
set UTILS "~/PolarHeightBinning/utilities"
set QWRAP "~/qwrap-master"
set VEC "~/utilities/vecexpr"

source $UTILS/BinTools.tcl
source $UTILS/assign_tilts.tcl
source $UTILS/leaflet_sorter_scripts.tcl

load ${QWRAP}/qwrap.so
load ${VEC}/vecexpr.so

 proc RtoD {r} {
    set pi 3.14159265358979323846
    return [expr $r*180.0/$pi]
}

proc get_theta {x y} {
    set pi 3.14159265358979323846
    set tmp  [expr {atan2($y,$x)}]
    if {$tmp < 0} {
        set theta [expr 2*$pi + $tmp]    
    } else {
        set theta $tmp
    }
    return [RtoD $theta]
}

;# Ouputs position of the centered protein in a membrane
;# accross both leaflets
;# only useful for analysis later - doesn't impact your polar density script

proc Protein_Position {name nframes hnames tnames} {
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
	set fout [open "${name}_helcoords_${eqtxt}.dat" w]
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

            if {$dist > 20} {
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


proc leaflet_sorter {} {

    set sel [atomselect top "name PO4" frame 1]
    set resids [$sel get resid]
    set indexs [$sel get index]
    $sel delete
    foreach resd $resids indx $indexs {
        set lipid [atomselect top "(not name BB SC1 to SC4 W and not resname ION) and (resid $resd)" frame 1]
        set po4 [atomselect top "index $indx" frame 1]
        set lipid_com [measure center $lipid weight mass]
        set lipid_com_z [lindex $lipid_com 2]
        set po4_z [$po4 get z]
        set dist [expr abs($po4_z-$lipid_com_z)] 
        if {$dist > 40} {
            $lipid set chain L
        } else {
            if {$po4_z > $lipid_com_z} {
                $lipid set chain U
            } else {
                $lipid set chain L
            }
        }
        $lipid delete
        $po4 delete
    }
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

    set sel [atomselect top "resname $species and name $tail_one"] 
    $sel set user 1
    $sel delete

    set sel [atomselect top "resname $species and name $tail_two"]
    $sel set user 2
    $sel delete 
    return [list $tail_one $tail_two]
}

;########################################################################################
;# polarHeight Function

proc polarHeightByShell {outfile} {


    set Rmin 0
    set Rmax 69
    set Rrange [expr $Rmax - $Rmin]
    set dr 2
    set Ntheta 30
    set sample_frame 200
    set dt 1
    set nframes [molinfo top get numframes]
    set delta_frame [expr ($nframes - $sample_frame) / $dt]
    set counter 0
    set num_subunits 5.0
    set headgrps [list "PC" "PG"]
    
    #create list of lipids used
    #will need to make this less breakable
    set tailnms $outfile
    set species ""
    foreach head $headgrps {
        set species "$species$tailnms$head "
    }
    
    #calculate dtheta from number of theta bins
    set dtheta [expr 360/$Ntheta]
    
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

    set t1H [lindex $tail_one 0]
    set t2H [lindex $tail_two 0]
    set t1T [lindex $tail_one end]
    set t2T [lindex $tail_two end]

    set headnames "$t1H $t2H"
    set tailnames "$t1T $t2T"

    #Helper scripts
    set_occupancy top ;#formats 5x29 to have separable chains and occupancies
    Center_System "occupancy 1 to 3 and name BB"
    Align "occupancy 1 to 3 and name BB"
    leaflet_sorter     ;#assigns lipids to chain U or L depending on leaflet based on 1st frame locations
    Protein_Position $outfile $nframes $headnames $tailnames ;#outputs a file that contains the location of the TMD helix of each monomer

    #need to calculate heights relative to some point on the protein
    #for 5x29 we chose the juncture between TMD and protein cap
    #because this corresponds to height zero in our elastic simulations
    set ref_bead [atomselect top "name BB and resid 15"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [expr [vecexpr $ref_height sum]/$num_subunits]

    #find top of protein for pbc crossing event
    set sel [atomselect top "name BB"]
    set prot_z [$sel get z]
    $sel delete
    set protein_top [::tcl::mathfunc::max {*}$prot_z]
    set protein_top [expr $protein_top - $ref_height]

    #outfiles setup
    set heights_up [open "${outfile}.zone.height.dat" w]
    set heights_down [open "${outfile}.ztwo.height.dat" w]
    set heights_zplus [open "${outfile}.zplus.height.dat" w]
    set heights_zzero [open "${outfile}.zzero.height.dat" w]

    puts "Helper scripts complete. Starting analysis now."	

    #position 0 is the hydrophobic interface bead; position end is the interleaflet interface bead (nominally)
    #position 0 is used for z1, z2, and zplus; position end is used for z_zero
    set heads [atomselect top "name $headnames"]
    set tails [atomselect top "((name $tailnames and chain U) and within 6 of (name $tailnames and chain L)) or ((name $tailnames and chain L) and within 6 of (name $tailnames and chain U))"]
        
    leaflet_flip_check_new $sample_frame "C1A C1B"

    #start frame looping here
    for {set frm $sample_frame} {$frm <= $nframes} {set frm [expr $frm + $dt]} {
        puts $frm
        set counter [expr $counter + 1]
        if {[expr $counter % 50] == 0} {
            leaflet_flip_check_new $frm "C1A C1B"
        }

        set box_height [molinfo top get c]
        set bead_counter 0

        set blist [list $heads $tails]

        foreach bead $blist {
            $bead frame $frm 
            $bead update

            set x_vals [$bead get x] 
            set y_vals [$bead get y]
            set z_vals [vecexpr [$bead get z] $ref_height sub]
            set chains [$bead get chain]
            set resids [$bead get resid]
            set indexs [$bead get index]

            #get theta values for all x,y pairs
            set theta_vals [vecexpr $y_vals $x_vals atan2 pi div 180 mult]  

            #atan2 gives values from -180 to 180; shifting to 0 to 360                            
            for {set i 0} {$i<[llength $theta_vals]} {incr i} {                                         
                if {[lindex $theta_vals $i] < 0} {                                                      
                    set theta_vals [lreplace $theta_vals $i $i [expr [lindex $theta_vals $i]+360]]
                }
            }

            #turn into bin numbers rather than theta values
            vecexpr [vecexpr $theta_vals $dtheta div] floor &theta_bins
            
            #calculate distance from origin for all x,y pairs
            set r_vals [vecexpr [vecexpr [vecexpr $x_vals sq] [vecexpr $y_vals sq] add] sqrt]
            
            #turn into bin numbers rather than r values
            vecexpr [vecexpr $r_vals $dr div] floor &r_bins

            #initialize arrays to zeros
            for {set m 0} {$m <= $Nr} {incr m} {
                for {set n 0} {$n <= $Ntheta} {incr n} {
                    set totals_up($m,$n) 0
                    set counts_up($m,$n) 0
                    set totals_down($m,$n) 0
                    set counts_down($m,$n) 0
                    set midpoint($m,$n) 0
                }
            }

            #fill in arrays with z sum and count sum
            for {set i 0} {$i < [llength $r_vals]} {incr i} {
                set m [lindex $r_bins $i]
                set n [lindex $theta_bins $i]
                if {[lindex $z_vals $i] > [expr $protein_top + 5]} {
                    set [lindex $z_vals $i] [expr [lindex $z_vals $i] - $box_height]
                }
                if {$m <= $Nr} {
                    if {$bead_counter == 0} {
                        if {[lindex $chains $i] == "U"} {
                            set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                            set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                        } elseif {[lindex $chains $i] == "L"} {
                            set totals_down($m,$n) [expr {$totals_down($m,$n) + [lindex $z_vals $i]}]
                            set counts_down($m,$n) [expr {$counts_down($m,$n) + 1}]
                        }
                    } elseif {$bead_counter == 1} {
                        set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                        set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                    }
                }
            }    

            #turn the z sum into a z avg
            for {set m 0} {$m <= $Nr} {incr m} {
                for {set n 0} {$n <= $Ntheta} {incr n} {
                    if {$counts_up($m,$n) != 0} {
                        set totals_up($m,$n) [expr $totals_up($m,$n) / $counts_up($m,$n)]
                    } else {
                        set totals_up($m,$n) "nan"
                    }
                    if {$bead_counter == 0} {
                        if {$counts_down($m,$n) != 0} {
                            set totals_down($m,$n) [expr $totals_down($m,$n) / $counts_down($m,$n)]
                            if {$counts_up($m,$n) != 0} {
                                set midpoint($m,$n) [expr [expr $totals_up($m,$n) + $totals_down($m,$n)]/2.0]
                            } else {
                                set midpoint($m,$n) "nan"
                            }
                        } else {
                            set totals_down($m,$n) "nan"
                            set midpoint($m,$n) "nan"
                        }
                    } elseif {$bead_counter == 1} {
                        set totals_down($m,$n) "nan"
                        set midpoint($m,$n) "nan"
                    }
                }
            }

            #output to files
            if { $bead_counter == 0 } {
                for {set m 0} {$m <= $Nr} {incr m} {
                    puts -nonewline $heights_up "[format {%0.2f} [expr $m * $dr + $Rmin]]  [format {%0.2f} [expr ($m+1) * $dr + $Rmin]]  "
                    puts -nonewline $heights_down "[format {%0.2f} [expr $m * $dr + $Rmin]]  [format {%0.2f} [expr ($m+1) * $dr + $Rmin]]  "
                    puts -nonewline $heights_zplus "[format {%0.2f} [expr $m * $dr + $Rmin]]  [format {%0.2f} [expr ($m+1) * $dr + $Rmin]]  "
                    for {set n 0} {$n < [expr $Ntheta - 1]} {incr n} {
                        puts -nonewline $heights_up " $totals_up($m,$n)" 
                        puts -nonewline $heights_down " $totals_down($m,$n)"
                        puts -nonewline $heights_zplus " $midpoint($m,$n)"
                    }
                    puts $heights_up " $totals_up($m,[expr $Ntheta-1])"
                    puts $heights_down " $totals_down($m,[expr $Ntheta-1])"
                    puts $heights_zplus " $midpoint($m,[expr $Ntheta-1])"
                } 
            } elseif {$bead_counter == 1} {
                for {set m 0} {$m <= $Nr} {incr m} {
                    puts -nonewline $heights_zzero "[format {%0.2f} [expr $m * $dr + $Rmin]]  [format {%0.2f} [expr ($m+1) * $dr + $Rmin]]  "
                    for {set n 0} {$n < [expr $Ntheta - 1]} {incr n} {
                        puts -nonewline $heights_zzero " $totals_up($m,$n)"
                    }
                    puts $heights_zzero " $totals_up($m,[expr $Ntheta-1])"
                }
            }
            set bead_counter 1
        }
    }

    close $heights_up
    close $heights_down
    close $heights_zplus
    close $heights_zzero
    $heads delete
    $tails delete
}

proc run_mult {list_of_systems} {
    foreach item $list_of_systems {
        set gro "/u1/home/js2746/Bending/PC/${item}.gro"
        set xtc "/u1/home/js2746/Bending/PC/${item}.xtc"
        mol new $gro
        mol addfile $xtc waitfor all
        puts $gro
        puts $xtc
        polarHeightByShell $item
        #set nframes [molinfo top get numframes]
        #puts "$item has $nframes frames"
        mol delete top
    }
}
