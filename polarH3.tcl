
 package require pbctools
 set UTILS "utilities" 
 set QWRAP "~/qwrap-master"
 
source $UTILS/BinTools.tcl
source $UTILS/assign_tilts.tcl
source $UTILS/leaflet_sorter_scripts.tcl

load ${QWRAP}/qwrap.so
load ~/utilities/vecexpr/vecexpr.so



proc lcount {list} {
    foreach x $list {lappend arr($x) {}}
    set res1 {}
    set res2 {}
    foreach name [array names arr] {
	lappend res1 $name
	lappend res2 [llength $arr($name)]
    }
    set res [list $res1 $res2]
    return $res
 }
 
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

proc Sum_list {list_in} {
    set list_out 0
    foreach li $list_in {
        set list_out [expr 1.0*$list_out+$li]
    }
    return $list_out
}

proc z_mid {init_frm nframes} {
    set z_list {}
    for {set frm ${init_frm}} {${frm} < ${nframes}} {incr frm} {
        set mid [atomselect top "name PO4 ROH C3 PO41" frame $frm]
        lappend z_list [lindex [measure center $mid weight mass] 2]
        $mid delete
    }
    return [expr 1.0*[vecsum $z_list]/([llength $z_list]) ]
}



;# Ouputs position of the centered protein in a membrane
;# accross both leaflets
;# only useful for analysis later - doesn't impact your polar density script

proc Protein_Position {name} {
    set chain_names [list "A" "B" "C" "D" "E"]
    set zed [z_mid 0 20]
    puts $zed
    foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
	set fout [open "${name}_helcoords_${eqtxt}.dat" w]
        puts $fout  "#These are the positions of your TMD helices in polar coords"
        foreach chnm $chain_names {
                set sel [atomselect top "(chain ${chnm}) and (name BB) and (occupancy 3) and (z ${eq} $zed)" frame 0]
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

proc output_bins {fl ri rf bins} {
    puts -nonewline $fl "[format {%0.2f} $ri]  [format {%0.2f} $rf]  " 
    puts $fl "$bins" 
}

proc shell_over_frames {shell sample_frame nframes dt} {
    set singleFrame_upper []
    set singleFrame_lower []
    for {set frm $sample_frame} {$frm < $nframes} {incr frm $dt} {
        #loop over frames
        if {[expr $frm % 50] == 0} {
            leaflet_flip_check_Height_binning $frm $shell
        }
    	set singleFrame_coords [shell_frame $shell $frm]
    	if {$singleFrame_coords!=""} {
    	    if {[lindex $singleFrame_coords 0]!=""} {
                lappend singleFrame_upper [lindex $singleFrame_coords 0] 
    	    }
            if {[lindex $singleFrame_coords 1]!=""} {
                lappend singleFrame_lower [lindex $singleFrame_coords 1] 
            }
        }
    }
    set newlower []
    set newupper []
    if {$singleFrame_lower != ""} { 
    	for {set indx 0} {$indx < [llength $singleFrame_lower]} {incr indx} {
    	    for {set indx2 0} {$indx2 < [llength [lindex $singleFrame_lower $indx]]} {incr indx2} {
    	        lappend newlower [lindex [lindex $singleFrame_lower $indx] $indx2]
    	    }
    	}
    }
    if {$singleFrame_upper != ""} {
        for {set indx 0} {$indx < [llength $singleFrame_upper]} {incr indx} {
    	    for {set indx2 0} {$indx2 < [llength [lindex $singleFrame_upper $indx]]} {incr indx2} {
                    lappend newupper [lindex [lindex $singleFrame_upper $indx] $indx2]
    	    }
        }
    }
    set theta_bin_high [theta_height_avg $newupper]
    set theta_bin_low [theta_height_avg $newlower]
    return [list $theta_bin_high $theta_bin_low]
}


proc theta_height_avg {input} {   
    set theta_list []
    set height_list []
    if {$input != ""} {
        set len [llength $input]
        for {set t 0} {$t<$len} {incr t} {
	    lappend theta_list [lindex [lindex $input $t] 0]
	    lappend height_list [lindex [lindex $input $t] 1]
        }
    
    }
    set out_list []
    for {set ti 0} {$ti<30} {incr ti 1} {
        set tindex [lsearch -all $theta_list $ti] 
	if {$tindex != ""} {
            set divisor [expr 1.0*[llength $tindex]]
            set sum 0
            foreach indx $tindex {
                set newnum [lindex $height_list $indx]
                set sum [expr $sum+$newnum]
            }
            set avg [expr $sum/$divisor]
            lappend out_list $avg 
        } else {
            lappend out_list "np.nan"
        }
    }
    return $out_list
}

;# captures polar coords and z values over a shell in a single frame
proc shell_frame {shell frm} {
    set theta_high_out []
    set shell_sel [atomselect top $shell frame $frm]
    set theta_low_out []
    set nshell [$shell_sel num]
    set box_height [molinfo top get c]
    if {$nshell != 0} {
        set indexs [$shell_sel get index]
	$shell_sel delete
	foreach indx $indexs {
            set thisPO4 [atomselect top "index $indx" frame $frm]
            set coordx [$thisPO4 get x]
	    set coordy [$thisPO4 get y]
            set coordz [$thisPO4 get z]
            if {$coordz > 30} {
		   set coordz [expr $coordz - $box_height]
	    }
	    set tpchain [$thisPO4 get chain]
            $thisPO4 delete
            set theta [get_theta $coordx $coordy]
            set ti [expr int($theta/12)]
            if {$tpchain == "L"} {
                lappend theta_low_out [list $ti $coordz]
            } elseif {$tpchain == "U"} {
                lappend theta_high_out [list $ti $coordz]
            } else {
                puts "shell_frame is broken"
		exit
            }
        }
	return [list $theta_high_out $theta_low_out] 
    } else {
	$shell_sel delete    
	return
    }
}
;# THIS IS FOR 5X29!
proc set_occupancy {molid} {

  set sel [atomselect $molid "type BB or type SC1 to SC4"]
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

     if {$dist > 10} {
        lappend list1 $i
        lappend list1 [expr $i+1]
     }
     $sel1 delete
     $sel2 delete
  }

  set chars "A B C D E"
  set k 0
  lappend list1 $resmax

  foreach {i j} $list1 {
     set sel [atomselect $molid "residue $i to $j"]
     $sel set chain [lindex $chars $k]
     set k [expr $k+1]
     $sel delete
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

proc avgHlower {} {
    set sel [atomselect top "name BB and resid 1"]
    set nframes [molinfo top get numframes]


    set tot 0
    set divisor 0
    
    for {set frm 100} {$frm < $nframes} {incr frm} {
        $sel frame $frm
        $sel update
        set zs [$sel get z]
        foreach item $zs {
            set tot [expr $tot + $item]
        }
        set divisor [expr $divisor + 1]
    }
    puts $divisor
    set divisor [expr $divisor*5]
    set tot [expr $tot/$divisor]
    puts $tot
}
;########################################################################################
;# polarHeight Function

proc polarHeightByShell {outfile} {

    set Rmin 0
    set Rmax 69
    set Rrange [expr $Rmax - $Rmin]
    set dr 4
    set Ntheta 30
    set sample_frame 200
    set dt 1
    set nframes [molinfo top get numframes]
    set delta_frame [expr ($nframes - $sample_frame) / $dt]
    set counter 0
    set nm "C1A C1B"


    
    #calculate dtheta from number of theta bins
    set dtheta [expr 360/$Ntheta]
    
    #calculate number of r bins from dr and Rrange
    if {[expr $Rrange % $dr] == 0} { 
        set Nr [expr $Rrange / $dr - 1] 
    } else {
        set Nr [expr $Rrange / $dr]
    }

    #Helper scripts
    set_occupancy top ;#formats 5x29 to have separable chains and occupancies
    Center_System "occupancy 1 to 3 and name BB"
    Align "occupancy 1 to 3 and name BB"
    Protein_Position $outfile  ;#outputs a file that contains the location of the TMD helix of each monomer
    leaflet_sorter     ;#assigns lipids to chain U or L depending on leaflet based on 1st frame locations

    #find top of protein for pbc crossing event
    set sel [atomselect top "name BB"]
    set prot_z [$sel get z]
    $sel delete
    set protein_top [::tcl::mathfunc::max {*}$prot_z]

    #need to calculate heights relative to some point on the protein
    #for 5x29 we chose the juncture between TMD and protein cap
    #because this corresponds to height zero in our elastic simulations
    set ref_bead [atomselect top "name BB and resid 15"]
    set ref_height [$ref_bead get z]
    $ref_bead delete
    set ref_height [expr [vecexpr $ref_height sum]/5.0]

    #outfiles setup
    set heights_up [open "${outfile}.up.height.dat" w]
    set heights_down [open "${outfile}.lw.height.dat" w]

    puts "Helper scripts complete. Starting analysis now."	

    set lipids [atomselect top "name $nm"]
    leaflet_flip_check_new $sample_frame $nm

    #start frame looping here
    for {set frm $sample_frame} {$frm < $nframes} {set frm [expr $frm + $dt]} {
        set counter [expr $counter + 1]
        if {[expr $counter % 50] == 0} {
            leaflet_flip_check_new $frm $nm
        }

        $lipids frame $frm
        $lipids update
        set box_height [molinfo top get c]

        set x_vals [$lipids get x] 
        set y_vals [$lipids get y]
        set z_vals [vecexpr [$lipids get z] $ref_height sub]
        set chains [$lipids get chain]

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
            }
        }

        #fill in arrays with z sum and count sum
        for {set i 0} {$i < [llength $r_vals]} {incr i} {
            set m [lindex $r_bins $i]
            set n [lindex $theta_bins $i]
            if {[lindex $z_vals $i] > [expr $protein_top + 10]} {
                set [lindex $z_vals $i] [expr [lindex $z_vals $i] - $box_height]
            }
            if {$m <= $Nr} {
                if {[lindex $chains $i] == "U"} {
                    set totals_up($m,$n) [expr {$totals_up($m,$n) + [lindex $z_vals $i]}]
                    set counts_up($m,$n) [expr {$counts_up($m,$n) + 1}]
                } elseif {[lindex $chains $i] == "L"} {
                    set totals_down($m,$n) [expr {$totals_down($m,$n) + [lindex $z_vals $i]}]
                    set counts_down($m,$n) [expr {$counts_down($m,$n) + 1}]
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
                if {$counts_down($m,$n) != 0} {
                    set totals_down($m,$n) [expr $totals_down($m,$n) / $counts_down($m,$n)]
                } else {
                    set totals_down($m,$n) "nan"
                }
            }
        }
        

        #output to files
        for {set m 0} {$m <= $Nr} {incr m} {
            puts -nonewline $heights_up "[format {%0.2f} [expr $m * $dr + $Rmin]]  [format {%0.2f} [expr ($m+1) * $dr + $Rmin]]  "
            puts -nonewline $heights_down "[format {%0.2f} [expr $m * $dr + $Rmin]]  [format {%0.2f} [expr ($m+1) * $dr + $Rmin]]  "
            for {set n 0} {$n < [expr $Ntheta - 1]} {incr n} {
                puts -nonewline $heights_up " $totals_up($m,$n)" 
                puts -nonewline $heights_down " $totals_down($m,$n)" 
            }
            puts $heights_up " $totals_up($m,[expr $Ntheta-1])"
            puts $heights_down " $totals_down($m,[expr $Ntheta-1])"
        } 
    }

    close $heights_up
    close $heights_down
    $lipids delete
}

proc run_mult {list_of_systems} {
    foreach item $list_of_systems {
        set gro "/u1/home/js2746/Bending/${item}.gro"
        set xtc "/u1/home/js2746/Bending/${item}.xtc"
        mol new $gro
        mol addfile $xtc waitfor all
        puts $gro
        puts $xtc
        polarHeightByShell $item
        mol delete top
    }
}