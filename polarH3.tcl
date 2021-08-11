
 package require pbctools
 set UTILS "~/Bending/scripts/" 
 set QWRAP "~/qwrap/"
 
source $UTILS/BinTools.tcl
source $UTILS/leaflet_sorter_scripts.tcl

load ${QWRAP}/qwrap.so



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

proc Protein_Position {} {
    set chain_names [list "A" "B" "C" "D" "E"]
    set zed [z_mid 0 20]
    puts $zed
    foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
	set fout [open "~/Bending/newanalysis/5x29_coords_${eqtxt}.dat" w]
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
    set dr 4
    set Ntheta 30

    set_occupancy top ;#formats 5x29 to have separable chains and occupancies
    Center_System "occupancy 1 to 3 and name BB"
    Align "occupancy 1 to 3 and name BB"
    Protein_Position   ;#outputs a file that contains the location of the TMD helix of each monomer
    leaflet_sorter
    set dt 1
    set nframes [molinfo top get numframes]

    set memb_up [open "${outfile}.up.height.dat" w]
    set memb_down [open "${outfile}.lw.height.dat" w]
    set sample_frame 100
    set delta_frame [expr ($nframes - $sample_frame) / $dt]
    puts "starting outer loop now"	
	for {set ri $Rmin} { $ri<=${Rmax}} { set ri [expr $ri + $dr]} {
		#loop over shells
		puts "loop $ri"
		set rf [expr $ri + $dr]
		set rf2 [expr $rf*$rf]
		set ri2 [expr $ri*$ri]
		set shell "(name PO4) and ((x*x + y*y < $rf2) and (x*x + y*y > $ri2))"
        	set shell_bin [shell_over_frames $shell $sample_frame $nframes $dt]
		set shell_up [lindex $shell_bin 0]
		set shell_down [lindex $shell_bin 1]
        	output_bins $memb_up $ri $rf $shell_up
		output_bins $memb_down $ri $rf $shell_down
		puts "loop $ri done"
	}
	close $memb_up
	close $memb_down
}
