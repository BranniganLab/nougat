proc count_frames {list_of_systems} {
    foreach directory "7k3g" {
        cd $directory
        foreach item $list_of_systems {
            cd $item 
            set gro "/u1/home/js2746/Bending/PC/latest_data/${directory}/${item}/${item}.gro"
            set xtc "/u1/home/js2746/Bending/PC/latest_data/${directory}/${item}/${item}.xtc"

            mol new $gro
            mol addfile $xtc waitfor all
            animate delete beg 0 end 0 skip 0 top
            puts "$item : [molinfo top get numframes]"
            mol delete top
            cd ..
        }
        cd ..
    }
}

proc measure_APM {list_of_systems} {
    foreach directory "7k3g" {
        cd $directory
        foreach item $list_of_systems {
            cd $item 
            set gro "/u1/home/js2746/Bending/PC/latest_data/${directory}/${item}/${item}.gro"
            set xtc "/u1/home/js2746/Bending/PC/latest_data/${directory}/${item}/${item}.xtc"
            puts $gro
            puts $xtc
            mol new $gro
            mol addfile $xtc waitfor all
            animate delete beg 0 end 0 skip 0 top
            set fout [open "${item}_boxarea.dat" w]
            if {[string match lg?? $item]} {
                set pruned_resname [string trimleft $item "lg"]
            } else {
                set pruned_resname $item
            }
            set sel [atomselect top "resname ${pruned_resname}PC"]
            set numlipids [llength [lsort -unique [$sel get resid]]]
            puts $fout "[expr $numlipids/2.0]"
            set frames [molinfo top get numframes]
            for {set frm 0} {$frm < $frames} {incr frm} {
                set area [expr [molinfo top get a frame $frm]*[molinfo top get b frame $frm]]
                puts $fout "$area"
            }
            mol delete top
            close $fout 
            cd ..
        }
        cd ..
    }
}

proc run_mult {list_of_systems} {
    foreach directory "5x29" {
        cd $directory
        foreach item $list_of_systems {
            cd $item 
            set gro "/u1/home/js2746/Bending/PC/whole_mols/${directory}/${item}/${item}.gro"
            set xtc "/u1/home/js2746/Bending/PC/whole_mols/${directory}/${item}/${item}.xtc"
            #set gro "/u1/home/js2746/Bending/Jam_test/nougattest/${item}/insane.gro"
            #set xtc "/u1/home/js2746/Bending/Jam_test/nougattest/${item}/md_reduced.xtc"
            
            #set gro "/home/jesse/Bending/sims/PG/${item}.gro"
            #set xtc "/home/jesse/Bending/sims/PG/${item}.xtc"
            mol new $gro
            mol addfile $xtc waitfor all
            puts $gro
            puts $xtc
            animate delete beg 0 end 0 skip 0 top
            start_nougat ${item} ~/PolarHeightBinning/nougat_config.txt 10 30 200 -1 1 1
            mol delete top
            cd ..
        }
        cd ..
    }
}

proc rotate_system {angle axis mol} {
    set nframes [molinfo $mol get numframes]
    set sel [atomselect $mol all]
    for {set i 0} {$i < $nframes} {incr i} {
        $sel frame $i
        $sel update
        set com [measure center $sel weight mass]
        set matrix [transaxis $axis $angle]
        $sel moveby [vecscale -1.0 $com]
        $sel move $matrix
        $sel moveby $com 
    }
}


proc transform_to_ref_height {ref} {
    set nframes [molinfo top get numframes]
    set refsel [atomselect top "$ref"]
    set totsel [atomselect top all]

    for {set i 0} {$i<$nframes} {incr i} {
        $refsel frame $i 
        $totsel frame $i 
        $refsel update
        $totsel update

        set zvals [$refsel get z]
        set avgz [expr -1.0* [vecexpr $zvals mean]]
        $totsel moveby "0 0 $avgz"
    }
    $refsel delete
    $totsel delete
}


proc read_order_params {} {
    set fp0 [open "~/Downloads/AVGPC1A.dat" r]
    set fp1 [open "~/Downloads/AVGPC1B.dat" r]
    set tail0 [read $fp0]
    set tail1 [read $fp1]
    close $fp0 
    close $fp1 

    set counter 0
    foreach dataset [list $tail0 $tail1] {
        
        set line [split $dataset "\n"]
        foreach row $line {
            set cols [split $row]
            set resd [lindex $cols 0]
            set order [lindex $cols 1]
            if {$resd eq "index"} {
                continue
            }
            if {$resd ne ""} {
                if {$counter == 0} {
                    set sel [atomselect top "resname POPC and resid $resd and name C1A D2A C3A C4A" frame 354]
                } elseif {$counter == 1} {
                    set sel [atomselect top "resname POPC and resid $resd and name C1B C2B C3B C4B" frame 354]
                } else {
                    puts "something broke"
                    break
                }
                $sel set user2 $order 

                $sel delete
            } else {
                continue
            }
        }
        set counter 1
    }
}

proc compare_lists {orders_from_nougat resids names} {
    read_order_params
    set shadow_sel [atomselect top "resname POPC and name C1A C1B" frame 354]
    set shadow_resid [$shadow_sel get resid]
    set shadow_order [$shadow_sel get user2]
    set shadow_name [$shadow_sel get name]
    $shadow_sel delete

    set orders [lindex $orders_from_nougat 0]
    set C1Alist []
    set C1Blist []

    for {set i 0} {$i < [llength $shadow_order]} {incr i} {
        set resid [lindex $shadow_resid $i]
        set order [lindex $shadow_order $i]
        set name [lindex $shadow_name $i]

        set resmatch [lsearch -all $resids $resid]
        set namematch [lsearch -all $names $name]
        foreach index $resmatch {
            if {$index in $namematch} {
                set matchindex $index 
            }
        }

        set diff [expr abs([expr $order - [lindex $orders $matchindex]])]
        if {$diff > 0.000001} {
            if {$name eq "C1B"} {
                lappend C1Blist $resid 
            } elseif {$name eq "C1A"} {
                lappend C1Alist $resid
            }
        }
    }
    puts $C1Alist
    puts $C1Blist
}

;# to use compare_lists insert the following into nougat.tcl:
                #if {$frm == 354} {
                #    compare_lists $orders [dict get $sel_info resid_list] [dict get $sel_info name_list]
                #}

proc count_lipids {sel} {
    return [llength [lsort -unique [$sel get resid]]]
}

proc track_asymmetry_over_traj {sys_name} {
    set outfile [open "${sys_name}.asymm.traj" w]
    set nframes [molinfo top get numframes]
    set sel1 [atomselect top "user 1"]
    set sel2 [atomselect top "user 2"]
    set sel3 [atomselect top "user 3"]
    set sel4 [atomselect top "user 4"]
    for {set i 0} {$i < $nframes} {incr i} {
        foreach sel [list $sel1 $sel2 $sel3 $sel4] {
            $sel frame $i 
            $sel update
        }
        set num1 [count_lipids $sel1]
        set num2 [count_lipids $sel2]
        set num3 [count_lipids $sel3]
        set num4 [count_lipids $sel4]
        puts $outfile "$i    $num1    $num2    $num3    $num4"
        puts $i
    }
    close $outfile
}

proc measure_area_APL_over_traj {sys_name lipids} {
    set outfile [open "${sys_name}.area.traj" w]
    set nframes [molinfo top get numframes]
    set sel [atomselect top "resname $lipids"]
    set num_lipids_per_leaf [expr [count_lipids $sel] / 2.0]
    for {set i 0} {$i < $nframes} {incr i} {
        set area [expr [molinfo top get a frame $i]*[molinfo top get b frame $i]]
        set apl [expr $area / $num_lipids_per_leaf]
        puts $outfile "$i    $area    $apl"
    }
    close $outfile
}





proc Select_Random_Lipids {seltext upper_final lower_final} {
    set final_list []
    for {set i 1} {$i <= 2} {incr i} {
        set sel [atomselect top "${seltext} and user ${i}"]
        set resids [lsort -unique [$sel get resid]]
        set Nl [llength $resids]
        set num_to_del [expr $Nl - $upper_final]
        if {$num_to_del < 1} {
            puts "you need to add lipids, not delete them!"
            return
        } else {
            set dellist []
            for {set j 0} {$j < $num_to_del} {incr j} {
                set resid_to_del [lindex $resids [expr {int(rand()*$Nl)}]]
                if {[lsearch -exact $dellist $resid_to_del] == -1} {
                    lappend dellist $resid_to_del
                } else {
                    set j [expr $j-1]
                }
            }
        }
        lappend final_list $dellist
    }
    return [join $final_list]
}


proc Delete_Random_Lipids {seltext upper_final lower_final outfile} {
    set list_to_delete [Select_Random_Lipids $seltext $upper_final $lower_final]
    set outputSelection [atomselect top "not ($seltext and resid ${list_to_delete})"]
    $outputSelection set beta 0
    $outputSelection writepdb $outfile 
}