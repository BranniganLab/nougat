proc leaflet_sorter2 {} {
    
    set lipidsel {(not resname W ION "NA+" "CL-") and (not name BB SC1 to SC4)}
    set sel [atomselect top "$lipidsel and (name PO4)" frame 1]             ;# need to make this smarter to include CHOL
    set resids [$sel get resid]
    set indexs [$sel get index]
    $sel delete
    foreach resd $resids indx $indexs {
        set lipid [atomselect top "$lipidsel and (resid $resd)" frame 1]
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



proc leaflet_flip_check_polar_binning {nframe ri rf} {
    set lipidsel {(not resname W ION "NA+" "CL-") and (not name BB SC1 to SC4) and (name PO4)}
    set chains {"U" "L"}
    set ri2 [expr $ri*$ri]
    set rf2 [expr $rf*$rf]
    foreach chn_nm $chains {
        set sel [atomselect top "$lipidsel and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2)) and (chain $chn_nm)" frame $nframe]
        set indexs [$sel get index]
        $sel delete 
        if {[llength $indexs] != 0} {
            foreach indx $indexs {
                if {$chn_nm == "U"} {
                    set sel [atomselect top "$lipidsel and (pbwithin 18 of index $indx)" frame $nframe]
                    set chns [$sel get chain]
                    set check1 [lsearch -all $chns "L"]
                    set check2 [lsearch -all $chns "U"]
                    if {[llength $check1] > [llength $check2]} {
                        set change_chain [atomselect top "index $indx"]
                        $change_chain set chain "L"
                        puts "$indx"
                    }
                    $sel delete
                } elseif {$chn_nm == "L"} {
                    set sel [atomselect top "$lipidsel and (pbwithin 18 of index $indx)" frame $nframe]
                    set chns [$sel get chain]
                    set check1 [lsearch -all $chns "U"]
                    set check2 [lsearch -all $chns "L"]
                    if {[llength $check1] > [llength $check2]} {
                        set change_chain [atomselect top "index $indx"]
                        $change_chain set chain "U"
                        puts "$indx"
                    }
                    $sel delete
                }
            }
        }
    }
}


proc leaflet_flip_check_Height_binning {nframe shell} {
    set lipidsel {(not resname W ION "NA+" "CL-") and (not name BB SC1 to SC4) and (name PO4)}
    set chains {"U" "L"}
    foreach chn_nm $chains {
        set sel [atomselect top "$shell and (chain $chn_nm)" frame $nframe]
        set indexs [$sel get index]
        $sel delete 
        if {[llength $indexs] != 0} {
            foreach indx $indexs {
                if {$chn_nm == "U"} {
                    set sel [atomselect top "$lipidsel and (pbwithin 18 of index $indx)" frame $nframe]
                    set chns [$sel get chain]
                    set check1 [lsearch -all $chns "L"]
                    set check2 [lsearch -all $chns "U"]
                    if {[llength $check1] > [llength $check2]} {
                        set change_chain [atomselect top "index $indx"]
                        $change_chain set chain "L"
                        puts "$indx"
                    }
                    $sel delete
                } elseif {$chn_nm == "L"} {
                    set sel [atomselect top "$lipidsel and (pbwithin 18 of index $indx)" frame $nframe]
                    set chns [$sel get chain]
                    set check1 [lsearch -all $chns "U"]
                    set check2 [lsearch -all $chns "L"]
                    if {[llength $check1] > [llength $check2]} {
                        set change_chain [atomselect top "index $indx"]
                        $change_chain set chain "U"
                        puts "$indx"
                    }
                    $sel delete
                }
            }
        }
    }
}

proc leaflet_flip_check_new {species nframe beadname} {
    set lipidsel "resname $species and name $beadname"
    set sel [atomselect top $lipidsel frame $nframe]
    set indexs [$sel get index] 
    $sel delete
    set inner {}
    set outer {}
    if {[llength $indexs] != 0} {
        foreach indx $indexs {
            set change_chain [atomselect top "same resid as index $indx" frame $nframe]
            set leaf [lindex [$change_chain get user] 0]
            set sel [atomselect top "$lipidsel and (pbwithin 18 of index $indx)" frame $nframe]
            if {$leaf == 1} {
                set chns [$sel get user]
                set check1 [lsearch -all $chns 0]
                set check2 [lsearch -all $chns 1]
                if {[llength $check1] > [llength $check2]} {
                    lappend inner $indx
                    puts "$indx flipped"
                } else {
                    lappend outer $indx
                }
            } elseif {$leaf == 0} {
                set chns [$sel get user]
                set check1 [lsearch -all $chns 1]
                set check2 [lsearch -all $chns 0]
                if {[llength $check1] > [llength $check2]} {
                    lappend outer $indx
                    puts "$indx flipped"
                } else {
                    lappend inner $indx
                }
            } else {
                puts "something unintended happened"
            }
            $change_chain delete
            $sel delete
        }
        set outersel [atomselect top "same resid as index $outer" frame $nframe]
        $outersel set user 1
        $outersel frame [expr $nframe + 1]
        $outersel update
        $outersel set user 1
        set innersel [atomselect top "same resid as index $inner" frame $nframe]
        $innersel set user 0
        $innersel frame [expr $nframe + 1]
        $innersel update
        $innersel set user 0
        $outersel delete
        $innersel delete
    } else {
        puts "something broke"
    }
}

#written in 2022 to be fast and dirty
#should rewrite at some point to use only user values and not chains
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
    #rewrite to make protein position passable value
    set sel [atomselect top "name BB and resid 30" frame $frm]
    set com [measure center $sel]
    set x [lindex $com 0]
    set y [lindex $com 1]
    $sel delete
    set porelipids [atomselect top "(resname $species and same resid as within 9 of resid 30) or (resname $species and same resid as ((x-$x)*(x-$x)+(y-$y)*(y-$y) <= 16))" frame $frm]
    $porelipids set chain "Z"
    $porelipids delete

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

#experimental - didn't work
proc rmv_outliers {frm species headnames lipidlength} {
    set change_list []

    set sel [atomselect top "name $headnames and user 1" frame $frm]
    set ids [$sel get resid]
    set z_vals [$sel get z]

    set avgz [vecexpr $z_vals mean]

    set cutoff [expr $avgz + $lipidlength]
    for {set ndx 0} {$ndx <= [llength $z_vals]} {incr ndx} {
        if {[lindex $z_vals $ndx] > $cutoff} {
            lappend change_list [lindex $ids $ndx]
        }
    }

    $sel delete        
    if {[llength $change_list] > 0} {
        set sel [atomselect top "resname $species and resid $change_list" frame $frm]
        $sel set user 4
        $sel delete
    }
}