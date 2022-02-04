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

proc leaflet_flip_check_new {nframe beadname} {
    set lipidsel {(not resname W ION "NA+" "CL-") and (not name BB SC1 to SC4) and (name $beadname)}
    set chains {"U" "L"}
    foreach chn_nm $chains {
        set sel [atomselect top "$lipidsel and (chain $chn_nm)" frame $nframe]
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
                        $change_chain delete
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
                        $change_chain delete
                    }
                    $sel delete
                }
            }
        }
    }
}