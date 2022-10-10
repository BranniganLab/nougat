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
    foreach directory "7k3g" {
        cd $directory
        foreach item $list_of_systems {
            cd $item 
            set gro "/u1/home/js2746/Bending/PC/latest_data/${directory}/${item}/${item}.gro"
            set xtc "/u1/home/js2746/Bending/PC/latest_data/${directory}/${item}/${item}.xtc"
            #set gro "/u1/home/js2746/Bending/Jam_test/nougattest/${item}/insane.gro"
            #set xtc "/u1/home/js2746/Bending/Jam_test/nougattest/${item}/md_reduced.xtc"
            
            #set gro "/home/jesse/Bending/sims/PG/${item}.gro"
            #set xtc "/home/jesse/Bending/sims/PG/${item}.xtc"
            mol new $gro
            mol addfile $xtc waitfor all
            puts $gro
            puts $xtc
            animate delete beg 0 end 0 skip 0 top
            cd newleaf_polar
            start_nougat $item 12 30 200 -1 1 1
            cd ..
            cd newleaf_cart
            start_nougat $item 12 30 200 -1 1 0
            cd ..
            mol delete top
            cd ..
        }
        cd ..
    }
}