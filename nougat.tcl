package require pbctools

# EDIT THE PATHS HERE
# TELL nougat WHERE TO FIND YOUR VERSIONS OF qwrap AND vecexpr
#set QWRAP "~/qwrap-master"
#set VEC "~/utilities/vecexpr"

set UTILS "~/PolarHeightBinning/utilities"

source ${UTILS}/helper_procs.tcl
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
    set excluded_sel "resname W CHOL ION 'CL-' 'NA+' lig"

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
    #set align_sel "name BB"

    ;# provide atomselection-style text that defines the reference point that should correspond with height 0 (on average) in your plots.
    ;# E.G. for 5x29 we decided resid 15 would be the 'zero-point' and all heights would be provided with reference to 
    ;# the average position of resid 15
    ;# IF YOU DO NOT WISH TO SET A REFERENCE POINT:
    ;# replace the text with "NULL"
    set reference_point "NULL"

    ;# provide the beadnames that you consider to form the surface of your membrane
    ;# we chose the top tail beads because they are what form the 'hydrophobic surface'
    ;# in our opinion
    set headnames "C1A C1B"

    ;# center, wrap, and align the system
    ;# if your inclusion 'tumbles' in the membrane (like a nanoparticle) comment out Align!
    Center_System "$wrap_sel"
    #Align "$align_sel"

    ;# custom proc to set my TMD helices to occupancy 1
    ;# this allows Protein_Position to work
    ;# comment this out or customize it for your inclusion
    #set_occupancy top 

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

    #need to calculate heights relative to some point (usually on the inclusion):
    if {$reference_point ne "NULL"} {
        set ref_bead [atomselect top "$reference_point"]
        set ref_height [$ref_bead get z]
        $ref_bead delete
        set ref_height [vecexpr $ref_height mean]
    } else {
        set ref_height "NULL"
    }

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
    set bindims [bin_prep $nframes $polar $min $d1 $N2]

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
        dict set selections z0 [atomselect top "resname $species and ((user 1 and within 6 of user 2) or (user 2 and within 6 of user 1))"]
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
        
        ;# height_density has two selections, so this will execute twice.
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
            if {$ref_height ne "NULL"} {
                set zvals_list [vecexpr [$sel get z] $ref_height sub]
            } else {
                set zvals_list [$sel get z]
            }   
            
            ;# user contains a 1 or 2 for outer or inner leaflet, respectively
            set leaflet_list [$sel get user]

            ;# user3 contains an int that describes which tail in the lipid this is
            ;# E.G. POPC will have 0 or 1 (it has two tails)
            ;# E.G. OANT will have 0, 1, 2, 3, 4, or 5 (it has 6 tails)
            set tail_list [$sel get user3]
            set tail_list [vecexpr $tail_list 1 sub]

            ;# calculate which bins each bead belongs in along both axes
            ;# and return as two lists of same length as the lists above
            set bins [bin_assigner $xvals_list $yvals_list $d1 $d2 $dthetadeg $polar]
            set dim1_bins_list [lindex $bins 0]
            set dim2_bins_list [lindex $bins 1]

            ;# Binning is controlled by the bead designated in $headnames.
            ;# Creates a dict that contains the bin and leaflet information linked to
            ;# a resid and index number. Facilitates easy binning later. 
            set res_dict [create_res_dict $species $headnames $lipid_list $name_list $resid_list $dim1_bins_list $dim2_bins_list $leaflet_list $selex]

            ;# Make necessary calculations (if any) and then bin them
            if {$quantity_of_interest eq "height_density"} {
                set outfiles [do_height_density_binning $res_dict $outfiles $leaflet_list $lipid_list $zvals_list]
            } elseif {$quantity_of_interest eq "tilt_order"} {
                set tilts [tilt_angles [dict keys $selections] $xvals_list $yvals_list $zvals_list]
                set orders [order_params [dict keys $selections] $xvals_list $yvals_list $zvals_list]
                set outfiles [do_tilt_order_binning $res_dict $outfiles $leaflet_list $lipid_list $tilts $orders $tail_list $selex]
            }

            ;# Now that all information has been binned, print it to files
            dict for {key val} [dict get $outfiles $selex] {
                print_frame $N1 $outfiles $key $d1 $min $N2 $polar $selex

                ;# cleanup before next step
                set outfiles [dict unset outfiles $selex $key bin]
            } 

            ;# cleanup before next step
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

    ;# output density normalization info 
    if {$quantity_of_interest eq "height_density"} {
        output_density_norm_info $start $nframes $step $species $system $headnames $coordsys
    }
}

;# Need to rewrite so that it works with all the new settings
proc run_mult {list_of_systems} {
    foreach directory "5x29" {
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
            start_nougat $item 12 30 200 -1 1 1
            start_nougat $item 12 30 200 -1 1 0
            mol delete top
            cd ..
        }
        cd ..
    }
}