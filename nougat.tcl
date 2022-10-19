package require pbctools

# EDIT THE PATH HERE
# TELL nougat WHERE TO FIND YOUR UTILITIES FOLDER
set UTILS "~/PolarHeightBinning/utilities"
source ${UTILS}/helper_procs.tcl

set CONFIG_PATH "~/PolarHeightBinning/nougat_config.txt"

proc cell_prep {config_path leaf_check} {

    set config_dict [read_config_file $config_path]

    load [dict get $config_dict qwrap_path]/qwrap.so
    load [dict get $config_dict vecexpr_path]/vecexpr.so

    ;# figures out which lipids and beads are in the system
    ;# no edits required
    set lipidsel [atomselect top "not [dict get $config_dict inclusion_sel] and not [dict get $config_dict excluded_sel]"]
    set species [lsort -unique [$lipidsel get resname]]
    set tail_info [tail_analyzer $species]
    set acyl_names [lindex $tail_info 0]
    set heads_and_tails [lindex $tail_info 1]
    set full_tails [lindex $tail_info 2]
    $lipidsel delete

    ;# need to manually substitute variables into the config_dict when applicable
    foreach key [dict keys $config_dict] {
        dict set config_dict $key [subst [dict get $config_dict $key]]
    }

    ;#****************************************************;#
    ;#          MAKE EDITS BELOW BEFORE STARTING          ;#
    ;#****************************************************;# 

    ;# center, wrap, and align the system
    ;# if your inclusion 'tumbles' in the membrane (like a nanoparticle) comment out Align!
    Center_System [dict get $config_dict wrap_sel]
    Align [dict get $config_dict align_sel]

    ;# custom proc to set my TMD helices to occupancy 1
    ;# this allows Protein_Position to work
    ;# comment this out or customize it for your inclusion
    #set_occupancy top 

    ;# this will only work if your TMD helices are set to occupancy 1
    ;# otherwise, comment it out
    ;# all it does is put a dot on the polar heatmap where a TMD helix should be, so not essential at all
    ;# this proc is currently broken, anyways (yes, there is a github issue)
    ;#Protein_Position $system $headnames $acyl_names ;# FIX ME

    ;# set user3 to hold a unique tail number for easy separation of tails later
    tail_numberer $species $acyl_names

    ;# sets user to 1 or 2 depending on if the lipid is in the outer or inner leaflet
    ;# sets user to 3 if the lipid is too horizontal to determine leaflet
    ;# sets user to 4 if you have pore_sorter turned on
    ;# TO USE PORE_SORTER YOU NEED TO CUSTOMIZE IT TO YOUR PROTEIN!!!
    if {$leaf_check == 1} {
        set end [molinfo top get numframes]
        for {set i 0} {$i < $end} {incr i} {
            puts $i
            leaflet_check $i $species $heads_and_tails 1.0
        }
    }

    ;#****************************************************;#
    ;#          MAKE EDITS ABOVE BEFORE STARTING          ;#
    ;#****************************************************;# 

    dict set config_dict species $species 
    dict set config_dict acyl_names $acyl_names
    dict set config_dict full_tails $full_tails
    dict set config_dict heads_and_tails $heads_and_tails
    
    return $config_dict  
}



;########################################################################################
;# polarHeight Functions

proc start_nougat {system CONFIG_PATH d1 N2 start end step polar} {

    ;# running cell_prep will do some important initial configuration based on user input. 
    ;# check the extensive documentation at the top of this file for instructions.
    set config_dict [cell_prep $CONFIG_PATH 0]

    ;# set nframes based on $end input
    set maxframes [molinfo top get numframes]
    if {$end == -1} {
        set nframes [expr $maxframes-1]
    } elseif {($end < $maxframes) || ($end > $start)} {
        set nframes $end
    } else {
        puts "you specified a frame number that is outside the allowable range"
        break
    }

    ;# change this value if you want to exclude an inner radius in polar coords
    set min 0 

    ;# determine number and size of bins
    set bindims [bin_prep $nframes $polar $min $d1 $N2]

    if {$polar == 0} {
        set bindims [update_dims $bindims 0]
    }

    ;# add all these new values to important_variables for easy transfer
    dict set config_dict start $start 
    dict set config_dict nframes $nframes 
    dict set config_dict step $step 
    dict set config_dict min $min

    ;# run nougat twice, once to compute height and density and once to compute
    ;# lipid tail vectors and order parameters
    run_nougat $system $config_dict $bindims $polar "height_density" 
    run_nougat $system $config_dict $bindims $polar "tilt_order" 
}

proc run_nougat {system config_dict bindims polar quantity_of_interest} {  
    
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
    set outfiles [create_outfiles $system $quantity_of_interest [concat_names [dict get $config_dict headnames]] [dict get $config_dict species] [dict get $config_dict acyl_names] $coordsys]
     
    ;#atomselections setup as dict
    if {$quantity_of_interest eq "height_density"} {
        dict set selections z1z2 [atomselect top "resname [dict get $config_dict species] and name [dict get $config_dict full_tails]"]
        dict set selections z0 [atomselect top "resname [dict get $config_dict species] and ((user 1 and within 6 of user 2) or (user 2 and within 6 of user 1))"]
    } elseif {$quantity_of_interest eq "tilt_order"} {
        set lists [tail_length_sorter [dict get $config_dict species] [dict get $config_dict acyl_names]]
        set sellist [lindex $lists 0]
        set lenlist [lindex $lists 1]
        foreach sel $sellist len $lenlist {
            dict set selections $len [atomselect top "$sel"] 
        }
    }

    puts "Setup complete. Starting frame analysis now."   

    ;# start frame looping here
    for {set frm [dict get $config_dict start]} {$frm <= [dict get $config_dict nframes]} {incr frm [dict get $config_dict step]} {

        if {$polar == 0} {
            set bindims [update_dims $bindims $frm]
        }
        
        ;# update leaflets in case lipids have flip-flopped
        leaflet_check $frm [dict get $config_dict species] [dict get $config_dict heads_and_tails] 1.0

        puts "$system $quantity_of_interest $frm"

        #need to calculate heights relative to some point (usually on the inclusion):
        set ref_height [calc_ref_height $config_dict $frm]
        
        ;# height_density has two selections, so this will execute twice.
        ;# tilt_order has different selections, one for each tail length present
        ;# in the system, so this will execute as many times as there are
        ;# unique tail lengths.
        foreach selex [dict keys $selections] {
            
            ;# $selex is a dict key that holds an atomselection as its value
            set sel [dict get $selections $selex]

            $sel frame $frm 
            $sel update

            ;# assemble all data (x,y,z,user, etc) into a dict of lists
            set sel_info [grab_sel_info $sel $ref_height]

            ;# calculate which bins each bead belongs in along both axes
            ;# and return as two lists of same length as the lists above
            set bins [bin_assigner [dict get $sel_info xvals_list] [dict get $sel_info yvals_list] [dict get $bindims d1] [dict get $bindims d2] [dict get $bindims dthetadeg] $polar]
            set dim1_bins_list [lindex $bins 0]
            set dim2_bins_list [lindex $bins 1]

            ;# Binning is controlled by the bead designated in $headnames.
            ;# Creates a dict that contains the bin and leaflet information linked to
            ;# a resid and index number. Facilitates easy binning later. 
            set res_dict [create_res_dict [dict get $config_dict species] [dict get $config_dict headnames] [dict get $sel_info lipid_list] [dict get $sel_info name_list] [dict get $sel_info resid_list] $dim1_bins_list $dim2_bins_list [dict get $sel_info leaflet_list] $selex]

            ;# Make necessary calculations (if any), then bin and average them
            if {$quantity_of_interest eq "height_density"} {
                set outfiles [height_density_averaging $res_dict $outfiles [dict get $sel_info leaflet_list] [dict get $sel_info lipid_list] [dict get $sel_info zvals_list] [dict get $sel_info name_list]]
            } elseif {$quantity_of_interest eq "tilt_order"} {
                set tilts [tilt_angles [dict keys $selections] [dict get $sel_info xvals_list] [dict get $sel_info yvals_list] [dict get $sel_info zvals_list]]
                set orders [order_params [dict keys $selections] [dict get $sel_info xvals_list] [dict get $sel_info yvals_list] [dict get $sel_info zvals_list]]
                set outfiles [tilt_order_averaging $res_dict $outfiles [dict get $sel_info leaflet_list] [dict get $sel_info lipid_list] $tilts $orders [dict get $sel_info tail_list] $selex]
            }

            ;# Now that all information has been binned, print it to files
            dict for {key val} [dict get $outfiles $selex] {
                print_frame [dict get $bindims N1] $outfiles $key [dict get $bindims d1] [dict get $config_dict min] [dict get $bindims N2] $polar $selex

                ;# cleanup before next step
                set outfiles [dict unset outfiles $selex $key bin]
            } 

            ;# cleanup before next step
            dict remove res_dict
            dict remove sel_info
        }
    }

    ;# output density normalization info 
    if {$quantity_of_interest eq "height_density"} {
        output_density_norm_info [dict get $config_dict start] [dict get $config_dict nframes] [dict get $config_dict step] [dict get $config_dict species] $system [dict get $config_dict headnames] $coordsys
    }

    ;# close all outfiles
    foreach channel [file channels "file*"] {
        close $channel
    }

    ;# delete all atomselections in scope
    ;# if the user defined atomselections in their VMD session,
    ;# these will be out of scope and cause an error.
    ;# catch will ignore this error, as it's not important to the user.
    foreach selection [atomselect list] {
        catch {$selection delete}
    }
}
