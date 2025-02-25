 package require pbctools

;# read in the config file and save settings to dictionary key/value pairs
proc read_config_file {path} {
    
    ;# read file and close it
    set fp [open $path r]
    set file_data [read $fp]
    close $fp

    set data [split $file_data "\n"]
    foreach line $data {
        ;# ignore comment lines
        if {[string match "#*" $line]} {
            continue
        }
        ;# ignore empty lines
        if {$line eq ""} {
            continue
        }
        ;# save what's left to config_dict
        set key_val [split $line "="]
        dict set config_dict [string trim [lindex $key_val 0]] [string trim [lindex $key_val 1]]
    }

    return $config_dict
}

proc cell_prep {config_path leaf_check} {

    set config_dict [read_config_file $config_path]

    source [dict get $config_dict utilities_path]/helper_procs.tcl
    
    if {[dict get $config_dict use_qwrap] == "yes"} {
        load [dict get $config_dict qwrap_path]/qwrap.so
    }
    
    if {[dict get $config_dict use_vecexpr] == "yes"} {
        load [dict get $config_dict vecexpr_path]/vecexpr.so
    }

    ;# figure out which lipids are in the system

    if {[dict get $config_dict inclusion_sel] ne "NULL"} {
        set lipidsel [atomselect top "not [dict get $config_dict inclusion_sel] and not [dict get $config_dict excluded_sel]"]
    } else {
        set lipidsel [atomselect top "not [dict get $config_dict excluded_sel]"]
    }

    set species [lsort -unique [$lipidsel get resname]]
    $lipidsel delete

    set tail_info [analyzeTails $species]

    set acyl_names [lindex $tail_info 0]
    set heads_and_tails [lindex $tail_info 1]
    set full_tails [lindex $tail_info 2]
    
    ;# need to manually substitute variables into the config_dict when applicable
    foreach key [dict keys $config_dict] {
        dict set config_dict $key [subst [dict get $config_dict $key]]
    }

    set end [molinfo top get numframes]

    ;# wrap and center the system
    if {[dict get $config_dict use_qwrap] == "yes"} {
        run_qwrap [dict get $config_dict wrap_sel] [dict get $config_dict inclusion_sel]
    } elseif {[dict get $config_dict inclusion_sel] != "NULL"} {
        pbc wrap -all -center bb -centersel [dict get $config_dict inclusion_sel] 
        center_system [dict get $config_dict inclusion_sel] [dict get $config_dict wrap_sel] $end
    } else {
        center_system [dict get $config_dict inclusion_sel] [dict get $config_dict wrap_sel] $end
    }

    # align the system, removing xy rotation
    if {[dict exists $config_dict align_sel] && [dict get $config_dict align_sel] != "NULL"} {
        Align [dict get $config_dict align_sel] 0
    }

    ;# custom proc to set my TMD helices to occupancy 1
    ;# this allows Protein_Position to work
    if {[dict exists $config_dict custom_occupancy] && ([dict get $config_dict custom_occupancy] != "NULL")} {
        separate_chains top 10
        set_occupancy top
    }

    ;# sets user to 1 or 2 depending on if the lipid is in the outer or inner leaflet
    ;# sets user to 3 if the lipid is too horizontal to determine leaflet
    ;# sets user to 4 if you have pore_sorter turned on
    if {$leaf_check == 1} {
        for {set i 0} {$i < $end} {incr i} {
            puts "leaflet check frame $i"
            assignLeaflet $i $species $heads_and_tails 1.0 [dict get $config_dict pore_sorter]
        }
    }

    ;# set user3 to hold a unique tail number for easy separation of tails later
    numberTails $species $acyl_names


    dict set config_dict species $species 
    dict set config_dict acyl_names $acyl_names
    dict set config_dict full_tails $full_tails
    dict set config_dict heads_and_tails $heads_and_tails
    
    puts "cell_prep complete"

    return $config_dict  
}



;########################################################################################
;# nougat main functions

proc start_nougat {system config_path dr_N1 N2 start end step polar {switch 0}} {
    ;# running cell_prep will do some important initial configuration based on user input. 
    set config_dict [cell_prep $config_path 0]

    set coordsys [readPolar $polar]
    set foldername "${system}_${coordsys}_${dr_N1}_${N2}_${start}_${end}_${step}"
    file mkdir $foldername
    file copy -force $config_path $foldername

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
    set bindims [prepareBins $nframes $polar $min $dr_N1 $N2]

    ;# add all these new values to important_variables for easy transfer
    dict set config_dict start $start 
    dict set config_dict nframes $nframes 
    dict set config_dict step $step 
    dict set config_dict min $min

    ;# run nougat twice, once to compute height and density and once to compute
    ;# lipid tail vectors and order parameters
    
    run_nougat $system $config_dict $bindims $polar "height" $foldername $switch
    #run_nougat $system $config_dict $bindims $polar "tilt_order" $foldername

}

proc run_nougat {system config_dict bindims polar quantity_of_interest foldername {switch 0}} {  
    
    set coordsys [readPolar $polar]
    set outfiles [createOutfiles $quantity_of_interest [dict get $config_dict species] [dict get $config_dict acyl_names] $foldername]
    set selections [createAtomSelections $quantity_of_interest $config_dict]
    if {$switch == 0} {
        puts "#### Writing Coordinates To File ####"
        set files [open "full_file.dat" w]
        writeHeader $files 
    } 

    puts "Setup complete. Starting frame analysis now."   

    ;# start frame looping here
    for {set frm [dict get $config_dict start]} {$frm <= [dict get $config_dict nframes]} {incr frm [dict get $config_dict step]} {
        if {$polar == 0} {
            set bindims [updateDimensions $bindims $frm]
            
        }

        ;# update leaflets in case lipids have flip-flopped

        assignLeaflet $frm [dict get $config_dict species] [dict get $config_dict heads_and_tails] 1.0 [dict get $config_dict pore_sorter]

        puts "$system $quantity_of_interest $frm"

        #need to calculate heights relative to some point (usually on the inclusion):
        set ref_height [calculateReferenceHeight $config_dict $frm]
       
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
            set sel_info [getSelInfo $sel $ref_height]
            
            if {$switch == 0 && $selex == "z1z2"} {
                puts $files [$sel get x] 
                puts $files [$sel get y]
                puts $files [$sel get z] 
                puts $files [$sel get resid]
                puts $files [$sel get index]
                puts $files [$sel get user]

            } else { 
            
            #### This is where binning happens, straddle here ####

            ;# calculate which bins each bead belongs in along both axes
            ;# and return as two lists of same length as the lists above
            set bins [assignBins [dict get $sel_info xvals_list] [dict get $sel_info yvals_list] [dict get $bindims d1] [dict get $bindims d2] [dict get $bindims dthetadeg] $polar $frm [dict get $config_dict use_vecexpr]]
            set dim1_bins_list [lindex $bins 0]
            set dim2_bins_list [lindex $bins 1]

            ;# Binning is controlled by the bead designated in $headnames.
            ;# Creates a dict that contains the bin and leaflet information linked to
            ;# a resid and index number. Facilitates easy binning later. 
            set res_dict [createResidueDictionaries [dict get $config_dict species] [dict get $config_dict headnames] [dict get $sel_info lipid_list] [dict get $sel_info name_list] $dim1_bins_list $dim2_bins_list [dict get $sel_info leaflet_list] $selex]


            ;# Make necessary calculations (if any), then bin and average them
            if {$quantity_of_interest eq "height"} {
                set outfiles [averageHeight $res_dict $outfiles [dict get $sel_info lipid_list] [dict get $sel_info zvals_list]]
            } elseif {$quantity_of_interest eq "tilt_order"} {
                set tilts [fitVecsToSel [dict keys $selections] [dict get $sel_info xvals_list] [dict get $sel_info yvals_list] [dict get $sel_info zvals_list]]
                set orders [calculateOrderParameters [dict keys $selections] [dict get $sel_info xvals_list] [dict get $sel_info yvals_list] [dict get $sel_info zvals_list]]
                set outfiles [averageTiltAndOrderParameter $res_dict $outfiles [dict get $sel_info lipid_list] $tilts $orders [dict get $sel_info tail_list] $selex]
            }

            ;# Now that all information has been binned, print it to files
            dict for {key val} [dict get $outfiles $selex] {
                if {[string match "counts*" $key]} {
                    continue
                }
                printFrame [dict get $bindims N1] $outfiles $key [dict get $bindims d1] [dict get $config_dict min] [dict get $bindims N2] $polar $selex
                
                ;# cleanup before next step
                set outfiles [dict unset outfiles $selex $key bin]
            } 

            ;# cleanup before next step
            dict remove $res_dict [dict keys $res_dict]
            dict remove $sel_info [dict keys $sel_info]
            }
        }
    }
    if {$switch == 0} {
        close $files
    }
    ;# output log info that nougat.py can read later 
    if {$quantity_of_interest eq "height"} {
        outputNougatLog [dict get $config_dict start] [dict get $config_dict nframes] [dict get $config_dict step] [dict get $config_dict species] [dict get $config_dict acyl_names] $system [dict get $config_dict headnames] $coordsys $foldername [dict get $bindims N1] [dict get $bindims N2] [dict get $bindims d1] [dict get $bindims d2]
        
        ;# this will only work if your TMD helices are set to occupancy 1
        ;# all it does is put a dot on the polar heatmap where a TMD helix should be, so not essential at all
        ;# this proc is currently broken, anyways (yes, there is a github issue)
        if {[dict exists $config_dict custom_position] && ([dict get $config_dict custom_position] != "NULL")} {
            puts "[dict get $config_dict custom_position]"
            Protein_Position $system [dict get $config_dict headnames] [dict get $config_dict custom_position] $foldername
        }
    }

    ;# close all outfiles
    foreach channel [file channels "file*"] {
        close $channel
    } 

    ;# delete all atomselections in scope.
    ;# if the user has previously defined atomselections in their VMD session,
    ;# these will be out of scope and cause an error.
    ;# catch will ignore this error, as it's not important to the user.
    foreach selection [atomselect list] {
        catch {$selection delete}
    }
}
