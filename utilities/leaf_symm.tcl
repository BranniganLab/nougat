;# A set of procs that will remove lipids from the bulk of one leaflet in order
;# to enforce starting conditions of a particular lipid number asymmetry across
;# membrane leaflets.

source ~/PolarHeightBinning/nougat.tcl

proc count {sel} {
	return [llength [lsort -unique [$sel get resid]]]
}

proc pick_rand {n} {
;# source: https://wiki.tcl-lang.org/page/rand

	return [expr {int(rand() * $n)}]
}

proc cleanup_sels {} {
	foreach selection [atomselect list] {
        catch {$selection delete}
    }
}


proc leaf_sym {desired_asymm} {
	cell_prep ~/PolarHeightBinning/nougat_config.txt 1
	
	;# this step is needed bc cell_prep renumbered the beta fields and will
	;# cause writepdb to fail as a result
	set sel [atomselect top all]
	$sel set beta 1
	$sel delete

	set upper_tot [atomselect top "user 1"]
	set lower_tot [atomselect top "user 2"]
	set delete_list []

	set current_asymm [expr [count $upper_tot] - [count $lower_tot]]
	set number_to_delete [expr $current_asymm - $desired_asymm]

	if {$number_to_delete < 0} {
		set userval 2
	} elseif {$number_to_delete > 0} {
		set userval 1
	} else {
		return
	}

	set bulk [atomselect top "user $userval and not within 20 of name BB"]
	set resids [lsort -unique [$bulk get resid]]

	for {set i 0} {$i < [expr {abs($number_to_delete)}]} {incr i} {
		set num_bulk [llength $resids]
		set ind [pick_rand $num_bulk]
		lappend delete_list [lindex $resids $ind]
		set resids [lreplace $resids $ind $ind]
	}

	puts $delete_list
	set sel [atomselect top "not (user $userval and resid $delete_list)"]
	$sel writepdb new.pdb

	cleanup_sels
}