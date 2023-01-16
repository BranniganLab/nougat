

proc separate_chains {molid cutoff} {
	;# Multichain proteins default to chain X in martinize.py.
	;# Separate_chains will jump from atom to atom measuring the distance
	;# between beads. 
	;# Anywhere the distance is greater than $cutoff 
	;# it will assume there is a new subunit.
	;# This makes the (large) assumption that consecutive indices will be 
	;# assigned to consecutive residues in the sequence!

	set sel [atomselect $molid "name BB SC1 to SC4"]
    set indxlist [$sel get index]
    set indmin [::tcl::mathfunc::min {*}$indxlist]
    set indmax [::tcl::mathfunc::max {*}$indxlist]
    $sel delete

    set list1 [list $indmin]

    for {set i $indmin} {$i < $indmax} {incr i} {
        set sel1 [atomselect $molid "index $i"]
        set sel2 [atomselect $molid "residue [expr $i+1]"]

        set loc1 [lindex [$sel1 get {x y z}] 0]
        set loc2 [lindex [$sel2 get {x y z}] 0]

        set dist [vecdist $loc1 $loc2]

            if {$dist > $cutoff} {
                lappend list1 $i
                lappend list1 [expr $i+1]
            }
        $sel1 delete
        $sel2 delete
    }

    lappend list1 $indmax

    set chars [list A B C D E]
    set k 0

    foreach {i j} $list1 {
        set sel [atomselect $molid "index $i to $j"]
        $sel set chain [lindex $chars $k]
        $sel delete
        set k [expr $k+1]
    }
}


proc contacts_by_residue {molid chain lipid beadname cutoff start step} {

	set sel [atomselect $molid "name BB and chain $chain"]
	set residlist [$sel get resid]
	set resnamelist [$sel get resname]
	$sel delete

	set outfile [open "tcl_output/contacts_and_heights_chain_${chain}.dat" w]

	foreach resid $residlist resname $resnamelist {
		set contactnumlist []
		set lipidsel [atomselect $molid "resname $lipid and name $beadname and (within $cutoff of (chain $chain and resid $resid))"]
		for {set frm $start} {$frm < [molinfo $molid get numframes]} {set frm [expr $frm + $step]} {
			$lipidsel frame $frm
			$lipidsel update
			lappend contactnumlist [llength [lsort -unique [$lipidsel get resid]]]
		}
		puts $outfile "$resid $resname $contactnumlist"
		$lipidsel delete
		puts "$resid finished"
	}

	close $outfile
	return
}


proc contact_mapper {molid lipid cutoff start step} {
	separate_chains $molid 25

	set sel [atomselect $molid "name BB"]
	set chain_list [lsort -unique [$sel get chain]]
	$sel delete

	foreach chain $chain_list {
		contacts_by_residue $molid $chain $lipid "C1A C1B" $cutoff $start $step
	}
}