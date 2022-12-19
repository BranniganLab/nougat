proc separate_chains {molid cutoff} {
	;# Multichain proteins default to chain X in martinize.py.
	;# Separate_chains will jump from atom to atom and anywhere the 
	;# distance is greater than $cutoff it will assume there is a new subunit.
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

    set chars [list A B C D E]
    set k 0
    lappend list1 $indmax

    foreach {i j} $list1 {
        set sel [atomselect $molid "index $i to $j"]
        $sel set chain [lindex $chars $k]
        $sel delete
        set k [expr $k+1]
    }
}


proc contacts_by_residue {residlist molid lipid beadname cutoff chain frame} {
	set resnamelist []
	set residlist []
	set contactnumlist []
	set heightlist []

	foreach resid $residlist {

		set incsel [atomselect $molid "name BB and resid $resid" frame $frame]
		set lipidsel [atomselect $molid "name $lipid and name $beadname and within $cutoff of (chain $chain and resid $resid)" frame $frame]

		lappend resnamelist [$incsel get resname]
		lappend residlist $resid 
		lappend contactnumlist [$lipidsel get num]
		lappend heightlist [$incsel get z]

		$incsel delete
		$lipidsel delete
	}
}