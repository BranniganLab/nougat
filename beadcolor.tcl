set scriptpath "~/PolarHeightBinning"
source $scriptpath/polarH3.tcl


proc read_in {path} {
	set fp [open $path r]
	set file_data [read $fp]
	close $fp
	
	set file_data [string map {" " ""} $file_data]
	
	set data_out []

	foreach row $file_data {
		set splitrow [split $row ","]
		lappend data_out $splitrow
	}

	return $data_out
}

proc bead_fetcher {system} {

	set syslist {DT DL DY DO PO lgPO DP DB DG DX}
	set beadlist {{C2A.C2B} {C2A.C2B C3A.C3B} {D2A.D2B C3A.C3B} {D2A.D2B C3A.C3B C4A.C4B} {D2A.C2B C3A.C3B C4A.C4B} {D2A.C2B C3A.C3B C4A.C4B} {C2A.C2B C3A.C3B C4A.C4B} {C2A.C2B C3A.C3B C4A.C4B C5A.C5B} {C2A.C2B D3A.D3B C4A.C4B C5A.C5B} {C2A.C2B C3A.C3B C4A.C4B C5A.C5B C6A.C6B}}

	set indx [lsearch $syslist $system]
	if {$indx == -1} {
		puts "Something is wrong with your system selection"
		return
	} else {
		return [lindex $beadlist $indx]
	}
}

proc colorize {data bead field qual} {
	set row [llength $data]
	set col [llength [lindex $data 0]]
	if {[string equal $field "zzero"]} {
		set fieldnm "zzer"
	} elseif {[string equal $field "zplus"]} {
		set fieldnm "zplu"
	} else {
		set fieldnm $field
	}
	for {set i 0} {$i < $row} {incr i} {
		for {set j 0} {$j < $col} {incr j} {
			if {[lindex [lindex $data $i] $j] ne "nan"} {
				set sel [atomselect top "resname $bead and segname $fieldnm and occupancy $i and beta $j"]
				if {[string equal $qual "avgcurvature"]} {
					$sel set user2 [lindex [lindex $data $i] $j]
				} elseif {[string equal $qual "avgKcurvature"]} {
					$sel set user3 [lindex [lindex $data $i] $j]
				} elseif {[string equal $qual "avgdensity"]} {
					$sel set user4 [lindex [lindex $data $i] $j]
				} else {
					puts "something is wrong"
					return
				}
				$sel delete
			}
		}
	}
}

proc color_by_user {system} {
	set field_list {zone ztwo zplus}
	set bead_list [bead_fetcher $system]
	set qual_list {avgcurvature avgKcurvature}
	foreach bead $bead_list {
		set beadnum [string index $bead 1]
		set beadnm C$beadnum
		foreach field $field_list {
			foreach qual $qual_list {
				set filename $system.$bead.$field.$qual.dat 
				set data [read_in $filename]
				colorize $data $beadnm $field $qual
			}
		}
	}
	lappend field_list "zzero"
	foreach field $field_list {
		foreach qual $qual_list {
			set filename $system.$field.$qual.dat 
			puts $filename
			set data [read_in $filename]
			colorize $data C1 $field $qual 
		}
	}
}






color_by_user "lgPO"