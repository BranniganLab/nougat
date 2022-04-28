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

proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.6 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.3
    graphics $mol cone $middle $end radius 0.6
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

proc draw_tilt_arrows {data bead field} {
	set row [llength $data]
	set col [expr [llength [lindex $data 0]] / 3]
	for {set i 0} {$i < $row} {incr i} {
		for {set j 0} {$j < $col} {incr j} {
			set indx1 [expr $j*3]
			if {[lindex [lindex $data $i] $indx1] ne "nan"} {
				set indx2 [expr $indx1 + 1]
				set indx3 [expr $indx2 + 1]
				set sel [atomselect top "resname $bead and segname $field and occupancy $i and beta $j"]
				set vecxyz "[lindex [lindex $data $i] $indx1] [lindex [lindex $data $i] $indx2] [lindex [lindex $data $i] $indx3]"
				set selxyz "[$sel get x] [$sel get y] [$sel get z]"
				set arrow_base [vecexpr $selxyz $vecxyz add]
				vmd_draw_arrow top $selxyz $arrow_base
				$sel delete
			}
		}
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
				set filename $system/$system.$bead.$field.$qual.dat 
				set data [read_in $filename]
				colorize $data $beadnm $field $qual
			}
		}
	}
	lappend field_list "zzero"
	foreach field $field_list {
		foreach qual $qual_list {
			set filename $system/$system.$field.$qual.dat 
			set data [read_in $filename]
			colorize $data C1 $field $qual 
		}
	}
	foreach field "zone ztwo" {
		set filename $system/$system.$field.avgtilt.dat 
		set data [read_in $filename]
		draw_tilt_arrows $data C1 $field 
	}
}






color_by_user "lgPO"