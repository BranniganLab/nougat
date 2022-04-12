set testpath "lgPO.D2A.C2B.zone.avgcurvature.dat"

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

proc color_by_user {path} {
	set data [read_in $path]

	set row [llength $data]
	set col [llength [lindex $data 0]]



	for {set i 0} {$i < $row} {incr i} {
		for {set j 0} {$j < $col} {incr j} {
			if {[lindex [lindex $data $i] $j] ne "nan"} {
				set sel [atomselect top "resname C2 and segname zone and occupancy $i and beta $j"]
				$sel set user2 [lindex [lindex $data $i] $j]
				$sel delete
			}
		}
	}
}

color_by_user $testpath