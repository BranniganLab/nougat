# readFile
#
# Reads file contents and splits them into individual lines.
#
# Arguments:
#   path: the path to the file you wish to read in.
#
# Returns:
#   List containing each line of the file.
#
# Example:
#   set data [readFile ~/Desktop/sample_file.txt]

proc readFile {path} {
	set f [open $path r]
	set file_data [read $f]
	close $f
	set lines [split $file_data "\n"]
	if {[lindex $lines end] eq ""} {
 	   set lines [lreplace $lines end end]
	}
	return $lines
}


# drawTriangles
#
# Reads in triangles file from $path and draws them on molecule number $molno.
# User must have used make_surface_triangles in nougat.py.
#
# Arguments:
#   path: the path to the output of make_surface_triangles.py
#   molno: the molecule number in vmd these triangles should be added to.
#   color: [OPTIONAL] draw all triangles as this color. Default is red.
#   material: [OPTIONAL] draw all triangles as this material. Default is Diffuse.
#
# Returns:
#   Nothing
#
# Example:
#   drawTriangles ~/Desktop/sample_file.txt top green AOChalky

proc drawTriangles {path molno {color "red"} {material "Diffuse"}} {
	set triangles [readFile $path]
	draw color $color
	draw material $material
	for {set line 0} {$line < [llength $triangles]} {set line [expr $line + 3]} {
		set coord1 [lindex $triangles $line]
		set coord2 [lindex $triangles [expr $line + 1]]
		set coord3 [lindex $triangles [expr $line + 2]]
		graphics $molno triangle $coord1 $coord2 $coord3
	}
}