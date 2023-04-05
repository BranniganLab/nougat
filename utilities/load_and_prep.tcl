proc load_and_prep {sysname trajpath configpath groname xtcname first last step d1 d2 polar} {
	set gro "${trajpath}/${groname}"
	set xtc "${trajpath}/${xtcname}"
	mol new $gro
	mol addfile $xtc first ${first} last ${last} step ${step} waitfor all
	animate delete beg 0 end 0 skip 0 top
	start_nougat $sysname $configpath $d1 $d2 0 -1 1 $polar
	mol delete top
}
