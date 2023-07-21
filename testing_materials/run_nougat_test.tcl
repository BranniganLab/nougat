source ~/PolarHeightBinning/nougat.tcl

proc load_and_run_test {trajpath groname xtcname test_name config d1 d2 polar} {
	set gro "${trajpath}/${groname}"
	set xtc "${trajpath}/${xtcname}"
	mol new $gro
	mol addfile $xtc waitfor all
	animate delete beg 0 end 0 skip 0 top
	start_nougat $test_name $config $d1 $d2 0 -1 1 $polar
	mol delete top
}

;# Put your tests below

;# E protein cartesian
load_and_run_test ~/PolarHeightBinning/testing_materials/E-protein_trajectory DT_test.gro DT_test.xtc E-protein ~/PolarHeightBinning/testing_materials/E-protein_trajectory/nougat_config_test_E-protein.txt 5 5 0

;# E protein polar
load_and_run_test ~/PolarHeightBinning/testing_materials/E-protein_trajectory DT_test.gro DT_test.xtc E-protein ~/PolarHeightBinning/testing_materials/E-protein_trajectory/nougat_config_test_E-protein.txt 3 12 1

;# Molecularly flat membrane cartesian



;# Molecularly flat membrane polar