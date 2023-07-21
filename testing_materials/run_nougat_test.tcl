source ~/PolarHeightBinning/nougat.tcl

proc load_and_run_test {trajpath groname xtcname test_name config d1 d2 polar} {
	cd $trajpath
	mol new $groname
	if {$xtcname != 0} {
		mol addfile $xtcname waitfor all
		animate delete beg 0 end 0 skip 0 top
	}
	start_nougat $test_name $config $d1 $d2 0 -1 1 $polar
	mol delete top
}

;# Put your tests below

;# E protein cartesian
set path "/home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory"
load_and_run_test $path DT_test.gro DT_test.xtc E-protein /home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory/nougat_config_test_E-protein.txt 5 5 0

;# E protein polar
set path "/home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory"
load_and_run_test $path DT_test.gro DT_test.xtc E-protein /home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory/nougat_config_test_E-protein.txt 3 12 1

;# Molecularly flat membrane cartesian
set path "/home/js2746/PolarHeightBinning/testing_materials/flat_surface_test"
load_and_run_test $path flat.gro 0 flat /home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory/nougat_config_test_E-protein.txt 5 5 0

;# Molecularly flat membrane polar
set path "/home/js2746/PolarHeightBinning/testing_materials/flat_surface_test"
load_and_run_test $path flat.gro 0 flat /home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory/nougat_config_test_E-protein.txt 3 12 1

;# exit
exit