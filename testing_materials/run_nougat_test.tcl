set script_path [pwd]
source $script_path/../nougat.tcl


proc load_and_run_test {trajpath groname xtcname test_name config d1 d2 polar home} {
	cd $trajpath
	mol new $groname
	if {$xtcname != 0} {
		mol addfile $xtcname waitfor all
		animate delete beg 0 end 0 skip 0 top
	}
	start_nougat $test_name $config $d1 $d2 0 -1 1 $polar
	mol delete top
	cd $home
}

;# Put your tests below

;# E protein cartesian
set path "${script_path}/E-protein_trajectory"
load_and_run_test $path DT_test.pdb DT_test.xtc test ${path}/nougat_config_test_E-protein.txt 5 5 0 $script_path

;# E protein polar
set path "${script_path}/E-protein_trajectory"
load_and_run_test $path DT_test.pdb DT_test.xtc test ${path}/nougat_config_test_E-protein.txt 3 12 1 $script_path

;# Molecularly flat membrane cartesian
set path "${script_path}/flat_surface_test"
load_and_run_test $path flat.gro 0 test ${script_path}/E-protein_trajectory/nougat_config_test_E-protein.txt 5 5 0 $script_path

;# Molecularly flat membrane polar
set path "${script_path}/flat_surface_test"
load_and_run_test $path flat.gro 0 test ${script_path}/E-protein_trajectory/nougat_config_test_E-protein.txt 3 12 1 $script_path

;# exit
exit
