#! /bin/sh

# Run unit tests
cd ../test/ 
vmd -dispdev none -e ./run_unit_tests.tcl 2>&1 | tee -a tcl_unit_test.log
cd ../testing_materials/

# Run nougat.tcl and nougat.py on test systems
vmd -dispdev none -e ./run_nougat_test.tcl
bash ./run_nougat_py_test.sh 2>&1 | tee -a nougpy.log

# Run test battery on resulting files; print to file and terminal (tee)
python3 -m pytest tests/ 2>&1 | tee -a pytest.log

# Sound the all clear
echo "Acceptance testing finished"
