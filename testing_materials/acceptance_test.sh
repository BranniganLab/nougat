#! /bin/sh

# Remove old log files to prevent accidental appending. || True is to suppress
# error if no such file exists.
rm nougpy.log || True
rm pytest.log || True
rm tcl_unit_test.log || True

# Run tcl unit tests and divert output to file AND to terminal (tee)
cd ../test/ 
vmd -dispdev none -e ./run_unit_tests.tcl 2>&1 | tee -a ../testing_materials/tcl_unit_test.log
cd ../testing_materials/

# Run python unit tests and divert output to file and terminal.

# Run nougat.tcl and nougat.py on test systems
vmd -dispdev none -e ./run_nougat_test.tcl
bash ./run_nougat_py_test.sh 2>&1 | tee -a nougpy.log

# Run acceptance tests on resulting files
python3 -m pytest tests/ 2>&1 | tee -a pytest.log

# Sound the all clear
echo "Acceptance testing finished"
