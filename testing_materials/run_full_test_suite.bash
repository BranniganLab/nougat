#! /bin/bash

# Remove old log files to prevent accidental appending. -f is to suppress
# error if no such file exists.
rm -f nougpy.log 
rm -f pytest.log 
rm -f tcl_unit_test.log 

# Run tcl unit tests and divert output to file
cd ../test/ 
echo "Starting TCL unit testing"
vmd -dispdev none -e ./run_unit_tests.tcl > ../testing_materials/tcl_unit_test.log
cd ../testing_materials/

# Check to make sure VMD didn't error
if grep "couldn't load file" tcl_unit_test.log
then
	echo "Something went wrong in VMD :("
	echo "Your unit tests likely did not run."
	echo "Do you want to continue with acceptance testing anyway?"
	echo "(Choose the number that corresponds to your choice)"
	select yn in "Yes, Continue" "No, Exit"; do
    	case $yn in
        	"Yes, Continue" ) break;;
        	"No, Exit" ) exit;;
    	esac
	done
fi

# Check to make sure tests all passed
if grep FAILED tcl_unit_test.log
then
	echo "TCL unit testing failed :("
	echo "Do you want to continue with acceptance testing anyway?"
	echo "(Choose the number that corresponds to your choice)"
	select yn in "Yes, Continue" "No, Exit"; do
    	case $yn in
        	"Yes, Continue" ) break;;
        	"No, Exit" ) exit;;
    	esac
	done
else
	echo "TCL unit testing passed :)"
fi

# Run python unit tests and divert output to file and terminal.
python3 -m pytest ../test/Unit_Test.py 2>&1 | tee -a pyunittest.log
# Check to make sure tests all passed

# Run nougat.tcl and nougat.py on test systems
vmd -dispdev none -e ./run_nougat_test.tcl
bash ./run_nougat_py_test.sh 2>&1 | tee -a nougpy.log

# Run acceptance tests on resulting files
python3 -m pytest tests/ 2>&1 | tee -a pytest.log

# Sound the all clear
echo "Acceptance testing finished"
