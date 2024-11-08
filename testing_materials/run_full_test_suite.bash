#! /bin/bash

# Allow bash aliases (specifically 'vmd')
shopt -s expand_aliases
source ~/.bashrc

# Remove old log files to prevent accidental appending. -f is to suppress
# error if no such file exists.
rm -f nougpy.log
rm -f pyunittest.log
rm -f pytest.log
rm -f tcl_unit_test.log
rm -f nougat_test_outputs.log

# Run tcl unit tests and divert output to file
echo "Starting TCL unit testing"
vmd -dispdev none -eofexit < tcltest/unit_test.test > ./tcl_unit_test.log

# Check to make sure VMD exists
if [ ! -s tcl_unit_test.log ]
then
	echo "VMD is not recognized or not installed."
	exit
fi

# Check to make sure VMD didn't error
if grep -E 'invalid \command|bad \option|nt \load \file|unknown \option' tcl_unit_test.log
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
echo "Starting pytest unit testing"
python3 -m pytest tcltest/Unit_Test.py 2>&1 | tee -a pyunittest.log

# Check to make sure testing did not fail
if grep -q "No module named" pyunittest.log
then
	echo "pytest unit testing failed :("
	echo "You need to install pytest before proceeding"
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
if grep -q FAILED pyunittest.log
then
	echo "pytest unit testing failed :("
	echo "Do you want to continue with acceptance testing anyway?"
	echo "(Choose the number that corresponds to your choice)"
	select yn in "Yes, Continue" "No, Exit"; do
    	case $yn in
        	"Yes, Continue" ) break;;
        	"No, Exit" ) exit;;
    	esac
	done
else
	echo "Pytest unit testing passed :)"
fi

# Run nougat.tcl and nougat.py on test systems
echo "Starting acceptance testing"
vmd -dispdev none -e ./run_nougat_test.tcl > nougat_test_outputs.log

# Check to make sure testing did not fail
if grep -q "can't" nougat_test_outputs.log
then
        echo "VMD didn't run the acceptance testing properly and it failed :("
        echo "Do you want to continue with acceptance testing anyway?"
        echo "(Choose the number that corresponds to your choice)"
        select yn in "Yes, Continue" "No, Exit"; do
        case $yn in
                "Yes, Continue" ) break;;
                "No, Exit" ) exit;;
        esac
        done
fi


# Run acceptance tests on resulting files
python3 -m pytest -rx tests/ 2>&1 | tee -a pytest.log

# Sound the all clear
echo "Acceptance testing finished"
