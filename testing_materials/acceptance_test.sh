#! /bin/sh

# Run nougat.tcl and nougat.py on test system(s)
vmd -e ~/PolarHeightBinning/testing_materials/run_nougat_test.tcl
python3 ~/PolarHeightBinning/testing_materials/run_nougat_test.py

# Compare new test results to saved test results
result='python3 ~/PolarHeightBinning/testing_materials/compare_test_results.py'

# Clean up if all tests passed
if [ "$result" == "True" ]; then
	rm -r ~/PolarHeightBinning/testing_materials/... 
fi

# Sound the all clear
echo "Acceptance testing finished"