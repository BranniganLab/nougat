#! /bin/sh

# Run nougat.tcl and nougat.py on test systems
vmd -e /home/js2746/PolarHeightBinning/testing_materials/run_nougat_test.tcl
bash /home/js2746/PolarHeightBinning/testing_materials/run_nougat_py_test.sh

# Compare new test results to saved test results
result=`python3 ~/PolarHeightBinning/testing_materials/compare_test_results.py`

# Clean up if all tests passed
if [ "$result" == "True" ]; then
	rm -r ~/PolarHeightBinning/testing_materials/E-protein_trajectory/test_*
	rm -r ~/PolarHeightBinning/testing_materials/flat_surface_test/test_*
fi

# Sound the all clear
echo "Acceptance testing finished"