#! /bin/sh

# Run nougat.tcl and nougat.py on test systems
vmd -e /home/js2746/PolarHeightBinning/testing_materials/run_nougat_test.tcl
bash /home/js2746/PolarHeightBinning/testing_materials/run_nougat_py_test.sh

# Run test battery on resulting files
python3 -m pytest tests/

# Clean up if all tests passed

#if [ "$cleanup" == "True" ]; then
#	rm -r ~/PolarHeightBinning/testing_materials/E-protein_trajectory/test_*
#	rm -r ~/PolarHeightBinning/testing_materials/flat_surface_test/test_*
#fi

# Sound the all clear
echo "Acceptance testing finished"