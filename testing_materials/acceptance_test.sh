#! /bin/sh

# Run nougat.tcl and nougat.py on test systems
vmd -e /home/js2746/PolarHeightBinning/testing_materials/run_nougat_test.tcl
bash /home/js2746/PolarHeightBinning/testing_materials/run_nougat_py_test.sh

# Run test battery on resulting files
cd /home/js2746/PolarHeightBinning/testing_materials
python3 -m pytest tests/

# Clean up if all tests passed
rm -r ~/PolarHeightBinning/testing_materials/E-protein_trajectory/test_*
rm -r ~/PolarHeightBinning/testing_materials/flat_surface_test/test_*

# Sound the all clear
echo "Acceptance testing finished"