#! /bin/sh

# Run nougat.tcl and nougat.py on test systems
vmd -dispdev none -e ./run_nougat_test.tcl
bash ./run_nougat_py_test.sh

# Run test battery on resulting files
python3 -m pytest tests/

# Clean up if all tests passed
rm -r ./E-protein_trajectory/test_*
rm -r ./flat_surface_test/test_*

# Sound the all clear
echo "Acceptance testing finished"