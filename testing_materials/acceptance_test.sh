#! /bin/sh

# Run unit tests
# cd ../test/ 
# vmd -dispdev none -e ./run.tcl > unit_tests.log
# cd ../testing_materials/

# Run nougat.tcl and nougat.py on test systems
vmd -dispdev none -e ./run_nougat_test.tcl
bash ./run_nougat_py_test.sh

# Run test battery on resulting files; print to file and terminal (tee)
python3 -m pytest tests/ 2>&1 | tee -a pytest.log

# Clean up
rm -r ./E-protein_trajectory/test_*
rm -r ./flat_surface_test/test_*

# Sound the all clear
echo "Acceptance testing finished"