#! /bin/sh/

git rm -r ./E-protein_trajectory/E-protein_*
git rm -r ./flat_surface_test/flat_*

git commit -m 'deleting old acceptance testing ref values'

vmd -dispdev none -e ./regen_acc_test_data.tcl

# Add any systems you want nougat to run below

python3 regen_acc_test_data.py ./E-protein_trajectory/E-protein_cart_5_5_0_-1_1
python3 regen_acc_test_data.py ./E-protein_trajectory/E-protein_polar_3_12_0_-1_1 -p
python3 regen_acc_test_data.py ./flat_surface_test/flat_cart_5_5_0_-1_1
python3 regen_acc_test_data.py ./flat_surface_test/flat_polar_3_12_0_-1_1 -p

git add ./E-protein_trajectory/E-protein_*
git add ./flat_surface_test/flat_*

git commit -m 'generating new acceptance testing ref values'
