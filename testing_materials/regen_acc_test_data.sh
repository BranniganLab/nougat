#! /bin/sh/

vmd -dispdev none -e ./regen_acc_test_data.tcl

# Add any systems you want nougat to run below

cd ./E-protein_trajectory/E-protein_cart_5_5_0_-1_1
python3 ../../../plotting/nougat.py E-protein ../../../plotting/nougat_plot_config.txt
cd ../..

cd ./E-protein_trajectory/E-protein_polar_3_12_0_-1_1
python3 ../../../plotting/nougat.py E-protein ../../../plotting/nougat_plot_config.txt -p
cd ../..

cd ./flat_surface_test/flat_cart_5_5_0_-1_1
python3 ../../../plotting/nougat.py flat ../../../plotting/nougat_plot_config.txt
cd ../..

cd ./flat_surface_test/flat_polar_3_12_0_-1_1
python3 ../../../plotting/nougat.py flat ../../../plotting/nougat_plot_config.txt -p
cd ../..