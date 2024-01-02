#! /bin/sh/

# Add any systems you want nougat to run below

cd ./E-protein_trajectory/test_cart_5_5_0_-1_1
python3 ../../../plotting/nougat.py
cd ../..

cd ./E-protein_trajectory/test_polar_3_12_0_-1_1
python3 ../../../plotting/nougat.py -p
cd ../..

cd ./flat_surface_test/test_cart_5_5_0_-1_1
python3 ../../../plotting/nougat.py
cd ../..

cd ./flat_surface_test/test_polar_3_12_0_-1_1
python3 ../../../plotting/nougat.py -p
cd ../..
