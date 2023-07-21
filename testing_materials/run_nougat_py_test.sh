#! /bin/sh/

cd /home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory/test_cart_5_5_0_-1_1
echo "I'm here"
python3 /home/js2746/PolarHeightBinning/plotting/nougat.py test /home/js2746/PolarHeightBinning/plotting/nougat_plot_config.txt

cd /home/js2746/PolarHeightBinning/testing_materials/E-protein_trajectory/test_polar_3_12_0_-1_1
python3 /home/js2746/PolarHeightBinning/plotting/nougat.py test /home/js2746/PolarHeightBinning/plotting/nougat_plot_config.txt -p