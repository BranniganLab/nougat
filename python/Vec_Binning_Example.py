from vecbinning import Multibinner, convert_cartesian_to_polar, convert_radians_and_degrees
DATA = [[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15]]
BINS = [5,10,6]
THETAS = [None,"theta",None]

DATA[0],DATA[1]=convert_cartesian_to_polar(DATA[0], convert_radians_and_degrees(DATA[1]))
Sample_data = Multibinner(DATA, BINS, theta=THETAS)
Sample_data.do_multi_binning()
print(Sample_data.multi_bin_dictionary)
print(Sample_data.raw_data[1].bins)