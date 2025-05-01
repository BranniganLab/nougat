from vecbinning import SimultaneousXYZBinner, convert_cartesian_to_polar, convert_radians_and_degrees
DATA = [[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15]]
BINS = [40,10,6]
THETAS = [None,"theta",None]

DATA[0],DATA[1]=convert_cartesian_to_polar(DATA[0], convert_radians_and_degrees(DATA[1]))
Sample_data = SimultaneousXYZBinner(DATA, BINS, theta=THETAS)

Sample_data.do_xyz_binning()
print(Sample_data.multi_bin_data[0])
print(Sample_data.xbin_indicies)
print(Sample_data.xbins)