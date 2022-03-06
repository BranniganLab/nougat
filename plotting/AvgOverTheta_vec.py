import matplotlib.pyplot as plt
import numpy as np
import warnings



#name_list = ["PO"]
#name_list = ["DL", "DT", "DG", "DX", "PO", "DB", "DY", "DO", "DP"]
#name_list = ["DL", "DT", "DX", "DB", "DP"]
#name_list = ["DG", "DY", "DO"]
#name_list = ["DP", "DO", "PO"]
#name_list = {"DL","DY"}
name_list = {"DB","../../scripts/PolarHeightBinning/results/DB"}
#field_list = ["zone","ztwo","zplus","zzero"]
field_list = ['zone', 'ztwo']

colordict = {
	"DT": "green",
	"DL": "blue",
	"DX": "purple",
	"DB": "red",
	"DY": "limegreen",
	"DO": "mistyrose",
	"PO": "darkorange",
	"DP": "deepskyblue",
	"DG": "orchid",
	"../../scripts/PolarHeightBinning/results/DB": "deepskyblue"
}

fig = plt.figure(figsize = (5,5))
plt.xlim(0,70)
plt.ylim(-50,20)
plt.gca().set_aspect('equal',adjustable='box')

x = np.arange(3,70,6)

for name in name_list:
	for field in field_list:
		data = np.load(name+"."+field+".avgheight.npy")

		for row in range(data.shape[0]):
			nonzerocount = np.count_nonzero(data[row,:])
			nancount = np.count_nonzero(np.isnan(data[row,:]))
			if (nancount/nonzerocount) >= 0.2:
				data[row,:] = np.nan
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", category=RuntimeWarning)
			z_vals_lower=np.nanmean(data, axis=1)
		
		plt.plot(x,z_vals_lower,color=colordict[name])



plt.show()