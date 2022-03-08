import matplotlib.pyplot as plt
import numpy as np
import warnings



allsys = ["DL", "DT", "DG", "DX", "PO", "DB", "DY", "DO", "DP"]
satsys = ["DL", "DT", "DX", "DB", "DP"]
monounsatsys = ["DG", "DY", "DO"]
fivebeads = ["DB","DG"]
fourbeads = ["DP", "DO", "PO"]
threebeads = {"DL","DY"}

sys_list = [allsys, satsys, monounsatsys, fivebeads, fourbeads, threebeads]
sys_name_list = ["allsys", "satsys", "monounsatsys", "fivebeads", "fourbeads", "threebeads"]

#field_list = ['zone', 'ztwo']
field_list = ['zone', 'ztwo', 'zzero', 'zplus']

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

fielddict = {
	"zone" : "solid",
	"ztwo" : "solid",
	"zzero" : "dotted",
	"zplus" : "dashed"
}


for name in allsys:
	fig = plt.figure(figsize = (5,5))
	plt.xlim(0,70)
	plt.ylim(-50,20)
	plt.gca().set_aspect('equal',adjustable='box')

	x = np.arange(3,70,6)
	for field in field_list:
		data = np.load(name+"."+field+".avgheight.npy")
		#for row in range(data.shape[0]):
		#	nonzerocount = np.count_nonzero(data[row,:])
		#	nancount = np.count_nonzero(np.isnan(data[row,:]))
		#	if (nancount/nonzerocount) >= 0.2:
		#		data[row,:] = np.nan
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", category=RuntimeWarning)
			z_vals=np.nanmean(data, axis=1)
		
		plt.plot(x,z_vals,color=colordict[name],linestyle=fielddict[field])
		X = [28.116,28.116]
		Y = [-1,5]
		plt.plot(X,Y,'k:')


	fig.set_size_inches(6,6)
	plt.savefig(name+"_avgovertheta.png", dpi = 700)
	plt.clf()
	plt.close()

'''
for system in sys_list:
	fig = plt.figure(figsize = (5,5))
	plt.xlim(0,70)
	plt.ylim(-50,20)
	plt.gca().set_aspect('equal',adjustable='box')

	x = np.arange(3,70,6)

	for name in system:
		for field in field_list:
			data = np.load(name+"."+field+".avgheight.npy")
			#for row in range(data.shape[0]):
			#	nonzerocount = np.count_nonzero(data[row,:])
			#	nancount = np.count_nonzero(np.isnan(data[row,:]))
			#	if (nancount/nonzerocount) >= 0.2:
			#		data[row,:] = np.nan
			with warnings.catch_warnings():
				warnings.simplefilter("ignore", category=RuntimeWarning)
				z_vals=np.nanmean(data, axis=1)
			
			plt.plot(x,z_vals,color=colordict[name])
			X = [28.116,28.116]
			Y = [-1,5]
			plt.plot(X,Y,'k:')


	fig.set_size_inches(6,6)
	plt.savefig(sys_name_list[counter]+"_avgovertheta.png", dpi = 700)
	plt.clf()
	plt.close()
	counter += 1
'''