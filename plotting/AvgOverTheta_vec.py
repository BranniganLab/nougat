import matplotlib.pyplot as plt
import numpy as np
import warnings
from utils import *


allsys = ["DL", "DT", "DG", "DX", "PO", "DB", "DY", "DO", "DP"]
satsys = ["DL", "DT", "DX", "DB", "DP"]
monounsatsys = ["DG", "DY", "DO"]
fivebeads = ["DB","DG"]
fourbeads = ["DP", "DO"]
threebeads = ["DL","DY"]
custom = ['DX']

sys_list = [custom]
sys_name_list = ['customDX']
#sys_list = [satsys, monounsatsys, fivebeads, fourbeads, threebeads]
#sys_name_list = ["satsys", "monounsatsys", "fivebeads", "fourbeads", "threebeads"]

field_list = ["zone","ztwo"]
#field_list = ['zone', 'ztwo', 'zzero']
#field_list = ['zone', 'ztwo', 'zzero', 'zplus']

colordict = {
	"DT": "red",
	"DL": "orange",
	"DX": "purple",
	"DB": "blue",
	"DY": "orange",
	"DO": "green",
	"PO": "green",
	"DP": "green",
	"DG": "blue"
}

fielddict = {
	"zone" : "solid",
	"ztwo" : "solid",
	"zzero" : "dotted",
	"zplus" : "dashed"
}

bead_dict = {
  "DT" : ['C2A.C2B'],
  "DL" : ['C2A.C2B', 'C3A.C3B'],
  "DY" : ['D2A.D2B', 'C3A.C3B'],
  "DO" : ['D2A.D2B', 'C3A.C3B', 'C4A.C4B'],
  "PO" : ['D2A.C2B', 'C3A.C3B', 'C4A.C4B'],
  "DP" : ['C2A.C2B', 'C3A.C3B', 'C4A.C4B'],
  "DB" : ['C2A.C2B', 'C3A.C3B', 'C4A.C4B', 'C5A.C5B'],
  "DG" : ['C2A.C2B', 'D3A.D3B', 'C4A.C4B', 'C5A.C5B'],
  "DX" : ['C2A.C2B', 'C3A.C3B', 'C4A.C4B', 'C5A.C5B', 'C6A.C6B']
}

beadcolordict = {
	"C2A.C2B" : "grey",
	"D2A.D2B" : "grey",
	"D2A.C2B" : "grey",
	"C3A.C3B" : "magenta",
	"D3A.D3B" : "magenta",
	"C4A.C4B" : "aquamarine",
	"C5A.C5B" : "aquamarine",
	"C6A.C6B" : "sienna"
}

'''
#plot zone, ztwo, zzero, zplus for all systems independently
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


#plot zone, ztwo, zzero, and all variations of zplus
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
	for bead in bead_dict[name]:
		data = np.load(name+"."+bead+".zplus.avgheight.npy")
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", category=RuntimeWarning)
			z_vals=np.nanmean(data, axis=1)
		plt.plot(x,z_vals,color=beadcolordict[bead],linestyle=fielddict[field])

	X = [28.116,28.116]
	Y = [-1,5]
	plt.plot(X,Y,'k:')


	fig.set_size_inches(6,6)
	plt.savefig(name+"_avgovertheta_allbeadszplus.png", dpi = 700)
	plt.clf()
	plt.close()

plot combined systems zone and ztwo
'''
max_scale_dict = {
	"height":10,
	"thickness":1.2,
	"curvature":0.03,
	"Kcurvature":0.002,
	"tail1":1.2,
	"tail0":1.2,
	"density":2
}
min_scale_dict = {
	"height":-30,
	"thickness":0,
	"curvature":-0.03,
	"Kcurvature":-0.002,
	"tail1":-0.2,
	"tail0":-0.2,
	"density":0
}

for measure in ["height", "curvature", "Kcurvature", "thickness", "tail1", "tail0"]:
	nmcount = 0 
	for system in sys_list:
		fig, axs = plt.subplots(2, sharex=True, sharey=True)
		counter = 0
		for field in field_list:
			axs[counter].set_xlim(0,18)
			axs[counter].set_ylim(min_scale_dict[measure],max_scale_dict[measure])
			for name in system:
				data = np.genfromtxt("lg"+name+"/lg"+name+"_polar_5_10_100_-1_1/dat/"+filename_generator("lg"+name, name+"PC", field, "C1A.C1B", "polar", measure, "dat"),delimiter=",",missing_values='nan',filling_values=np.nan)
				#if name == "DT":
				#	data = np.genfromtxt("dm1/lg"+name+"/lg"+name+"_polar_5_10_100_-1_1/dat/"+filename_generator("lg"+name, name+"PC", field, "C1A.C1B", "polar", measure, "dat"),delimiter=",",missing_values='nan',filling_values=np.nan)
				#else:
				#	data = np.genfromtxt("lg"+name+"/lg"+name+"_polar_5_10_100_-1_1/dat/"+filename_generator("lg"+name, name+"PC", field, "C1A.C1B", "polar", measure, "dat"),delimiter=",",missing_values='nan',filling_values=np.nan)
				#for row in range(data.shape[0]):
				#	nonzerocount = np.count_nonzero(data[row,:])
				#	nancount = np.count_nonzero(np.isnan(data[row,:]))
				#	if (nancount/nonzerocount) >= 0.2:
				#		data[row,:] = np.nan
				with warnings.catch_warnings():
					warnings.simplefilter("ignore", category=RuntimeWarning)
					z_vals=np.nanmean(data, axis=1)
				if measure == "height":
					first_val = find_first_val(z_vals)
					z_vals = z_vals - first_val
				#if measure == "tail0" or measure == "tail1":
				#	last_val = find_last_val(z_vals)
				#	z_vals = z_vals/last_val
				maxval = len(z_vals)
				x = np.arange(2.5,(maxval*5+2.5),5) / 10
				axs[counter].plot(x,z_vals,color=colordict[name])
				#X = [28.116,28.116]
				#Y = [-1,5]
				#plt.plot(X,Y,'k:')
			counter += 1

			#fig.set_size_inches(6,6)
		plt.savefig(sys_name_list[nmcount]+"_avg"+measure+"overtheta_combo.pdf", dpi = 700)
		plt.clf()
		plt.close()
		nmcount += 1

'''
#plot all bead variations of zone, ztwo, and zplus, as well as orig zzero
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
		for bead in bead_dict[name]:
			#if field != "zzero":
			data = np.load(name+"."+bead+"."+field+".avgheight.npy")
			with warnings.catch_warnings():
				warnings.simplefilter("ignore", category=RuntimeWarning)
				z_vals=np.nanmean(data, axis=1)
			plt.plot(x,z_vals,color=beadcolordict[bead],linestyle=fielddict[field])

	X = [28.116,28.116]
	Y = [-1,5]
	plt.plot(X,Y,'k:')


	fig.set_size_inches(6,6)
	plt.savefig(name+"_avgovertheta_allbeads.png", dpi = 700)
	plt.clf()
	plt.close()
	'''