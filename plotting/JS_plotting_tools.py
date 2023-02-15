import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 
from thickness import measure_t0
from utils import *

sys_dict = {
	"DT": "green",
	"lgDT" : "green",
	"DL": "blue",
	"lgDL": "blue",
	"DX": "purple",
	"lgDX": "purple",
	"DB": "red",
	"lgDB": "red",
	"DY": "limegreen",
	"lgDY" : "limegreen",
	"DO": "mistyrose",
	"lgDO": "mistyrose",
	"PO": "darkorange",
	"lgPO" : "darkorange",
	"DP": "deepskyblue",
	"lgDP": "deepskyblue",
	"DG": "orchid",
	"lgDG" : "orchid"
}

box_dict = {
	"large": "solid",
	"small" : "dotted"
}

field_dict = {
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

small_sys_list = ["DT", "DY", "DL", "DO", "DP", "PO", "DG", "DB", "DX"]
large_sys_list = ["lgDT", "lgDY", "lgDL", "lgDO", "lgDP", "lgPO", "lgDG", "lgDB", "lgDX"]

def plot_area_per_lipid(systems):
	
	lenlist = []
	fig = plt.figure(figsize = (5,5))
	for system in systems:
		data = np.loadtxt(system+'/'+system+'_boxarea.dat')
		numframes = len(data)-1
		lenlist.append(numframes)

	max_value = max(lenlist)
	X = np.linspace(0,max_value/100.0,max_value)


	for system in systems:
		
		data = np.loadtxt(system+'/'+system+'_boxarea.dat')
		areadata = data[1:]
		numlipidsperleaflet = data[0:1]
		apldata = areadata/(numlipidsperleaflet*1.0)
		numframes = len(areadata)
		window_size = 20 
		rollingavg = []
		for i in range(numframes-window_size):
			windowavg = np.sum(apldata[i:i+window_size])/window_size
			rollingavg.append(windowavg)
		rollingavg = np.array(rollingavg)
		fill_values = max_value-len(rollingavg)
		if fill_values > 0:
			fill_array = np.empty(fill_values)
			fill_array[:] = np.nan
			avgapldata = np.append(rollingavg,fill_array)
		elif fill_values < 0:
			print(system)
		
		plt.plot(X,avgapldata,color=colordict[system],linestyle=linetypedict[system])

	plt.show()

def zoom_in(systems):
	max_scale_dict = {
		"height":60,
		"thickness":2,
		"curvature":0.01,
		"Kcurvature":0.001,
		"tail1":0.6,
		"tail0":0.6,
		"density":2
	}
	min_scale_dict = {
		"height":-60,
		"thickness":0,
		"curvature":-0.01,
		"Kcurvature":-0.001,
		"tail1":0.0,
		"tail0":0.0,
		"density":0
	}
	values = ["height", "thickness", "curvature", "Kcurvature", "tail1", "tail0"]
	#values = ["height"]
	for system in systems:
		#os.chdir("dm1/"+system)
		os.chdir(system)
		os.chdir(system+'_polar_5_10_100_-1_1/')
		for field in ["zone", "ztwo"]:
			for value in values:
				print(value)
				if value == "tail0" or value == "tail1":
					data = np.genfromtxt('dat/'+system+'.'+system[2:]+'PC.'+value+'.'+field+'.polar.avgOrder.dat', delimiter=",", missing_values="nan", filling_values=np.nan)
				else:
					data = np.genfromtxt('dat/'+system+'.'+field+'.C1A.C1B.polar.avg'+value+'.dat', delimiter=",", missing_values="nan", filling_values=np.nan)
				dim1 = np.linspace(0,50,16)
				dim2 = np.linspace(0,2*np.pi,11)
				dim1vals,dim2vals=np.meshgrid(dim1, dim2, indexing='ij')
				print(data[:15,:].shape)
				plot_maker(dim1vals, dim2vals, data[:15,:], system, field, max_scale_dict[value], min_scale_dict[value], False, value, False, "polar")
		os.chdir('../..')



def diff_mid_interface(systems, mol, coordsys):
	for system in systems:
		os.chdir(system)
		filename_start = '/u1/home/js2746/Bending/PC/whole_mols/'+mol+'/dm1/'+system+'/'+system+'_polar_5_10_100_-1_1/npy/'+system+'.'
		filename_end = '.C1A.C1B.'+coordsys+'.height.npy'
		fig = plt.figure()
		ax = plt.subplot()
		zzero = np.load(filename_start+'zzero'+filename_end)
		zone = np.load(filename_start+'zone'+filename_end)
		ztwo = np.load(filename_start+'ztwo'+filename_end)
#		H1 = np.load(filename_start+'zone.C1A.C1B.'+coordsys+'.meancurvature.npy')
#		H2 = np.load(filename_start+'ztwo.C1A.C1B.'+coordsys+'.meancurvature.npy')
		zplus = np.load(filename_start+'zplus'+filename_end)
		diff = zplus-zzero
		avgdiff = np.nanmean(diff,axis=2)
#		H = H1+H2
#		avgH = np.nanmean(H,axis=2)
		t0 = measure_t0(zone, ztwo, coordsys)
		avgdiff = avgdiff/t0
#		avgcombo = avgdiff*avgH

		dims = bin_prep(system, "C1A.C1B", coordsys, "OFF")
		N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims
		plot_maker(dim1vals, dim2vals, avgdiff, system, 'comb', .1, -.1, False, "avgEpsilon", False, coordsys)
		np.save('/u1/home/js2746/Bending/PC/whole_mols/'+mol+'/dm1/'+system+'/npy/'+sys_name+'.epsilon_t0.npy',avgdiff)
		os.chdir('..')
#		c = plt.pcolormesh(dim1vals,dim2vals,avgdiff,cmap="RdBu_r",zorder=0,vmax=.001,vmin=-.001)
#		cbar = plt.colorbar(c)
#		plt.axis('off')
#		ax.set_xticklabels([])
#		ax.set_yticklabels([])
#		#fig.set_size_inches(6,6)
#		plt.savefig(system+"_epsilon_cart.pdf", dpi = 700)
#		plt.clf()
#		plt.close()


def measure_H_epsilon_corr(systems, mol):
	for system in systems:
		filename_start = '/u1/home/js2746/Bending/PC/whole_mols/'+mol+'/'+system+'/npy/'+system+'.'
		filename_end = '.C1A.C1B.cart.height.npy'
		fig = plt.figure()
		ax = plt.subplot()
		try:
			zzero = np.load(filename_start+'zzero'+filename_end)
			H1 = np.load(filename_start+'zone'+'.C1A.C1B.cart.meancurvature.npy')
			H2 = np.load(filename_start+'ztwo'+'.C1A.C1B.cart.meancurvature.npy')
		except:
			filename_end = '.C1A.C1B.D1A.D1B.cart.height.npy'
			zzero = np.load(filename_start+'zzero'+filename_end)
			H1 = np.load(filename_start+'zone'+'.C1A.C1B.D1A.D1B.cart.meancurvature.npy')
			H2 = np.load(filename_start+'ztwo'+'.C1A.C1B.D1A.D1B.cart.meancurvature.npy')
		dims = np.shape(zzero)
		dim1 = np.linspace(0,dims[0],dims[0]+1)
		dim2 = np.linspace(0,dims[1],dims[1]+1)
		dim1vals,dim2vals=np.meshgrid(dim1, dim2, indexing='ij')
		zplus = np.load(filename_start+'zplus'+filename_end)
		epsilon = zplus-zzero
		avgepsilon = np.nanmean(epsilon,axis=2)
		H = H1+H2
		avgH = np.nanmean(H,axis=2)
		t0 = measure_t0(system, mol)
		avgepsilon = avgepsilon/t0
		avgcombo = avgepsilon*avgH
		together = epsilon*H 
		avgtogether = np.nanmean(together,axis=2)
		correlation = avgtogether - avgcombo
		c = plt.pcolormesh(dim1vals,dim2vals,correlation,cmap="RdBu_r",zorder=0,vmax=.15,vmin=-.15)
		cbar = plt.colorbar(c)
		plt.axis('off')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		#fig.set_size_inches(6,6)
		plt.savefig(system+"_epsilonH_correlation_cart.pdf", dpi = 700)
		plt.clf()
		plt.close()

def sum_terms(systems, mol):
	for system in systems:
		filename_start = '/u1/home/js2746/Bending/PC/whole_mols/'+mol+'/'+system+'/npy/'+system+'.'
		filename_end = '.C1A.C1B.cart.height.npy'
		try:
			zzero = np.load(filename_start+'zzero'+filename_end)
			H1 = np.load(filename_start+'zone'+'.C1A.C1B.cart.meancurvature.npy')
			H2 = np.load(filename_start+'ztwo'+'.C1A.C1B.cart.meancurvature.npy')
		except:
			filename_end = '.C1A.C1B.D1A.D1B.cart.height.npy'
			zzero = np.load(filename_start+'zzero'+filename_end)
			H1 = np.load(filename_start+'zone'+'.C1A.C1B.D1A.D1B.cart.meancurvature.npy')
			H2 = np.load(filename_start+'ztwo'+'.C1A.C1B.D1A.D1B.cart.meancurvature.npy')
		zplus = np.load(filename_start+'zplus'+filename_end)
		epsilon = zplus-zzero
		epsilon2 = epsilon*epsilon 
		H = H1+H2
		H2 = H*H
		avgH2 = np.nanmean(H2,axis=2)
		sumavgH2 = np.nansum(avgH2)
		t0 = measure_t0(system, mol)
		print(t0)
		term6 = 2*H*epsilon/t0 
		term6 = np.nanmean(term6,axis=2)
		term7 = epsilon2/(2*t0**2) 
		term7 = np.nanmean(term7,axis=2)
		sumterm6 = np.nansum(term6)
		sumterm7 = np.nansum(term7)
		print(sumterm6)
		print(sumterm7)
		print(sumavgH2)

def sum_over_K(systems):

	for system in systems:
		fig,ax = plt.subplots()
		for field in ["zone", "ztwo"]:
			data = np.load(system+'/npy/'+system+'.'+field+'.C1A.C1B.cart.gausscurvature.npy')
			k_list = []
			for frame in range(np.shape(data)[2]):
				frame_data = data[:,:,frame]
				k_sum = np.nansum(frame_data)
				k_list.append(k_sum)
			avg_list = rollingavg(k_list,50)
		#	if "lg" in system:
		#		avg_list = avg_list/(40**2)
			#else:
			#	avg_list = avg_list/(12**2)
			if field == "zone":
				style = "solid"
			else: 
				style = "dashed"
			plt.plot(avg_list,linestyle=style)
			plt.ylim([0.0,.01])
		plt.savefig(system+".sum_over_K.pdf")

def rollingavg(list_in, window_size):
	rollingavg = []
	for i in range(window_size,len(list_in)-window_size):
		windowavg = np.sum(list_in[i-window_size:i+window_size])/(2*window_size)
		rollingavg.append(windowavg)
	rollingavg = np.array(rollingavg)
	return rollingavg

def plot_average_area_per_lipid(systems):
	
	avglist = []
	fig,ax = plt.subplots()


	for system in systems:
		
		data = np.loadtxt(system+'/'+system+'_boxarea.dat')
		areadata = data[1:]
		numlipidsperleaflet = data[0:1]
		apldata = areadata/(numlipidsperleaflet*1.0)
		
		average = np.mean(apldata)
		avglist.append(average)

	ax.bar(systems,avglist,label=systems)
	plt.show()































if __name__ == "__main__": 
	diff_mid_interface(["lgDT", "lgDY", "lgDG", "lgDO", "lgDP", "lgDL", "lgDX", "lgDB"], "5x29", "polar")
	#measure_H_epsilon_corr(["lgPO"], "empty")
	#measure_t0(["lgPO", "lgDG", "lgDY", "lgDT0", "lgDO", "lgDP", "lgDL", "lgDX", "lgDB"], "5x29")
	#diff_mid_interface(["lgPO"], "7k3g")
	#measure_H(["PO", "DG", "DY", "DT", "DL", "DO", "DP", "DX", "DB"], "5x29")
	#sum_over_K(["PO", "DG", "DY", "DT", "DL", "DO", "DP", "DX", "DB","lgPO", "lgDG", "lgDY", "lgDT"])
	#sum_over_K(["lgPO"])
	#sum_terms(["lgDG"], "5x29")
	#zoom_in(["lgDB","lgDL","lgDP","lgDX"])
	#zoom_in(["lgDT"])