import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 


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
	"lgPO" : "darkorange",
	"lgDT" : "green",
	"lgDY" : "limegreen",
	"lgDG" : "orchid"
}

linetypedict = {
	"DT": "solid",
	"DL": "solid",
	"DX": "solid",
	"DB": "solid",
	"DY": "solid",
	"DO": "solid",
	"PO": "solid",
	"DP": "solid",
	"DG": "solid",
	"lgPO" : "dotted",
	"lgDT" : "dotted",
	"lgDY" : "dotted",
	"lgDG" : "dotted"
}


sys_list = ["DT", "DY", "DL", "DO", "DP", "PO", "DG", "DB", "DX", "lgDT", "lgDY", "lgDG", "lgPO"]
#sys_list = ["DT", "DY", "DL", "DO", "DP", "PO", "DG", "DB", "DX", "lgPO"]

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

def diff_mid_interface(systems, mol):
	for system in systems:
		filename_start = '/u1/home/js2746/Bending/PC/whole_mols/'+mol+'/'+system+'/'+system+'.'
		filename_end = '.C1A.C1B.cart.height.npy'
		fig = plt.figure()
		ax = plt.subplot()
		zzero = np.load(filename_start+'zzero'+filename_end)
		dims = np.shape(zzero)
		dim1 = np.linspace(0,dims[0],dims[0]+1)
		dim2 = np.linspace(0,dims[1],dims[1]+1)
		dim1vals,dim2vals=np.meshgrid(dim1, dim2, indexing='ij')
		zplus = np.load(filename_start+'zplus'+filename_end)
		diff = zzero-zplus
		avgdiff = np.nanmean(diff,axis=2)
		print(avgdiff)
		t0 = measure_t0(filename_start,filename_end)
		print(t0)
		avgdiff = avgdiff/t0
		c = plt.pcolormesh(dim1vals,dim2vals,avgdiff,cmap="RdBu_r",zorder=0,vmax=.5,vmin=-.5)
		cbar = plt.colorbar(c)
		plt.axis('off')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		#fig.set_size_inches(6,6)
		plt.savefig(system+"_zzerozplusdiff_cart.pdf", dpi = 700)
		plt.clf()
		plt.close()


def measure_t0(fname_start, fname_end):
	z1 = np.load(fname_start+'zone'+fname_end)
	z2 = np.load(fname_start+'ztwo'+fname_end)

	thickness = z1-z2

	avgthickness = np.nanmean(thickness,axis=2)

	leftcol = np.mean(avgthickness[:,0])
	rightcol =  np.mean(avgthickness[:,-1])
	toprow =  np.mean(avgthickness[0,:])
	botrow =   np.mean(avgthickness[-1,:])


	avgt0 = (leftcol+rightcol+toprow+botrow)/4.0

	avgt0 = avgt0/2.0

	return avgt0

def sum_over_K(systems):

	for system in systems:
		fig,ax = plt.subplots()
		for field in ["zone", "ztwo"]:
			data = np.load(system+'/newleaf_cart/'+system+'.'+field+'.C1A.C1B.cart.gausscurvature.npy')
			k_list = []
			for frame in range(np.shape(data)[2]):
				frame_data = data[:,:,frame]
				k_sum = np.nansum(frame_data)
				k_list.append(k_sum)
			avg_list = rollingavg(k_list,50)
			if field == "zone":
				style = "solid"
			else: 
				style = "dashed"
			plt.plot(avg_list,linestyle=style)
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
	diff_mid_interface(["lgPO", "lgDG", "lgDY", "lgDT"], "5x29")
	#diff_mid_interface(["lgPO"], "7k3g")