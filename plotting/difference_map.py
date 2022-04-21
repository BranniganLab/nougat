import matplotlib.pyplot as plt
import numpy as np
import warnings

leaflet_list = ["up", "lw"]

for leaflet in leaflet_list:
	DX = np.load("DX."+leaflet+".avgcurvature.npy")
	DT = np.load("DT."+leaflet+".avgcurvature.npy")
	diff = DX - DT
	#prepare to plot
	rad = np.linspace(0,72,18,endpoint=True)
	the = np.linspace(0,2*np.pi,30)
	radius,theta=np.meshgrid(rad, the, indexing='ij')
    #plotting section
	fig = plt.figure()
	ax = plt.subplot(projection="polar")
	c = plt.pcolormesh(theta,radius,diff,cmap="RdBu_r",shading='auto',vmax=.5,vmin=-.5)
	cbar = plt.colorbar(c)
	plt.axis('off')

	ax.set_xticklabels([])
	ax.set_yticklabels([])

	#figure = plt.gcf()
	fig.set_size_inches(6,6)
	plt.savefig(leaflet+"_difference.png", dpi = 700)
	plt.clf()
	plt.close()