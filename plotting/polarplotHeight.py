import matplotlib.pyplot as plt
import numpy as np

protein = np.loadtxt("../newanalysis/5x29_coords_upr.dat",skiprows=1)
chain_up = []
for i in range(10):
	chain_up.append(protein[i])


protein = np.loadtxt("../newanalysis/5x29_coords_lwr.dat",skiprows=1)
chain_lw = []
for i in range(10):
  chain_lw.append(protein[i])

newfont = {'fontsize':'large',
 'fontweight':'bold',
 'color':'black',
 'verticalalignment': 'baseline'}


Ntheta = 30
data = np.genfromtxt('DB.lw.height.dat',missing_values='np.nan')

rad = data[:,0]
rad = np.append(rad, [72.]) 
the = np.linspace(0,2*np.pi,Ntheta+1)
theta,radius=np.meshgrid(the,rad)
height = data[:,2:] + 28.28





fig = plt.figure(figsize = (5,5))
ax = plt.subplot(projection="polar")
c = plt.pcolormesh(theta,radius,height,cmap="RdBu_r",zorder=0, vmax=0, vmin=-30)
#c = plt.pcolormesh(theta,radius,height,cmap="RdBu_r",zorder=0)
cbar = plt.colorbar(c)
#cbar.ax.set_title('Enrichment',fontdict=newfont)
#cbar.ax.set_xlabel('Depletion',fontdict=newfont,labelpad=14.0)
cbar.ax.set_ylabel('Time-averaged Deformation (A)',fontdict=newfont,labelpad=20)
for i in range(0,10,2):
  plt.scatter(np.deg2rad(chain_lw[i+1]),chain_lw[i],c="black",linewidth=4,zorder=2)
plt.axis('off')

circle1 = plt.Circle((0,0),28.116, transform=ax.transData._b, color='black',linestyle='dashed',linewidth=4,fill=False)
#ax.add_artist(circle1)
ax.set_xticklabels([])
ax.set_yticklabels([])

plt.show()

