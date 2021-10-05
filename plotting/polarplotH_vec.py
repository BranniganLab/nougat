import matplotlib.pyplot as plt
import numpy as np

protein = np.loadtxt("../5x29_coords_upr.dat",skiprows=1)
chain_up = []
for i in range(10):
	chain_up.append(protein[i])


protein = np.loadtxt("../5x29_coords_lwr.dat",skiprows=1)
chain_lw = []
for i in range(10):
  chain_lw.append(protein[i])

newfont = {'fontsize':'large',
 'fontweight':'bold',
 'color':'black',
 'verticalalignment': 'baseline'}


#read in from VMD traj
read_in_data = np.genfromtxt('../DL.lw.height.dat',missing_values='nan')

#figure out how many radial bins there are
counter = 1
flag = True
match_value = read_in_data[0,0]
while(flag==True):
  if read_in_data[counter,0] == match_value:
    flag = False
  else:
    counter = counter+1

N_r_bins = counter

#figure out how many azimuthal bins there are
Ntheta = len(read_in_data[0,:]) - 2

#figure out how many frames there are in the traj
Nframes = int(len(read_in_data[:,0])/N_r_bins)

#create a new array that has each frame in a different array level
newdata = np.zeros((N_r_bins, Ntheta, Nframes))
for x in range(Nframes):
  newdata[:,:,x] = read_in_data[x*N_r_bins:(x+1)*N_r_bins,2:]

avgHeight=np.nanmean(newdata, axis=2)

curvature = np.zeros((N_r_bins, Ntheta+2, Nframes))
curvature[:,1:31,:] = newdata
curvature[:,0,:] = curvature[:,30,:]
curvature[:,31,:] = curvature[:,1,:]





rad = read_in_data[0:N_r_bins,0]
rad = np.append(rad, read_in_data[N_r_bins-1,1])
the = np.linspace(0,2*np.pi,Ntheta+1)
theta,radius=np.meshgrid(the,rad)
avgHeight = avgHeight + 5.69





fig = plt.figure(figsize = (5,5))
ax = plt.subplot(projection="polar")
c = plt.pcolormesh(theta,radius,avgHeight,cmap="RdBu_r",zorder=0, vmax=0, vmin=-37)
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

