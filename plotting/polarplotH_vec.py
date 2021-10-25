import matplotlib.pyplot as plt
import numpy as np
import warnings

newfont = {'fontsize':'large',
 'fontweight':'bold',
 'color':'black',
 'verticalalignment': 'baseline'}

#name_list = ["DL"]
name_list = ["DL", "DT", "DG", "DX", "PO", "DY", "DB"]
leaflet_list = ["lw", "up"]



for name in name_list:
  for leaflet in leaflet_list:

    #read in protein helix coordinates
    protein = np.loadtxt(name+"_helcoords_"+leaflet+"r.dat",skiprows=1)
    chain = []
    for i in range(10):
    	chain.append(protein[i])

    #read in heights from VMD traj
    read_in_data = np.genfromtxt(name+'.'+leaflet+'.height.dat',missing_values='nan')

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
    dr = read_in_data[0,1] - read_in_data[0,0]

    #figure out how many azimuthal bins there are
    Ntheta = len(read_in_data[0,:]) - 2
    dtheta = (np.pi*2)/Ntheta

    #figure out how many frames there are in the traj
    Nframes = int(len(read_in_data[:,0])/N_r_bins)

    #create a new array that has each frame in a different array level
    newdata = np.zeros((N_r_bins, Ntheta, Nframes))
    for x in range(Nframes):
      newdata[:,:,x] = read_in_data[x*N_r_bins:(x+1)*N_r_bins,2:]


    for dtype in range(2):
      if dtype == 0:
        #take the average height over all frames
        with warnings.catch_warnings():
          warnings.simplefilter("ignore", category=RuntimeWarning)
          avgHeight=np.nanmean(newdata, axis=2)

        #save as file
        np.save(name+'.'+leaflet+'.avgheight.npy', avgHeight)

        #prepare to plot
        rad = read_in_data[0:N_r_bins,0]
        rad = np.append(rad, read_in_data[N_r_bins-1,1])
        the = np.linspace(0,2*np.pi,Ntheta+1)
        radius,theta=np.meshgrid(rad, the, indexing='ij')

        #plotting section
        fig = plt.figure()
        ax = plt.subplot(projection="polar")
        c = plt.pcolormesh(theta,radius,avgHeight,cmap="RdBu_r",zorder=0,vmax=0,vmin=-45)
        cbar = plt.colorbar(c)
        for i in range(0,10,2):
          plt.scatter(np.deg2rad(chain[i+1]),chain[i],c="black",linewidth=4,zorder=2)
        plt.axis('off')

        circle1 = plt.Circle((0,0),28.116, transform=ax.transData._b, color='black',linestyle='dashed',linewidth=4,fill=False)
        if leaflet == "up":
          ax.add_artist(circle1)

        ax.set_xticklabels([])
        ax.set_yticklabels([])

        #figure = plt.gcf()
        fig.set_size_inches(6,6)
        plt.savefig(name+"_"+leaflet+"_avgHeight.png", dpi = 700)
        plt.clf()
        plt.close()

        print(name+" height done!")

      elif dtype == 1:
        #create arrays for storing curvature data
        curvature_inputs = np.zeros((N_r_bins, Ntheta+2, Nframes))
        curvature_outputs = np.zeros((N_r_bins, Ntheta+2, Nframes))

        #wrap the inputs in the theta direction for calculating curvature
        curvature_inputs[:,1:31,:] = newdata
        curvature_inputs[:,0,:] = curvature_inputs[:,30,:]
        curvature_inputs[:,31,:] = curvature_inputs[:,1,:]

        #if a bin is empty, you can't measure its curvature
        nan_test = np.isnan(curvature_inputs)

        #if a bin is empty, you can't (nicely) measure the curvature of its neighbors
        nan_test2 = np.array(nan_test, copy=True)
        for frm in range(Nframes):
          for row in range(1,N_r_bins-1):
            for col in range(1,Ntheta+1):
              if nan_test2[row-1,col,frm] == True:
                nan_test[row,col,frm] = True
              elif nan_test2[row+1,col,frm] == True:
                nan_test[row,col,frm] = True
              elif nan_test2[row,col-1,frm] == True:
                nan_test[row,col,frm] = True
              elif nan_test2[row,col+1,frm] == True:
                nan_test[row,col,frm] = True

        nan_test[0,:,:] = True
        nan_test[N_r_bins-1,:,:] = True
        nan_test[:,0,:] = True
        nan_test[:,Ntheta+1,:] = True

        #polar laplacian = d2h/dr2 + 1/r dh/dr + 1/r^2 d2h/dtheta2

        for frm in range(Nframes):
          for row in range(N_r_bins):
            for col in range(Ntheta+2):
              if nan_test[row,col,frm] == False:
                
                #calculate d2h/dr2
                del2r = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
                del2r = del2r / dr**2

                #calculate dh/dr
                delr = curvature_inputs[row-1,col,frm] - curvature_inputs[row+1,col,frm]/(2*dr)
                
                #calculate d2h/dtheta2
                del2theta = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
                del2theta = del2theta / dr**2

                #calculate coefficients
                rad = (row*dr) + (dr/2)
                c1 = 1 / rad
                c2 = 1 / rad**2

                #calculate polar laplacian
                curvature_outputs[row,col,frm] = del2r + c1*delr + c2*del2theta

              else:
                curvature_outputs[row,col,frm] = np.nan

        curvature = curvature_outputs[:,1:Ntheta+1,:]
        with warnings.catch_warnings():
          warnings.simplefilter("ignore", category=RuntimeWarning)
          avgcurvature=np.nanmean(curvature, axis=2)

        np.save(name+'.'+leaflet+'.avgcurvature.npy',avgcurvature)

        #prepare to plot
        rad = read_in_data[0:N_r_bins,0]
        rad = np.append(rad, read_in_data[N_r_bins-1,1])
        the = np.linspace(0,2*np.pi,Ntheta+1)
        radius,theta=np.meshgrid(rad, the, indexing='ij')

        #plotting section
        fig = plt.figure()
        ax = plt.subplot(projection="polar")
        c = plt.pcolormesh(theta,radius,avgcurvature,cmap="RdBu_r",zorder=0,vmax=.1,vmin=-.1)
        #c = plt.pcolormesh(theta,radius,height,cmap="RdBu_r",zorder=0)
        cbar = plt.colorbar(c)
        for i in range(0,10,2):
          plt.scatter(np.deg2rad(chain[i+1]),chain[i],c="black",linewidth=4,zorder=2)
        plt.axis('off')

        circle1 = plt.Circle((0,0),28.116, transform=ax.transData._b, color='black',linestyle='dashed',linewidth=4,fill=False)
        if leaflet == "up":
          ax.add_artist(circle1)

        ax.set_xticklabels([])
        ax.set_yticklabels([])

        #figure = plt.gcf()
        fig.set_size_inches(6,6)
        plt.savefig(name+"_"+leaflet+"_curvature.png", dpi = 700)
        plt.clf()
        plt.close()
        print(name+" curvature done!")