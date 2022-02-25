import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings

name_list = ["DB"]
#name_list = ["DL", "DT", "DG", "DX", "PO", "DY", "DB", "DO", "DP"]
field_list = ["zone","ztwo","zplus","zzero"]
#field_list = ["test"]


def dimensions_analyzer(data):
  #figure out how many radial bins there are
  counter = 1
  flag = True
  match_value = data[0,0]
  while(flag==True):
    try:
      if data[counter,0] == match_value:
        flag = False
      else:
        counter = counter+1
    except IndexError:
      flag = False
  N_r_bins = counter
  dr = data[0,1] - data[0,0]

  #figure out how many azimuthal bins there are
  N_theta_bins = len(data[0,:]) - 2
  dtheta = (np.pi*2)/N_theta_bins

  #figure out how many frames there are in the traj
  Nframes = int(len(data[:,0])/N_r_bins)

  return N_r_bins, dr, N_theta_bins, dtheta, Nframes


def empty_neighbor_test(data, Nframes, N_r_bins, N_theta_bins):
  nan_test = np.array(data, copy=True)
  nan_test2 = np.array(data, copy=True)
  knan_test = np.array(data, copy=True)
  for frm in range(Nframes):
    for row in range(1,N_r_bins-1):
      for col in range(1,N_theta_bins+1):
        if nan_test2[row-1,col,frm] == True:
          nan_test[row,col,frm] = True
          knan_test[row,col,frm] = True
        elif nan_test2[row+1,col,frm] == True:
          nan_test[row,col,frm] = True
          knan_test[row,col,frm] = True
        elif nan_test2[row,col-1,frm] == True:
          nan_test[row,col,frm] = True
          knan_test[row,col,frm] = True
        elif nan_test2[row,col+1,frm] == True:
          nan_test[row,col,frm] = True
          knan_test[row,col,frm] = True
        elif nan_test2[row+1,col+1,frm] == True:
          knan_test[row,col,frm] = True 
        elif nan_test2[row-1,col+1,frm] == True:
          knan_test[row,col,frm] = True
        elif nan_test2[row+1,col-1,frm] == True:
          knan_test[row,col,frm] = True 
        elif nan_test2[row-1,col-1,frm] == True:
          knan_test[row,col,frm] = True 

  nan_test[0,:,:] = True
  nan_test[N_r_bins-1,:,:] = True
  nan_test[:,0,:] = True
  nan_test[:,N_theta_bins+1,:] = True

  knan_test[0,:,:] = True
  knan_test[N_r_bins-1,:,:] = True
  knan_test[:,0,:] = True
  knan_test[:,N_theta_bins+1,:] = True

  return nan_test, knan_test


def plot_maker(radius, theta, data, name, field, Vmax, Vmin, protein, dataname):
  fig = plt.figure()
  ax = plt.subplot(projection="polar")
  c = plt.pcolormesh(theta,radius,data,cmap="RdBu_r",zorder=0,vmax=Vmax,vmin=Vmin)
  cbar = plt.colorbar(c)
  for i in range(0,10,2):
    plt.scatter(np.deg2rad(protein[i+1]),protein[i],c="black",linewidth=4,zorder=2)
  plt.axis('off')

  circle1 = plt.Circle((0,0),28.116, transform=ax.transData._b, color='black',linestyle='dashed',linewidth=4,fill=False)
  if field == "zone":
    ax.add_artist(circle1)

  ax.set_xticklabels([])
  ax.set_yticklabels([])

  fig.set_size_inches(6,6)
  plt.savefig(name+"_"+field+"_"+dataname+".png", dpi = 700)
  plt.clf()
  plt.close()


def measure_curvature(Nframes, N_r_bins, N_theta_bins, knan_test, nan_test, curvature_inputs, curvature_outputs, kgauss_outputs, dr, dtheta):
  for frm in range(Nframes):
    for row in range(N_r_bins):
      for col in range(N_theta_bins+2):
        if knan_test[row,col,frm] == False:

          #calculate d2h/dr2
          del2r = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
          del2r = del2r / dr**2

          #calculate dh/dr
          delr = (curvature_inputs[row+1,col,frm] - curvature_inputs[row-1,col,frm])/(2*dr)
          
          #calculate d2h/drdtheta
          delrdeltheta = (curvature_inputs[row+1,col+1,frm] - curvature_inputs[row+1,col-1,frm] - curvature_inputs[row-1,col+1,frm] + curvature_inputs[row-1,col-1,frm])
          delrdeltheta = delrdeltheta/(4*dr*dtheta)

          #calculate dh/dtheta
          deltheta = (curvature_inputs[row,col+1,frm] - curvature_inputs[row,col-1,frm])/(2*dtheta)

          #calculate d2h/dtheta2
          del2theta = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
          del2theta = del2theta / dtheta**2

          #calculate coefficients
          r = (row*dr) + (dr/2)
          c1 = 1 / r
          c2 = 1 / r**2
          c3 = 1 / r**3
          c4 = 1 / r**4

          #calculate polar laplacian and gaussian curvature
          curvature_outputs[row,col,frm] = del2r + c1*delr + c2*del2theta
          kgauss_outputs[row,col,frm] = (-1*c1*del2r*delr) + (c2*(delrdeltheta**2 - (del2r*del2theta))) + (-2*c3*delrdeltheta*deltheta) + (c4*deltheta**2)

        elif nan_test[row,col,frm] == False:

          #calculate d2h/dr2
          del2r = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
          del2r = del2r / dr**2

          #calculate dh/dr
          delr = (curvature_inputs[row+1,col,frm] - curvature_inputs[row-1,col,frm])/(2*dr)
        
          #calculate d2h/dtheta2
          del2theta = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
          del2theta = del2theta / dtheta**2

          #calculate coefficients
          r = (row*dr) + (dr/2)
          c1 = 1 / r
          c2 = 1 / r**2

          curvature_outputs[row,col,frm] = del2r + c1*delr + c2*del2theta
          kgauss_outputs[row,col,frm] = np.nan

        else:
          curvature_outputs[row,col,frm] = np.nan
          kgauss_outputs[row,col,frm] = np.nan 

  return curvature_outputs, kgauss_outputs




#---------------------------------------------------------------------#

def main():
  for name in name_list:
    for field in field_list:

      #read in protein helix coordinates
      protein_coords = np.loadtxt(name+"_helcoords_"+field+".dat",skiprows=1)
      protein = []
      for i in range(10):
      	protein.append(protein_coords[i])

      #read in heights from VMD traj
      height_data = np.genfromtxt(name+'.'+field+'.height.dat',missing_values='nan',filling_values=np.nan)
      density_data = np.genfromtxt(name+'.'+field+'.density.dat')

      #strip r values from density info
      density = density_data[:,2:]

      #get bin info
      N_r_bins, dr, N_theta_bins, dtheta, Nframes = dimensions_analyzer(height_data)

      #create a new array that has each frame in a different array level
      height = np.zeros((N_r_bins, N_theta_bins, Nframes))
      for x in range(Nframes):
        height[:,:,x] = height_data[x*N_r_bins:(x+1)*N_r_bins,2:]

      #create arrays for storing curvature data
      curvature_inputs = np.zeros((N_r_bins, N_theta_bins+2, Nframes))
      curvature_outputs = np.zeros((N_r_bins, N_theta_bins+2, Nframes))
      kgauss_outputs = np.zeros((N_r_bins, N_theta_bins+2, Nframes))

      #wrap the inputs in the theta direction for calculating curvature
      curvature_inputs[:,1:31,:] = height
      curvature_inputs[:,0,:] = curvature_inputs[:,30,:]
      curvature_inputs[:,31,:] = curvature_inputs[:,1,:]

      #if a bin is empty, you can't measure its curvature
      nan_test = np.isnan(curvature_inputs)

      #if a bin is empty, you can't (nicely) measure the curvature of its neighbors
      nan_test, knan_test = empty_neighbor_test(nan_test, Nframes, N_r_bins, N_theta_bins)

      #prep plot dimensions
      rad = height_data[0:N_r_bins,0]
      rad = np.append(rad, height_data[N_r_bins-1,1])
      the = np.linspace(0,2*np.pi,N_theta_bins+1)
      radius,theta=np.meshgrid(rad, the, indexing='ij')

      #produce average height (dtype == 0) and curvature (dtype == 1) plots
      for dtype in range(3):
        if dtype == 0:

          #if a bin only has lipids in it <10% of the time, it shouldn't be considered part of the membrane
          for row in range(N_r_bins):
            for col in range(N_theta_bins):
              zerocount = np.count_nonzero(height[row,col,:])
              count = np.count_nonzero(np.isnan(height[row,col,:]))
              if (zerocount-count)/Nframes <= .1:
                height[row,col,:] = np.nan

          #take the average height over all frames
          with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            avgHeight=np.nanmean(height, axis=2)

          #save as file for debugging / analysis
          np.save(name+'.'+field+'.avgheight.npy', avgHeight)

          #plot and save
          plot_maker(radius, theta, avgHeight, name, field, 0, -45, protein, "avgHeight")
          print(name+" "+field+" height done!")

        elif dtype == 1:

          #measure the laplacian and gaussian curvatures
          curvature_outputs, kgauss_outputs = measure_curvature(Nframes, N_r_bins, N_theta_bins, knan_test, nan_test, curvature_inputs, curvature_outputs, kgauss_outputs, dr, dtheta)

          #unwrap along theta direction
          meancurvature = curvature_outputs[:,1:N_theta_bins+1,:]
          kcurvature = kgauss_outputs[:,1:N_theta_bins+1,:]

          #take the average curvatures over all frames
          with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            avgcurvature=np.nanmean(meancurvature, axis=2)
            avgkcurvature=np.nanmean(kcurvature, axis=2)

          #save as file for debugging / analysis
          np.save(name+'.'+field+'.avgcurvature.npy',avgcurvature)
          np.save(name+'.'+field+'.avgKcurvature.npy',avgkcurvature)

          #laplacian plotting section
          plot_maker(radius, theta, avgcurvature, name, field, .1, -.1, protein, "curvature")
          print(name+" "+field+" laplacian done!")

          #gaussian plotting section
          plot_maker(radius, theta, avgkcurvature, name, field, .1, -.1, protein, "gausscurvature")
          print(name+" "+field+" gaussian curvature done!")

        elif dtype == 2:

          #save as file for debuggging / analysis
          np.save(name+'.'+field+'.avgdensity.npy',density)

          #plot and save
          plot_maker(radius, theta, density, name, field, 0, 2, protein, "density")
          print(name+" "+field+" density done!")



if __name__ == "__main__": 
  main()