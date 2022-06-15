import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings

readbeads = 0
protein_onoff = 0
polar = True
name_list = ["PO"]
bead_dict = {
  "DT" : ['C2A.C2B'],
  "DL" : ['C2A.C2B', 'C3A.C3B'],
  "DY" : ['D2A.D2B', 'C3A.C3B'],
  "DO" : ['D2A.D2B', 'C3A.C3B', 'C4A.C4B'],
  "PO" : ['D2A.C2B', 'C3A.C3B', 'C4A.C4B'],
  "lgPO" : ['D2A.C2B', 'C3A.C3B', 'C4A.C4B'],
  "DP" : ['C2A.C2B', 'C3A.C3B', 'C4A.C4B'],
  "DB" : ['C2A.C2B', 'C3A.C3B', 'C4A.C4B', 'C5A.C5B'],
  "DG" : ['C2A.C2B', 'D3A.D3B', 'C4A.C4B', 'C5A.C5B'],
  "DX" : ['C2A.C2B', 'C3A.C3B', 'C4A.C4B', 'C5A.C5B', 'C6A.C6B']
}

# These determine the scale in your png files
# adjust as needed
height_min = -60
height_max = 60
mean_curv_min = 0.05
mean_curv_max = -0.05
gauss_curv_min = -0.005
gauss_curv_max = 0.005
density_min = 0
density_max = 2

field_list = ["zone","ztwo","zplus","zzero"]

def dimensions_analyzer(data, polar):
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
  N1_bins = counter
  d1 = data[0,1] - data[0,0]

  #figure out how many azimuthal bins there are
  N2_bins = len(data[0,:]) - 2
  if polar is True:
    d2 = (np.pi*2)/N2_bins
  elif polar is False:
    d2 = d1

  #figure out how many frames there are in the traj
  Nframes = int(len(data[:,0])/N1_bins)

  return N1_bins, d1, N2_bins, d2, Nframes


def empty_neighbor_test(data, Nframes, N1_bins, N2_bins):
  nan_test = np.array(data, copy=True)
  nan_test2 = np.array(data, copy=True)
  knan_test = np.array(data, copy=True)
  for frm in range(Nframes):
    for row in range(1,N1_bins-1):
      for col in range(1,N2_bins+1):
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
  nan_test[N1_bins-1,:,:] = True
  nan_test[:,0,:] = True
  nan_test[:,N2_bins+1,:] = True

  knan_test[0,:,:] = True
  knan_test[N1_bins-1,:,:] = True
  knan_test[:,0,:] = True
  knan_test[:,N2_bins+1,:] = True

  return nan_test, knan_test


def plot_maker(dim1vals, dim2vals, data, name, field, Vmax, Vmin, protein, dataname, bead, polar):
  fig = plt.figure()

  if polar is True:
    ax = plt.subplot(projection="polar")
  elif polar is False:
    ax = plt.subplot()
    
  c = plt.pcolormesh(dim1vals,dim2vals,data,cmap="RdBu_r",zorder=0,vmax=Vmax,vmin=Vmin)
  cbar = plt.colorbar(c)

  if protein != False:
    print(protein)
    for i in range(0,10,2):
      protein[i+1] = np.deg2rad(protein[i+1])
      if polar is False:
        protein[i], protein[i+1] = convert_to_cart(protein[i],protein[i+1])
      plt.scatter(protein[i+1],protein[i],c="black",linewidth=4,zorder=2)
    circle1 = plt.Circle((0,0),28.116, transform=ax.transData._b, color='black',linestyle='dashed',linewidth=4,fill=False)
    if field == "zone":
      ax.add_artist(circle1)

  plt.axis('off')
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  fig.set_size_inches(6,6)

  if polar is True:
    coordsys = "polar"
  elif polar is False:
    coordsys = "cart"

  if bead is False:
    plt.savefig(name+"_"+field+"_"+coordsys+"_"+dataname+".png", dpi = 700)
  else:
    plt.savefig(name+"_"+bead+"_"+field+"_"+coordsys+"_"+dataname+".png", dpi = 700)
  plt.clf()
  plt.close()

def convert_to_cart(rval, thetaval):
  xval = rval*np.cos(thetaval)
  yval = rval*np.sin(thetaval)
  return xval, yval

def measure_curvature_polar(Nframes, N1_bins, N2_bins, knan_test, nan_test, curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, d1, d2):
  for frm in range(Nframes):
    for row in range(N1_bins):
      for col in range(N2_bins+2):
        if knan_test[row,col,frm] == False:

          #calculate d2h/dr2
          del2r = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
          del2r = del2r / d1**2

          #calculate dh/dr
          delr = (curvature_inputs[row+1,col,frm] - curvature_inputs[row-1,col,frm])/(2*d1)
          
          #calculate d2h/drdtheta
          delrdeltheta = (curvature_inputs[row+1,col+1,frm] - curvature_inputs[row+1,col-1,frm] - curvature_inputs[row-1,col+1,frm] + curvature_inputs[row-1,col-1,frm])
          delrdeltheta = delrdeltheta/(4*d1*d2)

          #calculate dh/dtheta
          deltheta = (curvature_inputs[row,col+1,frm] - curvature_inputs[row,col-1,frm])/(2*d2)

          #calculate d2h/dtheta2
          del2theta = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
          del2theta = del2theta / d2**2

          #calculate coefficients
          r = (row*d1) + (d1/2)
          c1 = 1 / r
          c2 = 1 / r**2
          c3 = 1 / r**3
          c4 = 1 / r**4
          theta = ((col-1)*d2) + (d2/2)     ;# col-1 because this is wrapped in theta direction

          #calculate normal vector x,y components
          normalization_factor = np.sqrt(1 + c2*deltheta**2 + delr**2)
          norm_vec_x = (c1*np.sin(theta)*deltheta) - (np.cos(theta)*delr) / normalization_factor
          norm_vec_y = (-1*c1*np.cos(theta)*deltheta) - (np.sin(theta)*delr) / normalization_factor
          norm_vec_z = 1 / normalization_factor

          #calculate polar laplacian and gaussian curvature
          curvature_outputs[row,col,frm] = del2r + c1*delr ;# + c2*del2theta
          kgauss_outputs[row,col,frm] = (-1*c1*del2r*delr) + (c2*(delrdeltheta**2 - (del2r*del2theta))) + (-2*c3*delrdeltheta*deltheta) + (c4*deltheta**2)
          normal_vector_outputs[row,col*3,frm] = norm_vec_x
          normal_vector_outputs[row,col*3+1,frm] = norm_vec_y
          normal_vector_outputs[row,col*3+2,frm] = norm_vec_z

        elif nan_test[row,col,frm] == False:

          #calculate d2h/dr2
          del2r = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
          del2r = del2r / d1**2

          #calculate dh/dr
          delr = (curvature_inputs[row+1,col,frm] - curvature_inputs[row-1,col,frm])/(2*d1)

          #calculate dh/dtheta
          deltheta = (curvature_inputs[row,col+1,frm] - curvature_inputs[row,col-1,frm])/(2*d2)
        
          #calculate d2h/dtheta2
          del2theta = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
          del2theta = del2theta / d2**2

          #calculate coefficients
          r = (row*d1) + (d1/2)
          c1 = 1 / r
          c2 = 1 / r**2
          theta = ((col-1)*d2) + (d2/2)   ;# col-1 because this is wrapped in theta direction

          #calculate normal vector x,y components
          normalization_factor = np.sqrt(1 + c2*deltheta**2 + delr**2)
          norm_vec_x = (c1*np.sin(theta)*deltheta) - (np.cos(theta)*delr) / normalization_factor
          norm_vec_y = (-1*c1*np.cos(theta)*deltheta) - (np.sin(theta)*delr) / normalization_factor
          norm_vec_z = 1 / normalization_factor

          curvature_outputs[row,col,frm] = del2r + c1*delr ;# + c2*del2theta
          kgauss_outputs[row,col,frm] = np.nan
          normal_vector_outputs[row,col*3,frm] = norm_vec_x
          normal_vector_outputs[row,col*3+1,frm] = norm_vec_y
          normal_vector_outputs[row,col*3+2,frm] = norm_vec_z

        else:
          curvature_outputs[row,col,frm] = np.nan
          kgauss_outputs[row,col,frm] = np.nan 
          normal_vector_outputs[row,col*3,frm] = np.nan
          normal_vector_outputs[row,col*3+1,frm] = np.nan
          normal_vector_outputs[row,col*3+2,frm] = np.nan

  return curvature_outputs, kgauss_outputs, normal_vector_outputs


def coord_format(value):
  rounded = round(value,3)
  leftside,rightside = str(rounded).split('.')
  if len(rightside) < 3:
    rightside = rightside+(' '*(3-len(rightside)))
  if len(leftside) < 4:
    leftside = (' '*(4-len(leftside)))+leftside
  final_value = leftside+'.'+rightside
  return final_value

def bin_format(value):
  strval = str(value)
  length = len(strval)
  final_value = (' '*(3-length))+strval+'.00'
  return final_value

def Make_surface_PDB(data,name,field,d1,d2,f,serial,bead,polar):
  resseqnum = 1
  atom_name = 'SURF'
  chain = 'X'
  row,col = data.shape
  beadnum = str(bead[1])
  beadname = "C"+beadnum+"  "

  for rbin in range(row):
    for thetabin in range(col):
      if str(data[rbin][thetabin]) != "nan":
        seriallen = 5-(len(str(serial)))
        resseqlen = 4-(len(str(resseqnum)))
        x = (d1*rbin + .5*d1)*(np.cos(thetabin*d2 + 0.5*d2))
        x = coord_format(x)
        y = (d1*rbin + .5*d1)*(np.sin(thetabin*d2 + 0.5*d2))
        y = coord_format(y)
        z = coord_format(data[rbin][thetabin])
        rnum = bin_format(rbin)
        thetanum = bin_format(thetabin)
        print('HETATM'+(' '*seriallen)+str(serial)+' '+atom_name+' '+beadname+chain+(' '*resseqlen)+str(resseqnum)+'    '+x+y+z+rnum+thetanum+'      '+field[:4]+'  ',file=f)
        serial += 1
        resseqnum +=1
  return serial

def normalize_vectors_in_array(array, Nr, Ntheta):
  for row in range(Nr):
    for col in range(Ntheta):
      vector = array[row,col*3:((col+1)*3)]
      if vector[0] != "nan":
        norm_vec = 4*vector / np.sqrt(np.sum(vector**2))
      else:
        norm_vec = vector 
      array[row,col*3:((col+1)*3)] = norm_vec

  return array 

def gen_avg_tilt(name, field, polar):
  vector_data = np.genfromtxt(name+'.'+field+'.tilt.dat', missing_values='nan',filling_values=np.nan)
  height_data = np.genfromtxt(name+'.'+field+'.height.dat',missing_values='nan',filling_values=np.nan)
  normal_data = np.load(name+'.'+field+'.normal_vectors.npy')

  #match direction between lipid vector and normal vector
  if field == "zone":
    normal_data = normal_data * -1

  #strip r values from tilt info
  vector_data = vector_data[:,2:]

  #get bin info
  N1_bins, d1, N2_bins, d2, Nframes = dimensions_analyzer(height_data, polar)

  vector = np.zeros((N1_bins, N2_bins*3, Nframes))
  for x in range(Nframes):
    vector[:,:,x] = vector_data[x*N1_bins:(x+1)*N1_bins,:]

  #if lipids are not in the bin at least 10% of the time, exclude that bin
  for row in range(N1_bins):
    for col in range(N2_bins):
      zerocount = np.count_nonzero(vector[row,col,:])
      count = np.count_nonzero(np.isnan(vector[row,col,:]))
      if (zerocount-count)/Nframes <= .1:
        vector[row,col,:] = np.nan

  #tilt is the difference between lipid tail vectors and surface normal vectors
  tilt = vector - normal_data

  #take average across frames
  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    avgtilt=np.nanmean(tilt, axis=2)

  np.savetxt(name+'.'+field+'.avgtilt.dat',avgtilt,delimiter = ',',fmt='%10.7f')  

#---------------------------------------------------------------------#

def output_analysis(name, field, protein, data_opt, bead, surffile, serial, polar):
  
  if polar is True:
    coordsys = "polar"
  elif polar is False:
    coordsys = "cart"

  #read in heights from VMD traj
  if bead is False:
    height_data = np.genfromtxt(name+'.'+field+'.'+coordsys+'.height.dat',missing_values='nan',filling_values=np.nan)
    density_data = np.genfromtxt(name+'.'+field+'.'+coordsys+'.density.dat')
  else:
    height_data = np.genfromtxt(name+'.'+bead+'.'+field+'.'+coordsys+'.height.dat',missing_values='nan',filling_values=np.nan)
    density_data = np.genfromtxt(name+'.'+bead+'.'+field+'.'+coordsys+'.density.dat')

  #strip dim1 values from density info
  density = density_data[:,2:]

  #get bin info
  N1_bins, d1, N2_bins, d2, Nframes = dimensions_analyzer(height_data, polar)

  #create a new array that has each frame in a different array level
  height = np.zeros((N1_bins, N2_bins, Nframes))
  for x in range(Nframes):
    height[:,:,x] = height_data[x*N1_bins:(x+1)*N1_bins,2:]

  #create arrays for storing curvature data
  if polar is True:
    curvature_inputs = np.zeros((N1_bins, N2_bins+2, Nframes))
    curvature_outputs = np.zeros((N1_bins, N2_bins+2, Nframes))
    kgauss_outputs = np.zeros((N1_bins, N2_bins+2, Nframes))
    normal_vector_outputs = np.zeros((N1_bins, 3*(N2_bins+2), Nframes))
  elif polar is False:
    curvature_inputs = np.zeros((N1_bins+2, N2_bins+2, Nframes))
    curvature_outputs = np.zeros((N1_bins+2, N2_bins+2, Nframes))
    kgauss_outputs = np.zeros((N1_bins+2, N2_bins+2, Nframes))
    normal_vector_outputs = np.zeros((N1_bins+2, 3*(N2_bins+2), Nframes))

  #wrap the inputs in the theta direction for calculating curvature
  curvature_inputs[:,1:(N2_bins+1),:] = height
  curvature_inputs[:,0,:] = curvature_inputs[:,N2_bins,:]
  curvature_inputs[:,(N2_bins+1),:] = curvature_inputs[:,1,:]

  #if cartesian, wrap in both directions
  if polar == 0:
    curvature_inputs[1:(N1_bins+1),:,:] = height
    curvature_inputs[0,:,:] = curvature_inputs[N1_bins,:,:]
    curvature_inputs[(N1_bins+1),:,:] = curvature_inputs[1,:,:]

  #prep plot dimensions
  dim1 = height_data[0:N1_bins,0]
  dim1 = np.append(dim1, height_data[N1_bins-1,1])
  if polar is True:
    dim2 = np.linspace(0,2*np.pi,N2_bins+1)
  elif polar is False:
    dim2 = np.linspace(0,N2_bins+1,N2_bins+1)
  dim1vals,dim2vals=np.meshgrid(dim1, dim2, indexing='ij')

  #produce average height (dtype == 0) and curvature (dtype == 1) plots
  for dtype in range(data_opt):
    if dtype == 0:

      #if a bin only has lipids in it <10% of the time, it shouldn't be considered part of the membrane
      for row in range(N1_bins):
        for col in range(N2_bins):
          zerocount = np.count_nonzero(height[row,col,:])
          count = np.count_nonzero(np.isnan(height[row,col,:]))
          if (zerocount-count)/Nframes <= .1:
            height[row,col,:] = np.nan

      #take the average height over all frames
      with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avgHeight=np.nanmean(height, axis=2)

      #make plots!
      plot_maker(dim1vals, dim2vals, avgHeight, name, field, height_max, height_min, protein, "avgHeight", bead, polar)

      #save as file for debugging / analysis AND make PDB!
      if bead is False: 
        np.savetxt(name+'.'+field+'.'+coordsys+'.avgheight.dat', avgHeight,delimiter = ',',fmt='%10.5f')
        serial = Make_surface_PDB(avgHeight,name,field,d1,d2,surffile,serial,'C1  ',polar)
        print(name+" "+field+" height done!")
      else:
        np.savetxt(name+'.'+bead+'.'+field+'.'+coordsys+'.avgheight.dat', avgHeight,delimiter = ',',fmt='%10.5f')
        serial = Make_surface_PDB(avgHeight,name,field,d1,d2,surffile,serial,bead,polar)
        print(name+' '+bead+' '+field+" height done!")

    elif dtype == 1:
      #if a bin is empty, you can't measure its curvature
      nan_test = np.isnan(curvature_inputs)

      #if a bin is empty, you can't (nicely) measure the curvature of its neighbors
      nan_test, knan_test = empty_neighbor_test(nan_test, Nframes, N1_bins, N2_bins)

      #measure the laplacian and gaussian curvatures
      if polar is True:
        curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_polar(Nframes, N1_bins, N2_bins, knan_test, nan_test, curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, d1, d2)
      elif polar is False:
        curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_cart()

      #unwrap along dim2 direction
      meancurvature = curvature_outputs[:,1:N2_bins+1,:]
      kcurvature = kgauss_outputs[:,1:N2_bins+1,:]
      normal_vectors = normal_vector_outputs[:,3:3*(N2_bins+1),:]

      #if cartesian, unwrap along dim1 direction too
      if polar is False:
        meancurvature = meancurvature[1:N1_bins+1,:,:]
        kcurvature = kcurvature[1:N1_bins+1,:,:]
        normal_vectors = normal_vectors[1:N1_bins+1,:,:]

      #take the average curvatures over all frames
      with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avgcurvature=np.nanmean(meancurvature, axis=2)
        avgkcurvature=np.nanmean(kcurvature, axis=2)

      #make plots!
      plot_maker(dim1vals, dim2vals, avgkcurvature, name, field, gauss_curv_max, gauss_curv_min, protein, "gausscurvature", bead, polar)
      plot_maker(dim1vals, dim2vals, avgcurvature, name, field, mean_curv_max, mean_curv_min, protein, "curvature", bead, polar)

      #save as files for debugging / analysis
      if bead is False: 
        #save as file for debugging / analysis
        np.savetxt(name+'.'+field+'.'+coordsys+'.avgcurvature.dat',avgcurvature,delimiter = ',',fmt='%10.7f')
        np.savetxt(name+'.'+field+'.'+coordsys+'.avgKcurvature.dat',avgkcurvature,delimiter = ',',fmt='%10.7f')
        np.save(name+'.'+field+'.'+coordsys+'.normal_vectors.npy',normal_vectors)
        print(name+" "+field+" curvatures done!")
      else:
        np.savetxt(name+'.'+bead+'.'+field+'.'+coordsys+'.avgcurvature.dat',avgcurvature,delimiter = ',',fmt='%10.7f')
        np.savetxt(name+'.'+bead+'.'+field+'.'+coordsys+'.avgKcurvature.dat',avgkcurvature,delimiter = ',',fmt='%10.7f')
        np.save(name+'.'+bead+'.'+field+'.'+coordsys+'.normal_vectors.npy',normal_vectors)
        print(name+' '+bead+' '+field+" curvatures done!")

    elif dtype == 2:

      #make plot!
      plot_maker(dim1vals, dim2vals, density, name, field, density_max, density_min, protein, "density", bead, polar)

      #save as file for debugging / analysis
      if bead is False:
        np.savetxt(name+'.'+field+'.avgdensity.dat',density,delimiter = ',',fmt='%10.7f')
        print(name+" "+field+" density done!")
      else:
        np.savetxt(name+'.'+bead+'.'+field+'.avgdensity.dat',density,delimiter = ',',fmt='%10.7f')
        print(name+' '+bead+' '+field+" density done!")
  return serial



if __name__ == "__main__": 
  for name in name_list:
    if polar is True:
      f = open(name+"polar.avgheight.pdb","w")
    elif polar is False:
      f = open(name+"cart.avgheight.pdb","w")
    print('CRYST1  150.000  150.000  110.000  90.00  90.00  90.00 P 1           1', file=f)
    serial = 1
    for field in field_list:

      #read in protein helix coordinates
      if protein_onoff == 1:
        protein_coords = np.loadtxt(name+"_helcoords_"+field+".dat",skiprows=1)
        protein = []
        for i in range(10):
          protein.append(protein_coords[i])
      else:
        protein = False

      if readbeads == 0:
        serial = output_analysis(name, field, protein, 3, False, f, serial, polar)
      elif readbeads == 1:
        serial = output_analysis(name, field, protein, 3, False, f, serial, polar)
        if field != "zzero":
          for bead in bead_dict[name]:
            serial = output_analysis(name, field, protein, 3, bead, f, serial, polar)
    print('END', file=f)
    f.close()
  #for name in name_list:
  #  for field in ['zone', 'ztwo']:
  #    gen_avg_tilt(name, field, polar)
