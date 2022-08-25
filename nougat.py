import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob

sys_name = 'lgPOtest'
inclusion_drawn = 0
polar = False


# These determine the scale in your image files
# adjust as needed
height_min = -60
height_max = 60
mean_curv_min = 0.05
mean_curv_max = -0.05
gauss_curv_min = -0.05
gauss_curv_max = 0.05
density_min = 0
density_max = 2
thick_min = 5
thick_max = 15
order_min = -1
order_max = 1

field_list = ["zone","ztwo"]

def dimensions_analyzer(data, polar):
  # figure out how many radial bins there are
  counter = 1
  flag = True
  match_value = data[0,0]
  while(flag==True):
    try:
      if data[counter,0] == match_value:
        flag = False
      else:
        counter = counter+1
    # what if there is only 1 frame? Will raise IndexError 
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

  return N1_bins, d1, N2_bins, d2, Nframes, match_value


def empty_neighbor_test(curvature_inputs):
  data = np.isnan(curvature_inputs)
  nan_test = np.array(data, copy=True)
  nan_test2 = np.array(data, copy=True)
  knan_test = np.array(data, copy=True)
  
  shape = np.shape(data)
  dim1 = shape[0]
  dim2 = shape[1]
  dim3 = shape[2]

  for frm in range(dim3):
    for row in range(1,dim1-1):
      for col in range(1,dim2-1):
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
  nan_test[dim1-1,:,:] = True
  nan_test[:,0,:] = True
  nan_test[:,dim2-1,:] = True

  knan_test[0,:,:] = True
  knan_test[dim1-1,:,:] = True
  knan_test[:,0,:] = True
  knan_test[:,dim2-1,:] = True

  return nan_test, knan_test


def plot_maker(dim1vals, dim2vals, data, name, field, Vmax, Vmin, protein, dataname, bead, polar):
  fig = plt.figure()

  if polar is True:
    ax = plt.subplot(projection="polar")
    c = plt.pcolormesh(dim2vals,dim1vals,data,cmap="RdBu_r",zorder=0,vmax=Vmax,vmin=Vmin)
    coordsys = "polar"
  elif polar is False:
    ax = plt.subplot()
    c = plt.pcolormesh(dim1vals,dim2vals,data,cmap="RdBu_r",zorder=0,vmax=Vmax,vmin=Vmin)
    coordsys = "cart"
    
  
  cbar = plt.colorbar(c)

  if protein is not False:
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

  if bead is False:
    plt.savefig(name+"_"+field+"_"+coordsys+"_"+dataname+".pdf", dpi = 700)
  else:
    plt.savefig(name+"_"+bead+"_"+field+"_"+coordsys+"_"+dataname+".pdf", dpi = 700)
  plt.clf()
  plt.close()

def convert_to_cart(rval, thetaval):
  xval = rval*np.cos(thetaval)
  yval = rval*np.sin(thetaval)
  return xval, yval

def measure_curvature_cart(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims)
  
  #mean curvature: Hxx + Hyy
  #gaussian curvature: HxxHyy - Hxy^2

  for frm in range(Nframes):
    for row in range(N1_bins+2):
      for col in range(N2_bins+2):
        if knan_test[row,col,frm] == False:

          del2x = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
          del2x = del2x / (d1**2)

          del2y = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
          del2y = del2y / (d2**2)

          delxy = (curvature_inputs[row+1,col+1,frm] - curvature_inputs[row+1,col,frm] - curvature_inputs[row,col+1,frm] + curvature_inputs[row,col,frm] - curvature_inputs[row-1,col,frm] - curvature_inputs[row,col-1,frm] + curvature_inputs[row-1,col-1,frm])
          delxy = delxy/(2*d1*d2)

          #delxy = curvature_inputs[row+1,col+1,frm] - curvature_inputs[row+1,col-1,frm] - curvature_inputs[row-1,col+1,frm] + curvature_inputs[row-1,col-1,frm]
          #delxy = delxy / (4*d1*d2)

          delx = (curvature_inputs[row+1,col,frm] - curvature_inputs[row-1,col,frm])/(2*d1)

          dely = (curvature_inputs[row,col+1,frm] - curvature_inputs[row,col-1,frm])/(2*d2)

          normalization_factor = np.sqrt(1+delx**2+dely**2)
          norm_vec_x = -1*delx/normalization_factor
          norm_vec_y = -1*dely/normalization_factor
          norm_vec_z = 1/normalization_factor

          curvature_outputs[row,col,frm] = del2x + del2y
          kgauss_outputs[row,col,frm] = del2x*del2y - delxy**2
          normal_vector_outputs[row,col*3,frm] = norm_vec_x
          normal_vector_outputs[row,col*3+1,frm] = norm_vec_y
          normal_vector_outputs[row,col*3+2,frm] = norm_vec_z
        
        elif nan_test[row,col,frm] == False:
          del2x = curvature_inputs[row-1,col,frm] + curvature_inputs[row+1,col,frm] - 2*curvature_inputs[row,col,frm]
          del2x = del2x / d1**2

          del2y = curvature_inputs[row,col-1,frm] + curvature_inputs[row,col+1,frm] - 2*curvature_inputs[row,col,frm]
          del2y = del2y / d2**2

          delx = (curvature_inputs[row+1,col,frm] - curvature_inputs[row-1,col,frm])/(2*d1)

          dely = (curvature_inputs[row,col+1,frm] - curvature_inputs[row,col-1,frm])/(2*d2)

          normalization_factor = np.sqrt(1+delx**2+dely**2)
          norm_vec_x = -1*delx/normalization_factor
          norm_vec_y = -1*dely/normalization_factor
          norm_vec_z = 1/normalization_factor

          curvature_outputs[row,col,frm] = del2x + del2y
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

def measure_curvature_polar(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims)
  
  #mean curvature: h_rr + 1/r(h_r) + 1/r**2(h_thetatheta)
  #gaussian curvature: 1/r(h_r*h_rr) + 2/r**3(h_rtheta*h_theta) - 1/r**4(h_theta**2) - 1/r**2(h_rtheta**2 - h_rr*h_thetatheta)

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
          delrdeltheta = (curvature_inputs[row+1,col+1,frm] - curvature_inputs[row+1,col,frm] - curvature_inputs[row,col+1,frm] + curvature_inputs[row,col,frm] - curvature_inputs[row-1,col,frm] - curvature_inputs[row,col-1,frm] + curvature_inputs[row-1,col-1,frm])
          delrdeltheta = delrdeltheta/(2*d1*d2)

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
          curvature_outputs[row,col,frm] = del2r + c1*delr + c2*del2theta
          kgauss_outputs[row,col,frm] = c1*delr*del2r + 2*c3*delrdeltheta*deltheta - c4*deltheta**2 - c2*(delrdeltheta**2-del2r*del2theta)
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
  N1_bins, d1, N2_bins, d2, Nframes, min_val = dimensions_analyzer(height_data, polar)

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

def calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 
  zone = np.load(sys_name+'.zone.'+bead+'.'+coordsys+'.height.npy')
  ztwo = np.load(sys_name+'.ztwo.'+bead+'.'+coordsys+'.height.npy')

  zplus=(zone+ztwo)/2

  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    avgzplus=np.nanmean(zplus, axis=2)

  #make plots!
  plot_maker(dim1vals, dim2vals, avgzplus, sys_name, 'zplus', height_max, height_min, inclusion, "avgHeight", bead, polar)

  #save as file for debugging / analysis AND make PDB!
  np.save(sys_name+'.zplus.'+bead+'.'+coordsys+'.height.npy', zplus)
  np.savetxt(sys_name+'.zplus.'+bead+'.'+coordsys+'.avgheight.dat', avgzplus,delimiter = ',',fmt='%10.5f')
  serial = Make_surface_PDB(avgzplus, sys_name, 'zplus', d1, d2, pdb, serial, bead, polar)
  print(sys_name+' '+bead+" zplus height done!")

def calc_avg_over_time(matrix_data):
  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    avg=np.nanmean(matrix_data, axis=2)
    return avg

def init_curvature_data(height, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims)

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

  if polar is True:
    #wrap the inputs in the theta direction for calculating curvature
    curvature_inputs[:,1:(N2_bins+1),:] = height
    curvature_inputs[:,0,:] = curvature_inputs[:,N2_bins,:]
    curvature_inputs[:,(N2_bins+1),:] = curvature_inputs[:,1,:]
  elif polar is False:
    #if cartesian, wrap in both directions
    curvature_inputs[1:(N1_bins+1),1:(N2_bins+1),:] = height
    curvature_inputs[:,0,:] = curvature_inputs[:,N2_bins,:]
    curvature_inputs[:,(N2_bins+1),:] = curvature_inputs[:,1,:]
    curvature_inputs[0,:,:] = curvature_inputs[N1_bins,:,:]
    curvature_inputs[(N1_bins+1),:,:] = curvature_inputs[1,:,:]
    #and fill in the corners
    curvature_inputs[0,0,:] = curvature_inputs[N1_bins,N2_bins,:]
    curvature_inputs[N1_bins+1,N2_bins+1,:] = curvature_inputs[1,1,:]
    curvature_inputs[0,N2_bins+1,:] = curvature_inputs[N1_bins,1,:]
    curvature_inputs[N1_bins+1,0,:] = curvature_inputs[1,N2_bins,:]

  return curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs

def calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims)
  for field in ["zone", "ztwo", "zzero", "zplus"]: 
    field_height = np.load(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.height.npy')

    curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs = init_curvature_data(field_height, polar, dims)

    #if a bin is empty, you can't (nicely) measure the curvature of its neighbors
    nan_test, knan_test = empty_neighbor_test(curvature_inputs)
    
    #measure the laplacian and gaussian curvatures
    if polar is True:
      curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_polar(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims)
    elif polar is False:
      curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_cart(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims)

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
    avgcurvature = calc_avg_over_time(meancurvature)
    avgkcurvature = calc_avg_over_time(kcurvature)

    #make plots!
    plot_maker(dim1vals, dim2vals, avgkcurvature, sys_name, field, gauss_curv_max, gauss_curv_min, inclusion, "gausscurvature", bead, polar)
    plot_maker(dim1vals, dim2vals, avgcurvature, sys_name, field, mean_curv_max, mean_curv_min, inclusion, "curvature", bead, polar)

    #save as files for debugging / analysis
    np.savetxt(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.avgcurvature.dat',avgcurvature,delimiter = ',',fmt='%10.7f')
    np.savetxt(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.avgKcurvature.dat',avgkcurvature,delimiter = ',',fmt='%10.7f')
    np.save(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.meancurvature.npy',meancurvature)
    np.save(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.gausscurvature.npy',kcurvature)
    np.save(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.normal_vectors.npy',normal_vectors)
    
    print(sys_name+' '+bead+' '+field+" curvatures done!")

def calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 
  zone = np.load(sys_name+'.zone.'+bead+'.'+coordsys+'.height.npy')
  ztwo = np.load(sys_name+'.ztwo.'+bead+'.'+coordsys+'.height.npy')
  zzero = np.load(sys_name+'.zzero.'+bead+'.'+coordsys+'.height.npy')

  outer_leaflet = zone-zzero
  inner_leaflet = zzero-ztwo

  avgouter = calc_avg_over_time(outer_leaflet)
  avginner = calc_avg_over_time(inner_leaflet)

  #make plots!
  plot_maker(dim1vals, dim2vals, avgouter, sys_name, 'outer', thick_max, thick_min, inclusion, "avgThickness", bead, polar)
  plot_maker(dim1vals, dim2vals, avginner, sys_name, 'inner', thick_max, thick_min, inclusion, "avgThickness", bead, polar)

  #save as file for debugging / analysis AND make PDB!
  np.save(sys_name+'.outer.'+bead+'.'+coordsys+'.thickness.npy', outer_leaflet)
  np.save(sys_name+'.inner.'+bead+'.'+coordsys+'.thickness.npy', inner_leaflet)
  np.savetxt(sys_name+'.outer.'+bead+'.'+coordsys+'.avgthickness.dat', avgouter,delimiter = ',',fmt='%10.5f')
  np.savetxt(sys_name+'.inner.'+bead+'.'+coordsys+'.avgthickness.dat', avginner,delimiter = ',',fmt='%10.5f')

  print(sys_name+' '+bead+" thickness done!")

def calculate_density(sys_name, names_dict, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 
  areas = np.load(sys_name+".areas.npy")

  for species in names_dict['species_list']:
    zone = np.genfromtxt(sys_name+'.'+species+'.zone.'+coordsys+'.density.dat',missing_values='nan',filling_values="0")
    ztwo = np.genfromtxt(sys_name+'.'+species+'.ztwo.'+coordsys+'.density.dat',missing_values='nan',filling_values="0")

    #create a new array that has each frame in a different array level
    density_up = np.zeros((N1_bins, N2_bins, Nframes))
    density_down = np.zeros((N1_bins, N2_bins, Nframes))
    for frm in range(Nframes):
      density_up[:,:,frm] = zone[frm*N1_bins:(frm+1)*N1_bins,2:]
      density_down[:,:,frm] = ztwo[frm*N1_bins:(frm+1)*N1_bins,2:]

    avgouter = calc_avg_over_time(density_up)
    avginner = calc_avg_over_time(density_down)

    found = False
    with open(sys_name+'.density.normfactor.dat', 'r') as normfile:
      norms = normfile.readlines()
      while found is False:
        for line in norms:
          vals = line.strip().split(' ')
          if vals[0] == species:
            normfactor = float(vals[1])
            found = True

    #normalize
    avgouter = (avgouter * normfactor) / areas
    avginner = (avginner * normfactor) / areas 

    #make plots!
    plot_maker(dim1vals, dim2vals, avgouter, sys_name, species+'.outer', density_max, density_min, inclusion, "avgDensity", False, polar)
    plot_maker(dim1vals, dim2vals, avginner, sys_name, species+'.inner', density_max, density_min, inclusion, "avgDensity", False, polar)

    #save as file for debugging / analysis 
    np.save(sys_name+'.'+species+'.zone.'+coordsys+'.density.npy', density_up)
    np.save(sys_name+'.'+species+'.ztwo.'+coordsys+'.density.npy', density_down)
    np.savetxt(sys_name+'.'+species+'.zone.'+coordsys+'.avgdensity.dat', avgouter,delimiter = ',',fmt='%10.5f')
    np.savetxt(sys_name+'.'+species+'.ztwo.'+coordsys+'.avgdensity.dat', avginner,delimiter = ',',fmt='%10.5f')

    print(sys_name+' '+species+" density done!")

  if len(names_dict['species_list']) > 1:
    calculate_total_density(sys_name, names_dict, coordsys, inclusion, polar, dims)
  elif len(names_dict['species_list']) < 1:
    print("Something is wrong with species_list!")

def calculate_total_density(sys_name, names_dict, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 

  for leaflet in ['zone', 'ztwo']:
    tot_density = np.zeros((N1_bins, N2_bins, Nframes))
    
    for species in names_dict['species_list']:
      dens_per_species = np.load(sys_name+'.'+species+'.'+leaflet+'.'+coordsys+'.density.npy')
      tot_density = tot_density + dens_per_species
      #normalize?

    tot_avg_density = np.mean(tot_density, axis=2)

    #make plots!
    plot_maker(dim1vals, dim2vals, tot_avg_density, sys_name, species+'.'+leaflet, density_max, density_min, inclusion, "totAvgDensity", False, polar)

    #save as file for debugging / analysis
    np.save(sys_name+'.'+species+'.'+leaflet+'.'+coordsys+'.totalDensity.npy', tot_density)
    np.savetxt(sys_name+'.'+species+'.'+leaflet+'.'+coordsys+'.totalAvgDensity.dat', tot_avg_density, delimiter = ',', fmt='%10.5f')

  print("Total density done!")

def calculate_order(sys_name, names_dict, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 
  
  for species in names_dict['species_list']:
    for tail in names_dict[species]:
      zone = np.genfromtxt(sys_name+'.'+species+'.'+tail+'.zone.'+coordsys+'.order.dat',missing_values='nan',filling_values=np.nan)
      ztwo = np.genfromtxt(sys_name+'.'+species+'.'+tail+'.ztwo.'+coordsys+'.order.dat',missing_values='nan',filling_values=np.nan)

      #create a new array that has each frame in a different array level
      order_up = np.zeros((N1_bins, N2_bins, Nframes))
      order_down = np.zeros((N1_bins, N2_bins, Nframes))
      for frm in range(Nframes):
        order_up[:,:,frm] = zone[frm*N1_bins:(frm+1)*N1_bins,2:]
        order_down[:,:,frm] = ztwo[frm*N1_bins:(frm+1)*N1_bins,2:]

      avgouter = calc_avg_over_time(order_up)
      avginner = calc_avg_over_time(order_down)

      #make plots!
      plot_maker(dim1vals, dim2vals, avgouter, sys_name, species+'.'+tail+'.zone', order_max, order_min, inclusion, "avgOrder", False, polar)
      plot_maker(dim1vals, dim2vals, avginner, sys_name, species+'.'+tail+'.ztwo', order_max, order_min, inclusion, "avgOrder", False, polar)

      #save as file for debugging / analysis 
      np.save(sys_name+'.'+species+'.'+tail+'.zone.'+coordsys+'.order.npy', order_up)
      np.save(sys_name+'.'+species+'.'+tail+'.ztwo.'+coordsys+'.order.npy', order_down)
      np.savetxt(sys_name+'.'+species+'.'+tail+'.zone.'+coordsys+'.avgOrder.dat', avgouter,delimiter = ',',fmt='%10.5f')
      np.savetxt(sys_name+'.'+species+'.'+tail+'.ztwo.'+coordsys+'.avgOrder.dat', avginner,delimiter = ',',fmt='%10.5f')

      print(sys_name+' '+species+" "+tail+" order done!")

    #if len(names_dict[species]) > 1:
    #  calculate_total_order(sys_name, species, names_dict, coordsys, inclusion, polar, dims)
    #elif len(names_dict[species]) < 1:
    #  print("Something is wrong with your tails list!")

  #if len(names_dict['species_list']) > 1:
  #  calculate_total_order(sys_name, "all", names_dict, coordsys, inclusion, polar, dims)
  #elif len(names_dict['species_list']) < 1:
  #  print("Something is wrong with your species_list!")

def calculate_total_order(sys_name, species, names_dict, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 

  if species != "all":
    for leaflet in ['zone','ztwo']:
      tot_order = np.zeros((N1_bins, N2_bins, Nframes))

      for tail in names_dict[species]:
        order_per_tail = np.load(sys_name+'.'+species+'.'+tail+'.'+leaflet+'.'+coordsys+'.order.npy')
        tot_order = tot_order + order_per_tail
        #weight average by density!



def bin_prep(sys_name, names_dict, coordsys, polar):
  sample_data = np.genfromtxt(sys_name+'.zone.'+names_dict['beads_list'][0]+'.'+coordsys+'.height.dat',missing_values='nan',filling_values=np.nan)
  N1_bins, d1, N2_bins, d2, Nframes, min_val = dimensions_analyzer(sample_data, polar)
  
  #prep plot dimensions
  dim1 = sample_data[0:N1_bins,0]
  dim1 = np.append(dim1, sample_data[N1_bins-1,1])
  if polar is True:
    dim2 = np.linspace(0,2*np.pi,N2_bins+1)
  elif polar is False:
    dim2 = np.linspace(0,N2_bins+1,N2_bins+1)
  dim1vals,dim2vals=np.meshgrid(dim1, dim2, indexing='ij')

  #save an array that represents the area per bin for normalizing density later
  save_areas(N1_bins, d1, N2_bins, d2, min_val, polar, sys_name)

  return [N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals]

def save_areas(N1_bins, d1, N2_bins, d2, min_val, polar, sys_name):
  
  areas = np.ones([N1_bins,N2_bins])

  areas = areas * d1 * d2 
  if polar is True:
    #convert d2 from degrees to radians
    areas = areas * 2 * np.pi / 360
    for row in range(N1_bins):
      dist_to_center = min_val + row*d1 + d1/2.0
      areas[row,:] = areas[row,:]*dist_to_center
  np.save(sys_name+".areas.npy", areas)

def mostly_empty(data_array, N1_bins, N2_bins, Nframes):
  #if a bin only has lipids in it <10% of the time, it shouldn't be considered part of the membrane
  for row in range(N1_bins):
    for col in range(N2_bins):
      zerocount = np.count_nonzero(data_array[row,col,:])
      count = np.count_nonzero(np.isnan(data_array[row,col,:]))
      if (zerocount-count)/Nframes <= .1:
        data_array[row,col,:] = np.nan
  return data_array

def fetch_names(sys_name, coordsys):
  names_dict = {}
  names_dict['species_list'] = [] 
  names_dict['beads_list'] = []

  files = glob.glob(sys_name+'*.zone.*.tilt.dat')

  for filename in files:
    namefields = filename.split('.')
    species = namefields[1]
    tail = namefields[2]

    if species not in names_dict['species_list']:
      names_dict['species_list'].append(species)
      names_dict[species] = [tail]
    else:
      names_dict[species].append(tail)

  files = glob.glob(sys_name+'.zone.*.height.dat')
  
  for filename in files:
    beadstart = len(sys_name)+6
    beadend = filename.find('.'+coordsys)
    beadname = filename[beadstart:beadend]

    if beadname not in names_dict['beads_list']:
      names_dict['beads_list'].append(beadname)

  return names_dict

def unpack_dims(dims):
  return dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6]

def analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims):
  serial = 1

  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = unpack_dims(dims) 

  pdbname = sys_name+"."+coordsys+".avgheight.pdb"

  with open(pdbname,"w") as pdb:
    
    #print first line of pdb file
    print('CRYST1  150.000  150.000  110.000  90.00  90.00  90.00 P 1           1', file=pdb)

    #do height analysis
    for bead in names_dict['beads_list']:
      for field in ['zone', 'ztwo', 'zzero']:

        #import traj values
        height_data = np.genfromtxt(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.height.dat',missing_values='nan',filling_values=np.nan)

        #create a new array that has each frame in a different array level
        height = np.zeros((N1_bins, N2_bins, Nframes))
        for frm in range(Nframes):
          height[:,:,frm] = height_data[frm*N1_bins:(frm+1)*N1_bins,2:]

        #if a bin is occupied <10% of the time, it shouldn't be treated as part of the membrane
        height = mostly_empty(height, N1_bins, N2_bins, Nframes)

        #take the average height over all frames
        avgHeight = calc_avg_over_time(height)

        #make plots!
        plot_maker(dim1vals, dim2vals, avgHeight, sys_name, field, height_max, height_min, inclusion, "avgHeight", bead, polar)

        #save as file for debugging / analysis AND make PDB!
        np.save(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.height.npy', height)
        np.savetxt(sys_name+'.'+field+'.'+bead+'.'+coordsys+'.avgheight.dat', avgHeight,delimiter = ',',fmt='%10.5f')
        serial = Make_surface_PDB(avgHeight, sys_name, field, d1, d2, pdb, serial, bead, polar)
        print(sys_name+' '+bead+' '+field+" height done!")

      calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb)

    #print last line of pdb file
    print('END', file=pdb)

  return 

if __name__ == "__main__": 

  if inclusion_drawn is True:
    inclusion = add_inclusion(name, field_list)
  else:
    inclusion = False

  if polar is True:
    coordsys = 'polar'
  elif polar is False:
    coordsys = 'cart'

  #figure out all the file names that you'll need to fetch
  names_dict = fetch_names(sys_name, coordsys)

  #get data dimensions and prep plots from one of your trajectories
  dims = bin_prep(sys_name, names_dict, coordsys, polar)

  #analyze height
  #analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims)

  #for bead in names_dict['beads_list']:
  
    #calculate thickness
    #calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims)

    #calculate curvature
    #calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims)

  #analyze density
  calculate_density(sys_name, names_dict, coordsys, inclusion, polar, dims)

  #analyze order
  #calculate_order(sys_name, names_dict, coordsys, inclusion, polar, dims)

  #analyze tilts
  #calculate_tilt(sys_name, names_dict, coordsys, inclusion, polar, dims)
