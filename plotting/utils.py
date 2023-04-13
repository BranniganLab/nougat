import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 


def strip_blank_lines(f):
  for l in f:
    line = l.strip()
    if line:
      yield line


def read_config(path):
  config_dict = {}
  with open(path, "r+") as config_file:
    for line in strip_blank_lines(config_file):
      if line.startswith("#") is True:
        continue
      else:
        line = line.partition('#')[0]
        if line.rstrip():
          key = line.partition('=')[0].strip()
          value = line.partition('=')[2].strip()
          config_dict.update({key : value})

  return config_dict


def find_first_val(l):
  for value in l:
    if np.isnan(value):
      continue
    else:
      return value 
  return np.nan


def find_last_val(l):
  last_index = len(l)-1
  while last_index >= 0:
    if not np.isnan(l[last_index]):
      return l[last_index]
    else:
      last_index -= 1
  return "this didnt work"


def filename_generator(sys_name, lipid_name, field, beadname, coordsys, measure, dtype):
  if measure == "height" or measure == "curvature" or measure == "Kcurvature" or measure == "thickness":
    if dtype == "dat":
      if measure == "thickness":
        fullmeasure = "normthickness"
      else:
        fullmeasure = "avg"+measure
    elif dtype == "npy":
      if measure == "curvature":
        fullmeasure = "meancurvature"
      elif measure == "Kcurvature":
        fullmeasure = "gausscurvature"
      else:
        fullmeasure = measure 
    filename = sys_name+"."+field+"."+beadname+"."+coordsys+"."+fullmeasure+"."+dtype
  elif measure == "density":
    if dtype == "dat":
      filename = sys_name+"."+lipid_name+"."+field+"."+coordsys+".avg"+measure+"."+dtype
    elif dtype == "npy":
      filename = sys_name+"."+lipid_name+"."+field+"."+coordsys+"."+measure+"."+dtype
  elif measure == "tail1" or measure == "tail0":
    if dtype == "dat":
      filename = sys_name+"."+lipid_name+"."+measure+"."+field+"."+coordsys+".avgOrder."+dtype
    elif dtype == "npy":
      filename = sys_name+"."+lipid_name+"."+measure+"."+field+"."+coordsys+".order."+dtype
  if measure not in ["height", "curvature", "Kcurvature", "thickness", "density", "tail1", "tail0"]:
    raise RuntimeWarning('You used \"'+measure+'\" but that is not an allowed measurement')
  return filename


def calc_avg_over_time(matrix_data):
  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    avg=np.nanmean(matrix_data, axis=2)
    return avg


def bin_prep(sys_name, beadnames, coordsys, density):

  sample_data = np.genfromtxt('tcl_output/'+sys_name+'.zone.'+beadnames+'.'+coordsys+'.height.dat',missing_values='nan',filling_values=np.nan)

  N1_bins, d1, N2_bins, d2, Nframes, min_val = dimensions_analyzer(sample_data, coordsys)
  
  #prep plot dimensions
  dim1 = sample_data[0:N1_bins,0]
  dim1 = np.append(dim1, sample_data[N1_bins-1,1])
  if coordsys == "polar":
    dim2 = np.linspace(0,2*np.pi,N2_bins+1)
  elif coordsys == "cart":
    dim2 = np.linspace(0,N2_bins+1,N2_bins+1)
  dim1vals,dim2vals=np.meshgrid(dim1, dim2, indexing='ij')

  if density == "ON":
    #save an array that represents the area per bin for normalizing density later
    save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name)

  return [N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals]


def save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name):
  
  areas = np.ones([N1_bins,N2_bins])

  areas = areas * d1 * d2 
  if coordsys == "polar":
    for row in range(N1_bins):
      dist_to_center = min_val + row*d1 + d1/2.0
      areas[row,:] = areas[row,:]*dist_to_center
  np.save('npy/'+sys_name+"."+coordsys+".areas.npy", areas)


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
  os.chdir('tcl_output')
  names_dict = {}
  names_dict['species_list'] = [] 
  names_dict['beads_list'] = []

  files = glob.glob(sys_name+'*.zone.'+coordsys+'.tilt.dat')

  for filename in files:
    namefields = filename.split('.')
    species = namefields[1]
    tail = namefields[2]

    if species not in names_dict['species_list']:
      names_dict['species_list'].append(species)
      names_dict[species] = [tail]
    else:
      names_dict[species].append(tail)

  files = glob.glob(sys_name+'.zone.*.'+coordsys+'.height.dat')
  
  for filename in files:
    beadstart = len(sys_name)+6
    beadend = filename.find('.'+coordsys)
    beadname = filename[beadstart:beadend]
    if beadname not in names_dict['beads_list']:
      names_dict['beads_list'].append(beadname)

    os.chdir('..')

  return names_dict


def plot_maker(dim1vals, dim2vals, data, name, field, Vmax, Vmin, protein, dataname, bead, coordsys):
  fig = plt.figure()

  if coordsys == "polar":
    ax = plt.subplot(projection="polar")
    c = plt.pcolormesh(dim2vals,dim1vals,data,cmap="RdBu_r",zorder=0,vmax=Vmax,vmin=Vmin)
  elif coordsys == "cart":
    ax = plt.subplot()
    c = plt.pcolormesh(dim1vals,dim2vals,data,cmap="RdBu_r",zorder=0,vmax=Vmax,vmin=Vmin)
  else:
    print("something's wrong with coordsys")
  
  #cbar = plt.colorbar(c)

  if protein is not False:
    print(protein)
    for i in range(0,10,2):
      protein[i+1] = np.deg2rad(protein[i+1])
      if coordsys == "cart":
        protein[i], protein[i+1] = convert_to_cart(protein[i],protein[i+1])
      plt.scatter(protein[i+1],protein[i],c="black",linewidth=4,zorder=2)
    circle1 = plt.Circle((0,0),28.116, transform=ax.transData._b, color='black',linestyle='dashed',linewidth=4,fill=False)
    if field == "zone":
      ax.add_artist(circle1)

  plt.axis('off')
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  #fig.set_size_inches(6,6)

  if bead is False:
    plt.savefig('pdf/'+name+"_"+field+"_"+coordsys+"_"+dataname+".pdf", dpi = 700)
  else:
    plt.savefig('pdf/'+name+"_"+bead+"_"+field+"_"+coordsys+"_"+dataname+".pdf", dpi = 700)
  plt.clf()
  plt.close()


def convert_to_cart(rval, thetaval):
  xval = rval*np.cos(thetaval)
  yval = rval*np.sin(thetaval)
  return xval, yval


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


def dimensions_analyzer(data, coordsys):
  # figure out how many radial or x bins there are
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
  

  #figure out how many azimuthal or y bins there are
  N2_bins = len(data[0,:]) - 2
  
  #figure out how many frames there are in the traj
  Nframes = int(len(data[:,0])/N1_bins)

  if coordsys == "polar":
    d1 = data[0,1] - data[0,0]
    d2 = (np.pi*2)/N2_bins
  elif coordsys == "cart":
    #compute average d1, assume d2 is the same
    d1list = []
    for row in range(Nframes):
      d1list.append(data[row*N1_bins,1])
    d1 = np.mean(d1list)
    d2 = d1

  return N1_bins, d1, N2_bins, d2, Nframes, match_value

