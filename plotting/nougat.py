import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 
from curvature import *
from density import *
from tilt import *
from order import *
from utils import *
#from code_review2 import *

# These determine the scale in your image files
# adjust as needed
height_min = -60
height_max = 60
mean_curv_min = 0.01
mean_curv_max = -0.01
gauss_curv_min = -0.001
gauss_curv_max = 0.001
density_min = 0
density_max = 2
thick_min = 0
thick_max = 2
order_min = 0. 
order_max = .6

field_list = ["zone","ztwo", "zzero"]



def calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims 
  zone = np.load('npy/'+sys_name+'.zone.'+bead+'.'+coordsys+'.height.npy')
  ztwo = np.load('npy/'+sys_name+'.ztwo.'+bead+'.'+coordsys+'.height.npy')

  zplus=(zone+ztwo)/2

  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    avgzplus=np.nanmean(zplus, axis=2)

  #make plots!
  plot_maker(dim1vals, dim2vals, avgzplus, sys_name, 'zplus', height_max, height_min, inclusion, "avgHeight", bead, polar)

  #save as file for debugging / analysis AND make PDB!
  np.save('npy/'+sys_name+'.zplus.'+bead+'.'+coordsys+'.height.npy', zplus)
  np.savetxt('dat/'+sys_name+'.zplus.'+bead+'.'+coordsys+'.avgheight.dat', avgzplus,delimiter = ',',fmt='%10.5f')
  serial = Make_surface_PDB(avgzplus, sys_name, 'zplus', d1, d2, pdb, serial, bead, polar)
  print(sys_name+' '+bead+" zplus height done!")


def measure_t0(zone, ztwo):

  thickness = zone-ztwo

  avgthickness = calc_avg_over_time(thickness)

  leftcol = np.mean(avgthickness[:,0])
  rightcol =  np.mean(avgthickness[:,-1])
  toprow =  np.mean(avgthickness[0,:])
  botrow =   np.mean(avgthickness[-1,:])


  avgt0 = (leftcol+rightcol+toprow+botrow)/4.0

  avgt0 = avgt0/2.0

  return avgt0

def calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims 
  zone = np.load('npy/'+sys_name+'.zone.'+bead+'.'+coordsys+'.height.npy')
  ztwo = np.load('npy/'+sys_name+'.ztwo.'+bead+'.'+coordsys+'.height.npy')
  zzero = np.load('npy/'+sys_name+'.zzero.'+bead+'.'+coordsys+'.height.npy')

  outer_leaflet = zone-zzero
  inner_leaflet = zzero-ztwo

  avgt0 = measure_t0(zone, ztwo)

  avgouter = calc_avg_over_time(outer_leaflet)/avgt0
  avginner = calc_avg_over_time(inner_leaflet)/avgt0

  #make plots!
  plot_maker(dim1vals, dim2vals, avgouter, sys_name, 'outer', thick_max, thick_min, inclusion, "avgThickness", bead, polar)
  plot_maker(dim1vals, dim2vals, avginner, sys_name, 'inner', thick_max, thick_min, inclusion, "avgThickness", bead, polar)

  #save as file for debugging / analysis AND make PDB!
  np.save('npy/'+sys_name+'.outer.'+bead+'.'+coordsys+'.thickness.npy', outer_leaflet)
  np.save('npy/'+sys_name+'.inner.'+bead+'.'+coordsys+'.thickness.npy', inner_leaflet)
  np.savetxt('dat/'+sys_name+'.outer.'+bead+'.'+coordsys+'.avgthickness.dat', avgouter,delimiter = ',',fmt='%10.5f')
  np.savetxt('dat/'+sys_name+'.inner.'+bead+'.'+coordsys+'.avgthickness.dat', avginner,delimiter = ',',fmt='%10.5f')

  print(sys_name+' '+bead+" thickness done!")


def analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims):
  serial = 1

  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims 

  pdbname = sys_name+"."+coordsys+".avgheight.pdb"

  with open(pdbname,"w") as pdb:
    
    #print first line of pdb file
    print('CRYST1  150.000  150.000  110.000  90.00  90.00  90.00 P 1           1', file=pdb)

    #do height analysis
    for bead in names_dict['beads_list']:
      for field in field_list:

        #import traj values
        height_data = np.genfromtxt('tcl_output/'+sys_name+'.'+field+'.'+bead+'.'+coordsys+'.height.dat',missing_values='nan',filling_values=np.nan)

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
        np.save('npy/'+sys_name+'.'+field+'.'+bead+'.'+coordsys+'.height.npy', height)
        np.savetxt('dat/'+sys_name+'.'+field+'.'+bead+'.'+coordsys+'.avgheight.dat', avgHeight,delimiter = ',',fmt='%10.5f')
        serial = Make_surface_PDB(avgHeight, sys_name, field, d1, d2, pdb, serial, bead, polar)
        print(sys_name+' '+bead+' '+field+" height done!")

      calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb)

    #print last line of pdb file
    print('END', file=pdb)

  return 

def run_nougat(sys_name, polar, inclusion_drawn):
  cwd = os.getcwd()

  for filetype in ["npy", "dat", "pdf"]:
    dirname = os.path.join(cwd, filetype)
    try:
      os.mkdir(dirname) 
    except OSError as error:
      continue

  if inclusion_drawn is True:
    inclusion = add_inclusion(name, field_list) #this proc doesn't exist yet!
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
  analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims)

  for bead in names_dict['beads_list']:
    calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims)
    calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims)
  
  calculate_density(sys_name, names_dict, coordsys, inclusion, polar, dims)
  #calculate_order(sys_name, names_dict, coordsys, inclusion, polar, dims)
  #calculate_tilt(sys_name, names_dict, coordsys, inclusion, polar, dims)


if __name__ == "__main__": 
#  run_nougat("lgPO", False, False)
#  for system in ["lgDT0", "lgDP", "lgDX", "lgDB", "lgDL"]:
#  for system in ["lgDB"]:
#    run_nougat(system, False, False)
#  for system in ["DT", "DY", "DL", "DO", "DP", "PO", "DG", "DB", "DX"]: 
#    os.chdir(system)
  run_nougat("lgPO", False, False)
    #os.chdir('newleaf_polar')
    #run_nougat("lgPO", True, False)
#    os.chdir('..')