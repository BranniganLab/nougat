import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 
from utils import *
#from code_review2 import *

def measure_t0(zone, ztwo, coordsys):
  thickness = zone-ztwo

  avgthickness = calc_avg_over_time(thickness)

  if coordsys == "cart":
    leftcol = np.mean(avgthickness[:,0])
    rightcol =  np.mean(avgthickness[:,-1])
    toprow =  np.mean(avgthickness[0,:])
    botrow =   np.mean(avgthickness[-1,:])
    avgt0 = (leftcol+rightcol+toprow+botrow)/4.0
  elif coordsys == "polar":
    avgt0 = np.mean(avgthickness[-1:])

  avgt0 = avgt0/2.0

  return avgt0

def calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims, scale_dict):
  N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims 
  zone = np.load('npy/'+sys_name+'.zone.'+bead+'.'+coordsys+'.height.npy')
  ztwo = np.load('npy/'+sys_name+'.ztwo.'+bead+'.'+coordsys+'.height.npy')
  zzero = np.load('npy/'+sys_name+'.zzero.'+bead+'.'+coordsys+'.height.npy')

  outer_leaflet = zone-zzero
  inner_leaflet = zzero-ztwo

  avgt0 = measure_t0(zone, ztwo, coordsys)

  avgouter = calc_avg_over_time(outer_leaflet)/avgt0
  avginner = calc_avg_over_time(inner_leaflet)/avgt0

  #make plots!
  plot_maker(dim1vals, dim2vals, avgouter, sys_name, 'zone', scale_dict["thick_max"], scale_dict["thick_min"], inclusion, "avgThickness", bead, coordsys)
  plot_maker(dim1vals, dim2vals, avginner, sys_name, 'zone', scale_dict["thick_max"], scale_dict["thick_min"], inclusion, "avgThickness", bead, coordsys)

  #save as file for debugging / analysis AND make PDB!
  np.save('npy/'+sys_name+'.zone.'+bead+'.'+coordsys+'.thickness.npy', outer_leaflet)
  np.save('npy/'+sys_name+'.ztwo.'+bead+'.'+coordsys+'.thickness.npy', inner_leaflet)
  np.savetxt('dat/'+sys_name+'.zone.'+bead+'.'+coordsys+'.avgthickness.dat', avgouter,delimiter = ',',fmt='%10.5f')
  np.savetxt('dat/'+sys_name+'.ztwo.'+bead+'.'+coordsys+'.avgthickness.dat', avginner,delimiter = ',',fmt='%10.5f')

  print(sys_name+' '+bead+" thickness done!")

