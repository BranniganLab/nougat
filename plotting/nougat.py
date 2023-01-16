import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 
from height import *
from thickness import *
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