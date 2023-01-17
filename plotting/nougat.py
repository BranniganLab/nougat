import argparse
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


def run_nougat(sys_name, polar, inclusion_drawn, scale_dict):
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
  elif polar is None:
    coordsys = 'cart'

  field_list = ["zone","ztwo", "zzero"]

  #figure out all the file names that you'll need to fetch
  names_dict = fetch_names(sys_name, coordsys)

  #get data dimensions and prep plots from one of your trajectories
  dims = bin_prep(sys_name, names_dict, coordsys, polar)

  #analyze height
  analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims, field_list, scale_dict)

  for bead in names_dict['beads_list']:
    calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims, scale_dict)
    calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims, field_list, scale_dict)
  
  calculate_density(sys_name, names_dict, coordsys, inclusion, polar, dims, scale_dict)
  calculate_order(sys_name, names_dict, coordsys, inclusion, polar, dims, scale_dict)
  #calculate_tilt(sys_name, names_dict, coordsys, inclusion, polar, dims, scale_dict)


if __name__ == "__main__": 
  parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
  parser.add_argument("sys_name", help="what system do you want to run nougat.py on?")
  parser.add_argument("-p", "--polar", help="add this flag if you ran nougat.tcl in polar coordinates")
  parser.add_argument("-i", "--inclusion", help="add this flag if you want your inclusion to show up in images")
  args = parser.parse_args()

  if args.polar != True and args.polar != None:
    print("You tried to specify something in your polar flag. Is that what you meant to do?")
    exit()
  if args.inclusion != True and args.inclusion != None:
    print("You tried to specify something in your inclusion flag. Is that what you meant to do?")
    exit()

  # These determine the scale in your image files
  # adjust as needed
  scale_dict = {
    "height_min" : -60,
    "height_max" : 60,
    "mean_curv_min" : 0.01,
    "mean_curv_max" : -0.01,
    "gauss_curv_min" : -0.001,
    "gauss_curv_max" : 0.001,
    "density_min" : 0,
    "density_max" : 2,
    "thick_min" : 0,
    "thick_max" : 2,
    "order_min" : 0., 
    "order_max" : .6
  }

  run_nougat(args.sys_name, args.polar, args.inclusion, scale_dict)

  print("Thank you for using nougat!")