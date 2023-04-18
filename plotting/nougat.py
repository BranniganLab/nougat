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


def run_nougat(sys_name, polar, inclusion_drawn, config_dict):
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

  field_list = ["zone","ztwo", "zzero"]

  #figure out all the file names that you'll need to fetch
  names_dict = read_log(sys_name, coordsys)

  #get data dimensions and prep plots from one of your trajectories
  dims = bin_prep(sys_name, names_dict['beads_list'][0], coordsys, "ON")

  #analyze height
  analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims, field_list, config_dict)

  for bead in names_dict['beads_list']:
    calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims, config_dict)
    calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims, field_list, config_dict)
  
  calculate_density(sys_name, names_dict, coordsys, inclusion, polar, dims, config_dict)
  #calculate_order(sys_name, names_dict, coordsys, inclusion, polar, dims, config_dict)
  #calculate_tilt(sys_name, names_dict, coordsys, inclusion, polar, dims, config_dict)

  calc_epsilon_and_H_terms(sys_name, ".", coordsys)


if __name__ == "__main__": 
  parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
  parser.add_argument("sys_name", help="what system did you name this?")
  parser.add_argument("config", help="what config file should nougat use?")
  parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl in polar coordinates")
  parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
  args = parser.parse_args()

  config_dict = read_config(args.config)

  run_nougat(args.sys_name, args.polar, args.inclusion, config_dict)

  print("Thank you for using nougat!")