import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os 
from utils import *
#from code_review2 import *


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

def analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims, field_list, scale_dict):
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
        plot_maker(dim1vals, dim2vals, avgHeight, sys_name, field, scale_dict["height_max"], scale_dict["height_min"], inclusion, "avgHeight", bead, polar)

        #save as file for debugging / analysis AND make PDB!
        np.save('npy/'+sys_name+'.'+field+'.'+bead+'.'+coordsys+'.height.npy', height)
        np.savetxt('dat/'+sys_name+'.'+field+'.'+bead+'.'+coordsys+'.avgheight.dat', avgHeight,delimiter = ',',fmt='%10.5f')
        serial = Make_surface_PDB(avgHeight, sys_name, field, d1, d2, pdb, serial, bead, polar)
        print(sys_name+' '+bead+' '+field+" height done!")

      calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb)

    #print last line of pdb file
    print('END', file=pdb)

  return 