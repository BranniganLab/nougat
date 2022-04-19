import numpy as np
import warnings
import subprocess as subp
import sys
import os
import glob


#sys_list = ["DL", "DT", "DG", "DX", "PO", "DB", "DY", "DO", "DP", "lgPO"]
sys_list = ["DT"]
#field_list = ['zone', 'ztwo', 'zzero', 'zplus']
field_list = ['zone']

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

dth = (np.pi*2)/30


def coord_format(value):
	rounded = round(value,3)
	leftside,rightside = str(rounded).split('.')
	if len(rightside) < 3:
		rightside = rightside+(' '*(3-len(rightside)))
	if len(leftside) < 4:
		leftside = (' '*(4-len(leftside)))+leftside
	final_value = leftside+'.'+rightside
	return final_value

def triangulate_surface(name,field,dr,dtheta):
	data = np.load(name+"."+field+".avgheight.npy")
	row,col = data.shape
	
	#wrap in theta direction
	wrapdata = np.zeros((row,col+2))
	wrapdata[:,1:col+1] = data 
	wrapdata[:,0] = wrapdata[:,col]
	wrapdata[:,col+1] = wrapdata[:,1]
	
	wrow,wcol = wrapdata.shape

	for rbin in range(wrow):
		for thetabin in range(wcol):
			if str(wrapdata[rbin][thetabin] != "nan"):
				i1x = (dr*rbin + .5*dr)*(np.cos(thetabin*dtheta + 0.5*dtheta))
				i1y = (dr*rbin + .5*dr)*(np.sin(thetabin*dtheta + 0.5*dtheta))


def Make_surface_PDB(name,field,dr,dtheta):
	data = np.load(name+"."+field+".avgheight.npy")
	serial = 1
	resseqnum = 1
	atom_name = 'SURF'
	resname = 'SURF'
	chain = 'X'
	row,col = data.shape

	with open(name+'.'+field+'.avgheight.pdb', 'w') as f:
		print('CRYST1  150.000  150.000  110.000  90.00  90.00  90.00 P 1           1', file=f)
		for rbin in range(row):
			for thetabin in range(col):
				if str(data[rbin][thetabin]) != "nan":
					seriallen = 5-(len(str(serial)))
					resseqlen = seriallen - 1
					x = (dr*rbin + .5*dr)*(np.cos(thetabin*dtheta + 0.5*dtheta))
					x = coord_format(x)
					y = (dr*rbin + .5*dr)*(np.sin(thetabin*dtheta + 0.5*dtheta))
					y = coord_format(y)
					z = coord_format(data[rbin][thetabin])
					print('HETATM'+(' '*seriallen)+str(serial)+' '+atom_name+' '+resname+chain+(' '*resseqlen)+str(resseqnum)+'    '+x+y+z+'  3.00  0.00              ',file=f)
					serial += 1
					resseqnum +=1
		print('END', file=f)

for system in sys_list:
	for field in field_list:
		triangulate_surface(system,field,6,dth)