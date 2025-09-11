#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:42:51 2025

@author: jje63
"""

import numpy as np
import argparse
import os
from pathlib import Path

def read_nougat_data(filepath,filename, seclen):
    """
    Creates a numpy array from nougat txt document 

    Parameters
    ----------
    filepath : str
        Path for directory of the file.
    filename : str
        Name of the file.
    seclen : int
        Number of output variables in each frame.

    Returns
    -------
    systemArray : arr
        Nougat text document converted to a numpy array.

    """
    fullpath = Path(filepath+"/"+filename)
    loadFile = np.genfromtxt(fullpath,dtype='str')
    systemArray = loadFile.reshape(int(loadFile.shape[0]/seclen), seclen, int(loadFile.shape[1]))
    return systemArray

def save_nougat_data(filepath, nparry, boole):
    """
    Can save the numpy array as a text file

    Parameters
    ----------
    filepath : str
        Path for directory of the file.
    nparry : arr
        numpy array.
    boole : bool
        If true then save.

    Returns
    -------
    None.

    """
    if boole:
        np.save(Path(filepath+"/"+"NGArrayData.npy"), nparry)
    

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, default=os.getcwd(), help="Path to nougat output file directory")
    parser.add_argument("-n", "--name", type=str, default="full_file.dat", help="Name of file")
    parser.add_argument("-l", "--seclen", type=int, default=6, help="Number of output varibles in each frame")
    parser.add_argument("-s", "--save", type=bool, default=False, help="Save Frame")
    args = parser.parse_args()
    ngArray = read_nougat_data(args.path, args.name, args.seclen)
    save_nougat_data(args.path, ngArray, args.save)
    