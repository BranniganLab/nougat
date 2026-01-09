#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:42:51 2025

@author: jje63
"""

import numpy as np
import argparse
from pathlib import Path

def read_nougat_data(filepath, seclen):
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
        Nougat text document converted to a Ndimensional numpy array.

    """
    
    loadFile = np.genfromtxt(filepath, dtype='str')
    array_shape = loadFile.shape[0]/seclen
    assert isinstance(array_shape, int), "Section length gave non-integer array shape"
    systemArray = loadFile.reshape(int(loadFile.shape[0]/seclen), seclen, int(loadFile.shape[1]))
    return systemArray
    

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-ip", "--inputpath",  default=Path.cwd().joinpath("full_file.dat").resolve(), help="Path to nougat input file")
    parser.add_argument("-op", "--outputpath",  default=Path.cwd().joinpath("NGArrayData.npy"), help="Path to nougat output file")
    parser.add_argument("-l", "--seclen", type=int, default=6, help="Number of output varibles in each frame")
    args = parser.parse_args()
    if isinstance(args.inputpath, str):
        args.inputpath = Path(args.inputpath).resolve()
    if isinstance(args.outputpath, str):
        args.inputpath = Path(args.outputpath).resolve()
    ngArray = read_nougat_data(args.inputpath, args.seclen)
    np.save(args.outputpath, ngArray)
    