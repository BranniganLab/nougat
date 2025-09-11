#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 12:16:49 2025

@author: jje63
"""
import subprocess
from pathlib import Path
import numpy as np


class parseTCL:
    def __init__(self, data_path, row=6, column=5216,frames=2):
        self.dataset = np.genfromtxt(data_path, comments="#",delimiter=(8,)*column)
        #print(len((8,)*column))
        zstack = int(self.dataset.size/row/column)
        #print(zstack)
        self.frames = self.dataset.reshape(zstack,row,column)
data1 = Path("/home/jje63/Documents/SingleNP_Sim/testnougatfolder/full_file.dat")
pdata = parseTCL(data1)
print(pdata.frames[1][5][5000])
print(pdata.frames.shape)
