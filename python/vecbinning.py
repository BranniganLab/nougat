#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 16:58:32 2025

@author: jje63
"""
from pathlib import Path

import numpy as np
import itertools

class VectorBinning:
    def __init__(self, binningdata:[list], N_bins:int):
        self.prebinneddata = np.array(binningdata)
        self.N_bins = N_bins
        
        # Calculate Bins
        self.bin_interval = None
        self.bins = None
        self.calculate_bins()
        
        # Update Binning
        self.bin_indicies = None
        self.update_bin_indicies()
        
    def calculate_bins(self):
        max_bin_interval = np.ceil(max(self.prebinneddata))
        min_bin_interval = np.floor(min(self.prebinneddata))
        self.bin_interval = [min_bin_interval, max_bin_interval]
        self.bins = np.linspace(min_bin_interval, max_bin_interval, self.N_bins)
        
    def update_bin_indicies(self):
        self.bin_indicies = np.digitize(self.prebinneddata, self.bins)

class Multibinner:
    def __init__(self, multi_bin_data:[list [list]], N_bins:[list]):
        self.multi_bin_data = multi_bin_data
        self.N_bins = N_bins
        self.raw_data = []
        self.raw_bin_indicies = []
        self.raw_binning_data = []
        self.multi_bin_dictionary = {}

    def do_multi_binning(self):
        for i, data_array in enumerate(self.multi_bin_data):
            self.raw_data.append(VectorBinning(data_array, self.N_bins[i]))
            self.raw_bin_indicies.append(VectorBinning(data_array, self.N_bins[i]).bin_indicies)
            self.raw_binning_data.append(VectorBinning(data_array, self.N_bins[i]).prebinneddata)
        self.raw_binning_data = np.stack(self.raw_binning_data)
        for i, bindata in enumerate(zip(*self.raw_bin_indicies)):
            if bindata in self.multi_bin_dictionary:
                self.multi_bin_dictionary[bindata].append(self.raw_binning_data[:,i])
            else:
                self.multi_bin_dictionary[bindata] = [self.raw_binning_data[:,i]]
            
data = [[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15]]

l = Multibinner(data, [5,10,6])
l.do_multi_binning()
print(l.raw_data[0].bins)
        
    