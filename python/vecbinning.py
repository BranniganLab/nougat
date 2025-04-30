#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 16:58:32 2025

@author: jje63
"""
from pathlib import Path

import numpy as np

class VecBinning:
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
        step_size = (max_bin_interval - min_bin_interval)/self.N_bins
        self.bin_interval = [min_bin_interval, max_bin_interval]
        self.bins = np.arange(min_bin_interval, max_bin_interval+step_size, step_size, dtype=int)
        
    def update_bin_indicies(self):
        print(self.bins)
        self.bin_indicies = np.digitize(self.prebinneddata, self.bins)
        
data = [1, 5, 2, 8, 12, 6, 15]     
data = VecBinning(data, 4)

        
    