#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 16:58:32 2025

@author: jje63
"""

import numpy as np

def convert_cartesian_to_polar(xvalues:[list],yvalues:[list]):
    """
    Converts cartesian x and y coordinates to polar coordinates

    Parameters
    ----------
    xvalues : [list]
        Values associated with the x component.
    yvalues : [list]
       Values associated with the y component.

    Returns
    -------
    list
        Returns list of r and theta values associted with the x and y.

    """
    xval = np.array(xvalues)
    yval = np.array(yvalues)
    
    r = np.sqrt(xval**2 + yval**2)
    theta = np.arctan(np.divide(yval,xval))
    return [r, theta]

def convert_radians_and_degrees(values:[list], specifier:str="rad_to_deg"):
    """
    Converts radians to degrees or vice versa 

    Parameters
    ----------
    values : [list]
        Values to convert.
    specifier : str, optional
        Specify if you want to go from radians to degrees 
        or degrees to radians. The default is "rad_to_deg".

    Returns
    -------
    list
        returns the list of numbers converted to radians or degrees.

    """
    if specifier == "deg_to_rad":
        return np.array(values)*(np.pi/180)
    elif specifier == "rad_to_deg":
        return np.array(values)*(180/np.pi)
        
class VectorBinning:
    """
    Class for doing binning along a single axis 
    """
    def __init__(self, binningdata:[list], N_bins:int, calc_theta:str=None):
        """
        Data for initializing class attributes

        Parameters
        ----------
        binningdata : [list]
            list of floats of ints to bin.
        N_bins : int
            Number of bins to be created.

        Returns
        -------
        None.

        """
        self.prebinneddata = np.array(binningdata)
        self.N_bins = N_bins
        
        # Calculate Bins
        self.bin_interval = None
        self.bins = None
        if calc_theta == "theta":
            self.calculate_theta_bins()
        else:
            self.calculate_bins()
        
        # Update Binning
        self.bin_indicies = None
        self.update_bin_indicies()
        
    def calculate_bins(self):
        """
        Calulates the approiate bin spacing for 
        the specified number of bins 

        Returns
        -------
        None.

        """
        max_bin_interval = np.ceil(max(self.prebinneddata))
        min_bin_interval = np.floor(min(self.prebinneddata))
        self.bin_interval = [min_bin_interval, max_bin_interval]
        self.bins = np.linspace(min_bin_interval, max_bin_interval, self.N_bins)
        
    def calculate_theta_bins(self):
        self.bin_interval = [0, 360]
        self.bins = np.linspace(0, 360, self.N_bins)
        
        
    def update_bin_indicies(self):
        """
        Determines which bin each of 
        the values falls into and assigns 
        it a bin index 

        Returns
        -------
        None.

        """
        self.bin_indicies = np.digitize(self.prebinneddata, self.bins)

class Multibinner:
    """
    Class for doing binning along multiple axis 
    """
    def __init__(self, multi_bin_data:[list [list]], N_bins:[list], theta:[list]=None):
        """
        Data for initializing class attributes

        Parameters
        ----------
        multi_bin_data : [list [list]]
            list of list of floats of ints to bin.
        N_bins : [list]
            list of the number of bins to be created.

        Returns
        -------
        None.

        """
        self.multi_bin_data = multi_bin_data
        self.N_bins = N_bins
        if theta == None:
            self.theta = [None] * len(self.multi_bin_data)
        else:
            self.theta = theta
        self.raw_data = []
        self.raw_bin_indicies = []
        self.raw_binning_data = []
        self.multi_bin_dictionary = {}

    def do_multi_binning(self):
        """
        Does the binning for multi arrays, 
        the data is thought of as multiple coordinates to bin 

        Returns
        -------
        None.

        """
        for i, data_array in enumerate(self.multi_bin_data):
            vec_bin = VectorBinning(data_array, self.N_bins[i], calc_theta = self.theta[i])
            self.raw_data.append(vec_bin)
            self.raw_bin_indicies.append(vec_bin.bin_indicies)
            self.raw_binning_data.append(vec_bin.prebinneddata)
        self.raw_binning_data = np.stack(self.raw_binning_data)
        for i, bindata in enumerate(zip(*self.raw_bin_indicies)):
            if bindata in self.multi_bin_dictionary:
                self.multi_bin_dictionary[bindata].append(self.raw_binning_data[:,i])
            else:
                self.multi_bin_dictionary[bindata] = [self.raw_binning_data[:,i]]
            
data = [[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15],[1, 5, 2, 8, 12, 6, 15]]

#l = Multibinner(data, [5,10,6])
#l.do_multi_binning()
data[0],data[1]=convert_cartesian_to_polar(data[0], convert_radians_and_degrees(data[1]))

l = Multibinner(data, [5,10,6], theta=[None,"theta",None])
l.do_multi_binning()
print(l.multi_bin_dictionary)
print(l.raw_data[1].bins)
        
    