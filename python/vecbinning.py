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
    theta = np.arctan2(yval,xval)
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
    else: 
        ValueError("incorrect input for specifier put 'deg_to_rad' or 'rad_to_deg'")
        
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
        self.unbinneddata = np.array(binningdata)
        self.N_bins = N_bins
        
        # Calculate Bins
        self.bin_interval = None
        self.bins = None
        if calc_theta == "theta":
            self.calculate_theta_bins()
        else:
            self.calculate_bins()
        
        # Update Binning
        self.bin_indices = None
        self.update_bin_indices()
        
    def calculate_bins(self):
        """
        Calulates the approiate bin spacing for 
        the specified number of bins 

        Returns
        -------
        None.

        """
        max_bin_interval = np.ceil(max(self.unbinneddata))
        min_bin_interval = np.floor(min(self.unbinneddata))
        self.bin_interval = [min_bin_interval, max_bin_interval]
        self.bins = np.linspace(min_bin_interval, max_bin_interval, self.N_bins)
        
    def calculate_theta_bins(self):
        self.bin_interval = [0, 360]
        self.bins = np.linspace(0, 360, self.N_bins)
        
        
    def update_bin_indices(self):
        """
        Determines which bin each of 
        the values falls into and assigns 
        it a bin index 

        Returns
        -------
        None.

        """
        self.bin_indices = np.digitize(self.unbinneddata, self.bins)


class XYZBinner():
    """
    Class for doing binning soley along the XYZ axis 
    """
    def __init__(self, multi_bin_data:[[list], [list], [list]], N_bins:[list], theta:[list]=None):
        if len(multi_bin_data) != 3:
            raise ValueError("Data must only have an X,Y, and Z components")
        self.multi_bin_data = multi_bin_data
        self.N_bins = N_bins
        if theta == None:
            self.theta = [None] * len(self.multi_bin_data)
        else:
            self.theta = theta
        self.xbin_indices = None
        self.ybin_indices = None
        self.zbin_indices = None
        
        self.xbins = None
        self.ybins = None
        self.zbins = None
        
    def do_xyz_binning(self):
        """
        Does the binning for multi arrays, 
        the data is thought of as x y z and 
        class is the container

        Returns
        -------
        None.

        """
        for i, pre_binned_data in enumerate(self.multi_bin_data):
            if i == 0:
                vec_bin = VectorBinning(pre_binned_data, self.N_bins[i], calc_theta = self.theta[i])
                self.xbin_indices = vec_bin.bin_indices
                self.xbins = vec_bin.bins
            if i == 1:
                vec_bin = VectorBinning(pre_binned_data, self.N_bins[i], calc_theta = self.theta[i])
                self.ybin_indices = vec_bin.bin_indices
                self.ybins = vec_bin.bins
            if i == 2:
                vec_bin = VectorBinning(pre_binned_data, self.N_bins[i], calc_theta = self.theta[i])
                self.zbin_indices = vec_bin.bin_indices
                self.zbins = vec_bin.bins
            
            