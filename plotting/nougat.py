"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from height import *
from curvature import *
from quantity_analysis import *
from tilt import *
from utils import *


class Field(Membrane):
    """A field contains x y and z.

    Attributes
    ----------
    polar  :  bool
        A switch for using cylindrical versus Cartesian coordinates.
    field_data  :  ndarray
        A two- or three-dimensional array containing data. Could be height\
        values, curvature values, etc. If three-dimensional, assume zero-th\
        dimension to be time (frames in trajectory).
    grid_dims  :  tuple
        The shape of the field_data ndarray.
    avg  :  ndarray
        If this Field is a trajectory, then avg is the the 2D ndarray that\
        represents the average over time.
    avg_over_theta : ndarray
        If polar coordinates were used, this is the 1D ndarray that represents\
        the average over time in each radial bin.
    """

    def __init__(self, path, polar, quantity=None, leaflet=None):
        """
        Construct a Field object.

        Parameters
        ----------
        path  :  Path, str, or ndarray
            Either contains a path to TCL output data that needs to be parsed,\
            or contains a numpy ndarray that should just be saved into the\
            Field's field_data attribute.
        polar  :  bool
            If true, use cylindrical coordinates. Otherwise, use Cartesian.
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height', 'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo', or 'zzero'.

        """
        self.polar = polar

        # read in the data
        if isinstance(path, (Path, str)):
            assert quantity is not None
            assert leaflet is not None
            self.field_data, self.grid_dims = parse_tcl_output(path, quantity, leaflet)
        elif isinstance(path, np.ndarray):
            self.field_data = path
            self.grid_dims = np.shape(self.field_data)
        else:
            raise ValueError("path must either be a numpy ndarray or a path")

        assert len(self.grid_dims) >= 2

        # calculate averages if appropriate
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if len(self.grid_dims) == 3:
                self.avg = np.nanmean(self.field_data, axis=0)
                if self.polar:
                    self.avg_over_theta = np.nanmean(self.avg, axis=0)
            elif (len(self.grid_dims) == 2) and self.polar:
                self.avg_over_theta = np.nanmean(self.field_data, axis=0)

    def __add__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, Field):
            return self.field_data + other.field_data
        elif isinstance(other, (np.ndarray, int, float)):
            return self.field_data + other
        else:
            return NotImplemented

    def __sub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, Field):
            return self.field_data - other.field_data
        elif isinstance(other, (np.ndarray, int, float)):
            return self.field_data - other
        else:
            return NotImplemented

    def __mul__(self,other):
        """Use numpy to multiply things."""
        if isinstance(other, Field):
            return self.field_data * other.field_data
        elif isinstance(other, (np.ndarray, int, float)):
            return self.field_data * other
        else:
            return NotImplemented

    def __div__(self,other):
        """Use numpy to divide things."""
        if isinstance(other, Field):
            return self.field_data / other.field_data
        elif isinstance(other, (np.ndarray, int, float)):
            return self.field_data / other
        else:
            return NotImplemented

    def __pow__(self, exponent):
        """Use numpy's power() on the array stored in this Field."""
        return np.power(self.field_data, exponent)


class Field_set(Membrane):
    def __init__(self, path, polar, quantity):
        if quantity is "height":
            self.outer = field(path, polar, quantity, "zone")
            self.inner = field(path, polar, quantity, "ztwo")
            self.plus = field(self.outer+self.inner, polar)
            self.minus = field(self.outer-self.inner, polar)
        elif quantity is "curvature":
            


class Membrane:
    def __init__(self, path, polar, list_of_options):
        if "height" in list_of_options:
            self.height = Field_set(path, polar, "height")
            self.zzero = Field(path, polar, "height", "zzero")
            if "curvature" in list_of_options:
                mean_curvature = Field_set()


def run_nougat(polar, inclusion_drawn):
    """
    Run nougat's averaging and image processing routines.

    Parameters
    ----------
    polar : boolean
        True for polar coordinate system; False for cartesian
    inclusion_drawn : boolean
        This feature currently isn't implemented, but would be how you \
        include an inclusion in your graphics

    Returns
    -------
    None.

    """
    cwd = Path.cwd()

    # make necessary folders
    create_outfile_directories(cwd)

    # define inclusion if present
    if inclusion_drawn is True:
        inclusion = add_inclusion(name, field_list)  # this proc doesn't exist yet!
    else:
        inclusion = False

    # figure out all the important info about the system you'll need
    system_dict = read_log()

    # read in height files and calculate surface heights
    # parse_height returns system_dict bc it adds nframes to the dictionary
    system_dict = parse_height(system_dict, polar, cwd)

    calculate_thickness(polar, cwd)

    calculate_curvature(polar, system_dict, cwd)

    # $$$$$$$$$$ UNTESTED FEATURES BELOW $$$$$$$$$$$

    # save_areas(system_dict["bin_info"], 0, polar)

    # calculate_density(polar, system_dict, cwd)

    # calculate_order(polar, system_dict, cwd)

    # calculate_tilt(sys_name, system_dict, coordsys, inclusion, cwd)

    # calc_elastic_terms(".", coordsys, config_dict, system_dict['bin_info'])

    plot_all_quantities(polar, system_dict, cwd, inclusion)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl in polar coordinates")
    # parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
    args = parser.parse_args()

    run_nougat(args.polar, False)

    print("Thank you for using nougat!")
