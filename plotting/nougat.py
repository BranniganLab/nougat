"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import warnings
from height import *
from curvature import *
from quantity_analysis import *
from tilt import *
from utils import *


class Membrane:
    """
    The basic class for nougat. A membrane is comprised of multiple Fields and\
    Field_sets. Each field contains a 3D numpy ndarray that corresponds to some\
    measurement of interest (E.G. the height of the outer leaflet, or the mean\
    curvature of the bilayer midplane) over the course of an MD trajectory.

    Attributes
    ----------
    active_list  :  list
        The list of all Fields and Field_sets that have been computed. This\
        list is updated any time Fields are turned into Field_sets so that\
        there is no duplication.
    todo_list  :  list
        The list of quantities that the user has selected for analysis.
    composition  :  str
        The composition of the membrane.
    t0  :  float
        The equilibrium thickness of the membrane.

    List of all default Fields and Field_sets
    ----------
    height  :  Field_set
        The height (in z) of the outer and inner leaflets, as well as the\
        symmetric/antisymmetric variables "z plus" and "z minus". Z plus is the\
        bilayer midplane. Z minus is the bilayer thickness. Heights may be\
        relative to some reference point on a protein or could be the absolute\
        height of the membrane in VMD.
    zzero  :  Field
        The height (in z) of the interface between the outer and inner leaflet.\
        "Z zero" is commonly thought to be equal to the bilayer midplane but\
        this is not always the case, especially around inclusions.
    """

    def __init__(self, polar, todo_list, composition=None, t0=None):
        """
        Create a Membrane object.

        Parameters
        ----------
        polar  :  bool
            A switch for using cylindrical versus Cartesian coordinates.
        todo_list  :  list
            The list of quantities that the user has selected for analysis.
        """
        self.active_list = []
        self.polar = polar
        self.todo_list = todo_list
        self.composition = composition
        self.t0 = t0

    def __iter__(self):
        """Iterate through active_list."""
        for item in self.active_list:
            if isinstance(item, Field):
                yield item
            elif isinstance(item, Field_set):
                for field in item:
                    yield field

    def create_Field_set(self, outer, inner, name):
        """
        Create a Field_set object by supplying the inner and outer leaflet\
        Fields. These will be incorporated into a Field_set that then calculates\
        the symmetric "plus" and anti-symmetric "minus" fields, respectively.

        Parameters
        ----------
        outer : Field
            The Field object that corresponds to the outer leaflet quantity.
        inner : Field
            The Field object that corresponds to the inner leaflet quantity.

        Returns
        -------
        new_Field_set  :  Field_set
            The new Field_set object you just created.

        """
        new_Field_set = Field_set(self.polar, outer, inner, name)
        self.active_list.remove(outer)
        self.active_list.remove(inner)
        self.active_list.append(new_Field_set)
        return new_Field_set

    def create_Field(self, path, name, quantity=None, leaflet=None):
        """
        Create a Field object. Use this method before attempting to create a\
        Field_set.

        Parameters
        ----------
        path  :  Path, str, or ndarray
            Either contains a path to TCL output data that needs to be parsed,\
            or contains a numpy ndarray that should just be saved into the\
            Field's field_data attribute.
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height', 'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo', or 'zzero'.

        Returns
        -------
        new_Field  :  Field
            The new Field object you just created.

        """
        new_Field = Field(path, self.polar, name, quantity, leaflet)
        self.active_list.append(new_Field)
        return new_Field

    def create_Vector_field(self):
        """Not implemented yet."""
        return NotImplemented


class Field:
    """A field contains a measurement of some surface over the course of an MD\
    trajectory. This could be the height of the outer leaflet, the mean curvature\
    of the bilayer midplane, etc...

    Attributes
    ----------
    polar  :  bool
        A switch for using cylindrical versus Cartesian coordinates.
    field_data  :  ndarray
        A three-dimensional array containing data. Could be height values,\
        curvature values, etc. Assume zero-th dimension to be time (frames\
        in trajectory), 1st and 2nd dims to be x/r and y/theta bins. Individual\
        cells can contain int or float.
    grid_dims  :  tuple
        The shape of the field_data ndarray.
    avg  :  ndarray
        The 2D ndarray that represents the average over time.
    avg_over_theta : ndarray
        If polar coordinates were used, this is the 1D ndarray that represents\
        the average over time in each radial bin.
    """

    def __init__(self, path, polar, name, quantity=None, leaflet=None):
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
        name  :  str
            The name of this Field
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height', 'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo', or 'zzero'.

        """
        self.polar = polar
        self.name = name

        # read in the data
        if isinstance(path, (Path, str)):
            assert quantity is not None, "quantity is required in order to use a path"
            assert leaflet is not None, "leaflet name is required in order to use a path"
            self.field_data, self.grid_dims = parse_tcl_output(path, quantity, leaflet)
        elif isinstance(path, np.ndarray):
            self.field_data = path
            self.grid_dims = np.shape(path)
        else:
            raise ValueError("path must either be a numpy ndarray or a path")

        assert len(self.grid_dims) == 3, "This Field should contain a 3D array"
        if self.polar is False:
            # remove this if you ever allow for rectangular boxes in Cart. coords.
            assert self.grid_dims[0] == self.grid_dims[1], "Your box is not square"

        # calculate averages if appropriate
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.avg = np.nanmean(self.field_data, axis=0)
            if self.polar:
                self.avg_over_theta = np.nanmean(self.avg, axis=0)

    def __iter__(self):
        """Return self so that Membrane.active_list can be looped easily."""
        return self

    def __str__(self):
        """Say your name, rather than your address."""
        return self.name

    # BASIC MATH MAGIC METHODS BELOW #

    def __add__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, Field):
            return self.field_data + other.field_data
        elif isinstance(other, (np.ndarray, int, float)):
            return self.field_data + other
        else:
            return NotImplemented

    def __radd__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, (np.ndarray, int, float)):
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

    def __mul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, Field):
            return self.field_data * other.field_data
        elif isinstance(other, (np.ndarray, int, float)):
            return self.field_data * other
        else:
            return NotImplemented

    def __rmul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.field_data * other
        else:
            return NotImplemented

    def __div__(self, other):
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


class Field_set:
    def __init__(self, polar, outer, inner, name):
        self.outer = outer
        self.inner = inner
        self.name = name
        self.plus = Field((outer + inner) / 2., polar, self.name + "_plus")
        self.minus = Field((outer - inner) / 2., polar, self.name + "_minus")

    def __iter__(self):
        """Iterate through the four Fields in a Field_set."""
        for f in [self.outer, self.inner, self.plus, self.minus]:
            yield f

    def __str__(self):
        """Say your name, rather than your address."""
        return self.name


class Vector_field(Field):
    """Not implemented yet."""
    pass


'''
def run_nougat(polar, inclusion_drawn):
    """
Run nougat's averaging and image processing routines.

 Parameters
  ----------
   polar: boolean
     True for polar coordinate system
      False for cartesian
    inclusion_drawn: boolean
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
'''

m = Membrane(True, None, "100% POPC", 3.5)
test1 = np.ones((10, 5, 5))
test2 = test1 * 2
outer = m.create_Field(test2, "z_one")
inner = m.create_Field(test1, "z_two")
height = m.create_Field_set(outer, inner, "z")
for f in m:
    print(f)
print(m.active_list)
