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
    grid_dims  :  dict
        A dictionary that contains information about the grid dimensions used.\
        N1 is the number of bins in the first dimension (x/r) and d1 is the\
        distance between bin centers. N2 and d2 are similar, but along the\
        second dimension (y/theta). Nframes is the number of frames in the\
        trajectory.
    composition  :  str
        The composition of the membrane.
    t0  :  float
        The equilibrium thickness of the membrane.
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
        composition  :  str
            The composition of the membrane.
        t0  :  float
            The equilibrium thickness of the membrane.
        """
        self.active_list = []
        self.polar = polar
        self.todo_list = todo_list
        self.composition = composition
        self.t0 = t0
        self.grid_dims = {
            "N1": None,
            "N2": None,
            "Nframes": None,
            "d1": None,
            "d2": None
        }

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
        name  :  str
            The name you want to give this Field_set.

        Returns
        -------
        new_Field_set  :  Field_set
            The new Field_set object you just created.

        """
        new_Field_set = Field_set(outer, inner, name, self)
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
        name  :  str
            The name you want to give this Field.
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height', 'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo', or 'zzero'.

        Returns
        -------
        new_Field  :  Field
            The new Field object you just created.

        """
        new_Field = Field(path, name, self, quantity, leaflet)
        self.active_list.append(new_Field)
        return new_Field

    def create_Vector_field(self):
        """Not implemented yet."""
        return NotImplemented

    def measure_correlation(self, field1, field2):
        """
        Measure the correlation between two Fields.

        Parameters
        ----------
        field1 : Field
            One of the two Fields.
        field2 : Field
            The other Field.

        Returns
        -------
        corr  :  Field
            The correlation between Fields 1 and 2.
        """
        together = calc_avg_over_time(field1 * field2)
        apart = field1.avg * field2.avg
        corr = together - apart
        return self.create_Field(corr, "corr_" + field1.name + "_" + field2.name)

    def measure_rms_of_field(self, field, eq2=None):
        """
        Measure the root-mean-squared value of some Field. Optionally, supply\
        a squared equilibrium value in order to generate a fold-enrichment score.

        Parameters
        ----------
        field : Field
            The Field you want analyzed.
        eq2 : float, optional
            The squared equilibrium value of the Field. The default is None.

        Returns
        -------
        Field
            The rms or enrichment (if eq2 used) values for the Field.
        """
        squared = field**2
        if eq2 is not None:
            squared = squared / eq2
        mean_squared = calc_avg_over_time(squared)
        rms = np.sqrt(mean_squared)
        if eq2 is not None:
            return self.create_Field(rms, "rmsTilde_" + field.name)
        else:
            return self.create_Field(rms, "rms_" + field.name)


class Field:
    """A field contains a measurement of some surface over the course of an MD\
    trajectory. This could be the height of the outer leaflet, the mean curvature\
    of the bilayer midplane, etc...

    Attributes
    ----------
    field_data  :  ndarray
        A two- or three-dimensional array containing data. Could be height values,\
        curvature values, etc. Assume zero-th dimension to be time (frames\
        in trajectory), 1st and 2nd dims to be x/r and y/theta bins. Individual\
        cells can contain int or float.
    avg  :  ndarray
        The 2D ndarray that represents the average over time.
    avg_over_theta : ndarray
        If polar coordinates were used, this is the 1D ndarray that represents\
        the average over time in each radial bin.
    parent  :  Membrane
        The Membrane object to which this field belongs.
    """

    def __init__(self, path, name, parent, quantity=None, leaflet=None):
        """
        Construct a Field object.

        Parameters
        ----------
        path  :  Path, str, or ndarray
            Either contains a path to TCL output data that needs to be parsed,\
            or contains a numpy ndarray that should just be saved into the\
            Field's field_data attribute.
        name  :  str
            The name of this Field
        parent  :  Membrane
            The Membrane object that his Field belongs to.
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height', 'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo', or 'zzero'.

        """
        self.parent = parent
        self.name = name

        # read in the data
        if isinstance(path, (Path, str)):
            assert quantity is not None, "quantity is required in order to use a path"
            assert leaflet is not None, "leaflet name is required in order to use a path"
            self.field_data = self.parse_tcl_output(path, quantity, leaflet)
        elif isinstance(path, np.ndarray):
            self.field_data = path
            if len(np.shape(self.field_data)) == 2:
                if self.parent.grid_dims["N1"] is not None:
                    err_msg = "This ndarray doesn't have the same dimensions as its parent Membrane."
                    assert self.parent.grid_dims["N1"] == np.shape(path)[0], err_msg
                    assert self.parent.grid_dims["N2"] == np.shape(path)[1], err_msg
                else:
                    self.parent.grid_dims["N1"] = np.shape(path)[0]
                    self.parent.grid_dims["N2"] = np.shape(path)[1]
            elif len(np.shape(self.field_data)) == 3:
                if self.parent.grid_dims["N1"] is not None:
                    err_msg = "This ndarray doesn't have the same dimensions as its parent Membrane."
                    assert self.parent.grid_dims["N1"] == np.shape(path)[1], err_msg
                    assert self.parent.grid_dims["N2"] == np.shape(path)[2], err_msg
                else:
                    self.parent.grid_dims["N1"] = np.shape(path)[1]
                    self.parent.grid_dims["N2"] = np.shape(path)[2]
                    self.parent.grid_dims["Nframes"] = np.shape(path)[0]
            else:
                raise Exception("Field must contain a 2- or 3D array.")
        else:
            raise ValueError("path must either be a numpy ndarray or a path")

        assert len(np.shape(self.field_data)) == 3, "This Field should contain a 3D array"

        # calculate averages if appropriate
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if self.parent.grid_dims["Nframes"] > 1:
                self.avg = np.nanmean(self.field_data, axis=0)
            elif self.parent.grid_dims["Nframes"] == 1:
                self.avg = self.field_data
            else:
                raise Exception("You have created a Field with no data in it somehow")
            if self.parent.polar:
                self.avg_over_theta = np.nanmean(self.avg, axis=1)

    def __iter__(self):
        """Return self so that Membrane.active_list can be looped easily."""
        return self

    def __str__(self):
        """Say your name, rather than your address."""
        return self.name

    def __repr__(self):
        """Say your name, rather than your address."""
        return self.name

    def parse_tcl_output(self, path, quantity, leaflet):
        """
        Read in the tcl output data, update the parent Membrane's grid_dims,\
        and generate the field_data array.

        Parameters
        ----------
        path : Path or str
            The path to the nougat.tcl results folder.
        quantity : str
            A valid nougat.tcl output quantity, E.G. "height", "order", etc.
        leaflet : str
            A valid nougat.tcl leaflet name, I.E. "zone", "ztwo", or "zzero".

        Returns
        -------
        field_data : ndarray
            A 3D array containing nougat.tcl output data. The zero-th dimension\
            is time (frames), the first dimension is x or r, and the second\
            dimension is y or theta.
        """
        # import traj values
        input_file_path = path.joinpath("tcl_output", quantity, leaflet + ".dat")
        unrolled_data = np.genfromtxt(input_file_path, missing_values='nan', filling_values=np.nan)

        err_msg = "This ndarray doesn't have the same dimensions as its parent Membrane."

        # determine grid_dims along second dimension
        N2 = np.shape(unrolled_data)[1] - 2
        d2 = (2 * np.pi) / N2
        if self.parent.grid_dims["N2"] is not None:
            assert self.parent.grid_dims["N2"] == N2, err_msg
            assert self.parent.grid_dims["d2"] == d2, err_msg
        else:
            self.parent.grid_dims["N2"] = N2
            self.parent.grid_dims["d2"] = d2

        # determine grid_dims along first dimension
        d1 = unrolled_data[0, 1] - unrolled_data[0, 0]
        """nougat.tcl's output is structured such that the starting value of\
        x or r will be repeated each time there is a new frame. Look for the\
        first repeat and you will know how many x/r bins there are"""
        match_value = unrolled_data[0, 0]
        index = np.where(unrolled_data[1:, 0] == match_value)
        if index[0].size != 0:
            N1 = index[0][0] + 1
        else:
            # this would happen if there was only one frame in the trajectory
            N1 = np.shape(unrolled_data[0])
        assert np.shape(unrolled_data)[0] % N1 == 0, "N1 incorrectly calculated, or error in nougat.tcl write-out stage."
        if self.parent.grid_dims["N1"] is not None:
            assert self.parent.grid_dims["N1"] == N1, err_msg
            assert self.parent.grid_dims["d1"] == d1, err_msg
        else:
            self.parent.grid_dims["N1"] = N1
            self.parent.grid_dims["d1"] = d1

        # determine Nframes
        Nframes = int(np.shape(unrolled_data)[0] / N1)
        if self.parent.grid_dims["Nframes"] is not None:
            assert self.parent.grid_dims["Nframes"] == Nframes, err_msg
        else:
            self.parent.grid_dims["Nframes"] = Nframes

        # create a new array that has each frame in a different array level
        field_data = np.zeros((Nframes, N1, N2))
        for frm in range(Nframes):
            field_data[frm, :, :] = unrolled_data[frm * N1: (frm + 1) * N1, 2:]

        return field_data

    # BASIC MATH MAGIC METHODS BELOW #
    # These make it so that you can do math on the field object, rather than\
    # having to specify the object's .field_data attribute every time.

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
    """A Field_set contains four Fields: the outer and inner leaflets, plus the\
    symmetric and anti-symmetric variables <Field>_plus and <Field>_minus. This\
    object serves as a constructor for the _plus and _minus variables.

    Attributes
    ----------
    outer  :  Field
        The outer leaflet Field.
    inner  :  Field
        The inner leaflet Field.
    plus  :  Field
        Outer leaflet + inner leaflet, divided by 2.
    minus  :  Field
        Outer leaflet - inner leaflet, divided by 2.
    name  :  str
        The name of this Field set, and the prefix that will be given to the\
        _plus and _minus fields.
    parent  :  Membrane
        The Membrane object to which the Fields belong.
    """

    def __init__(self, outer, inner, name, parent):
        """
        Construct a Field_set.

        Parameters
        ----------
        outer  :  Field
            The outer leaflet Field.
        inner  :  Field
            The inner leaflet Field.
        name  :  str
            The name of this Field set, and the prefix that will be given to the\
            _plus and _minus fields.
        parent  :  Membrane
            The Membrane object to which the Fields belong.

        Returns
        -------
        None.
        """
        self.outer = outer
        self.inner = inner
        self.name = name
        self.parent = parent
        self.plus = Field((outer + inner) / 2., self.name + "_plus", self.parent)
        self.minus = Field((outer - inner) / 2., self.name + "_minus", self.parent)

    def __iter__(self):
        """Iterate through the four Fields in a Field_set."""
        for f in [self.outer, self.inner, self.plus, self.minus]:
            yield f

    def __str__(self):
        """Say your name, rather than your address."""
        return self.name

    def __repr__(self):
        """Say your name, rather than your address."""
        return self.name


class Vector_field(Field):
    """Not implemented yet."""

    pass


def run_nougat(polar, quantities):
    """
    Run nougat's averaging and image processing routines.

    Parameters
    ----------
    polar: boolean
        True for cylindrical coordinate system, False for Cartesian.
    quantities: str
        A string that specifies which quantities to carry out analysis on. If\
        None, assume all quantities should be analyzed.

    Returns
    -------
    None.

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
    thickness  :  Field_set
        The thickness of the outer and inner leaflets, as well as the symmetric\
        and anti-symmetric variables "t plus" and "t minus." T plus is the\
        average leaflet thickness and t minus is the leaflet thickness asymmetry.
    epsilon  :  Field
        An alternative measurement of leaflet thickness asymmetry, defined in\
        [Watson & Brown, PRL, 2012].
    mean_curvature  :  Field_set
        The mean curvature (H) of the membrane outer and inner leaflets, as well as\
        the symmetric/anti-symmetric variables "H plus" and "H minus".
    gaussian_curvature  :  Field_set
        The Gaussian curvature (K) of the membrane outer and inner leaflets, as well\
        as the symmetric/anti-symmetric variables "K plus" and "K minus".
    """
    cwd = Path.cwd()  # TO-DO: change this so that path is supplied by user!
    todo_list = make_todo_list(quantities)

    m = Membrane(polar, todo_list)
    if "height" in m.todo_list:
        zone = m.create_Field(cwd, "z_one", "height", "zone")
        ztwo = m.create_Field(cwd, "z_two", "height", "ztwo")
        zzero = m.create_Field(cwd, "z_zero", "height", "zzero")
        height = m.create_Field_set(ztwo, zone, "z")
    if "thickness" in m.todo_list:
        tone = m.create_Field(zone - zzero, "t_one")
        ttwo = m.create_Field(zzero - ztwo, "t_two")
        thickness = m.create_Field_set(tone, ttwo, "t")

    print(thickness.minus.field_data[0, 0, 0])
    print(tone.field_data)
    print(m.active_list)
    print(m)
    for f in m:
        print(f)
    '''
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

    calc_elastic_terms(cwd, polar, system_dict['bin_info'])

    plot_all_quantities(polar, system_dict, cwd, inclusion)
    '''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze output from nougat.tcl")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl with polar coordinates")
    parser.add_argument("-q", "--quantities", help="Specify the quantities you want to calculate: height=h, thickness=t, curvature=c, order=o")
    # parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
    args = parser.parse_args()

    run_nougat(args.polar, args.quantities)

    print("Thank you for using nougat!")
