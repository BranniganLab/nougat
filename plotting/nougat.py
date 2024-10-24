"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

from pathlib import Path
import numpy as np
import warnings
from utils import calc_avg_over_time, bin_prep, plot_maker, mostly_empty, read_log


class Membrane:
    """
    The basic class for nougat. A membrane is comprised of multiple Fields and\
    Field_sets. Each field contains a 3D numpy ndarray that corresponds to\
    some measurement of interest (E.G. the height of the outer leaflet, or the\
    mean curvature of the bilayer midplane) over the course of an MD\
    trajectory.

    Attributes
    ----------
    children  :  list
        The list of all Fields and Field_sets that have been computed. This\
        list is updated any time Fields are turned into Field_sets so that\
        there is no duplication.
    grid_dims  :  dict
        Contains information about the grid dimensions used.\
        N1 is the number of bins in the first dimension (x/r) and d1 is the\
        distance between bin centers. N2 and d2 are similar, but along the\
        second dimension (y/theta). Nframes is the number of frames in the\
        trajectory.
    polar  :  bool
        If True, use cylindrical coordinates.
    composition  :  dict
        The composition of the membrane where keys are lipid names and values \
        are composition percentages.
    t0  :  float
        The equilibrium thickness of the membrane.
    """

    def __init__(self, polar):
        """
        Create a Membrane object.

        Parameters
        ----------
        polar  :  bool
            A switch for using cylindrical versus Cartesian coordinates.
        """
        self.children = {}
        self.polar = polar
        self.composition = {}
        self.t0 = None
        self.grid_dims = {
            "N1": None,
            "N2": None,
            "Nframes": None,
            "d1": None,
            "d2": None
        }

    def __iter__(self):
        """Iterate through children."""
        for item in self.children:
            if isinstance(item, Field):
                yield item
            elif isinstance(item, Field_set):
                for field in item:
                    yield field

    def create_Field_set(self, outer, inner, name):
        """
        Create a Field_set object by supplying the inner and outer leaflet\
        Fields. These will be incorporated into a Field_set that then\
        calculates the symmetric "plus" and anti-symmetric "minus" fields,\
        respectively.

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
        for key, value in dict(self.children).items():
            if value == outer:
                del self.children[key]
            if value == inner:
                del self.children[key]
        self.children[name] = new_Field_set
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
            Field's traj attribute.
        name  :  str
            The name you want to give this Field.
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height',\
            'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo',\
            or 'zzero'.

        Returns
        -------
        new_Field  :  Field
            The new Field object you just created.

        """
        new_Field = Field(path, name, self, quantity, leaflet)
        self.children[name] = new_Field
        return new_Field

    def plot2d(self, obj):
        """
        Plot a two-dimensional heatmap. If Frame supplied, plots Frame values.\
        If Field or Trajectory supplied, plots average over time.

        Parameters
        ----------
        obj : Field, Trajectory, Frame
            The object whose data you want to plot.

        Returns
        -------
        None.

        """
        if isinstance(obj, Field):
            data = obj.traj.avg()
        elif isinstance(obj, Trajectory):
            data = obj.avg()
        elif isinstance(obj, Frame):
            data = obj.bins
        elif isinstance(obj, np.ndarray):
            if len(np.shape(obj)) == 3:
                assert np.shape(obj)[1] == m.grid_dims["N1"], "unexpected array size"
                assert np.shape(obj)[2] == m.grid_dims["N2"], "unexpected array size"
                data = calc_avg_over_time(obj)
            elif len(np.shape(obj)) == 2:
                assert np.shape(obj)[0] == m.grid_dims["N1"], "unexpected array size"
                assert np.shape(obj)[1] == m.grid_dims["N2"], "unexpected array size"
                data = obj
            else:
                raise Exception("This is not a 2D or 3D array.")
        else:
            raise Exception("I don't recognize this dtype")

        hmap_dims = bin_prep(self.grid_dims, self.polar)
        fig, ax = plot_maker(hmap_dims, data, False, False, self.polar)
        return fig, ax

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
        apart = field1.traj.avg() * field2.traj.avg()
        corr = together - apart
        return self.create_Field(corr, "corr_" + field1.name + "_" + field2.name)

    def measure_rms(self, field, eq2=None):
        """
        Measure the root-mean-squared value of some Field. Optionally, supply\
        a squared equilibrium value in order to generate a fold-enrichment\
        score.

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

    def dump(self, path):
        """
        Print all trajectories and averages to file.

        Parameters
        ----------
        path : Path or str
            The path to the directory where you want to dump your files.

        Returns
        -------
        None.

        """
        assert self.children.keys()
        for file_type in ["trajectory", "average"]:
            dir_name = path.joinpath(file_type)
            dir_name.mkdir(parents=True, exist_ok=True)
            for obj in self.children.values():
                obj_name = obj.name
                if isinstance(obj, Field_set):
                    subdir_name = dir_name.joinpath(obj_name)
                    subdir_name.mkdir(parents=True, exist_ok=True)
                    for field in obj:
                        field.save_to_file(subdir_name, file_type)
                elif isinstance(obj, Field):
                    obj.save_to_file(dir_name, file_type)
                else:
                    raise Exception(f"{obj_name} is not a Field or a Field_set.")
        return


class Field:
    """A field contains a measurement of some surface over the course of an MD\
    trajectory. This could be the height of the outer leaflet, the mean\
    curvature of the bilayer midplane, etc...

    Attributes
    ----------
    traj  :  Trajectory
        A 1D array containing Frame data from an MD trajectory. Could be \
        height values, curvature values, etc.
    name  :  str
        The name of the Field, as it is listed in its parent's children list.
    """

    def __init__(self, data, name, parent, quantity=None, leaflet=None):
        """
        Construct a Field object.

        Parameters
        ----------
        data  :  Path, str, or ndarray
            Either contains a path to TCL output data that needs to be parsed,\
            or contains a numpy ndarray that should just be saved into the\
            Field's Trajectory object.
        name  :  str
            The name of this Field
        parent  :  Membrane
            The Membrane object that this field belongs to.
        quantity  :  str
            If used, must contain a valid nougat quantity I.E. 'height',\
            'order', etc.
        leaflet  :  str
            If used, must contain a valid nougat leaflet I.E. 'zone', 'ztwo',\
            or 'zzero'.

        """
        self.name = name

        err_msg = "This ndarray doesn't have the same dimensions as its parent Membrane."

        # read in the data
        if isinstance(data, (Path, str)):
            assert quantity is not None, "quantity is required in order to use a path"
            assert leaflet is not None, "leaflet name is required in order to use a path"
            self.traj = Trajectory(self._parse_tcl_output(data, parent, quantity, leaflet), parent.polar)
        elif isinstance(data, np.ndarray):
            if len(np.shape(data)) == 1:
                # this is a Trajectory object
                data = np.stack(data)
            if len(np.shape(data)) == 2:
                # this is a single frame
                data = data[None, :, :]
            assert len(np.shape(data)) == 3, "A trajectory must have 3 dimensions."
            if parent.grid_dims["N1"] is not None:
                assert parent.grid_dims["N1"] == np.shape(data)[1], err_msg
                assert parent.grid_dims["N2"] == np.shape(data)[2], err_msg
            else:
                parent.grid_dims["N1"] = np.shape(data)[1]
                parent.grid_dims["N2"] = np.shape(data)[2]
                parent.grid_dims["Nframes"] = np.shape(data)[0]
            self.traj = Trajectory(data, parent.polar)
        else:
            raise ValueError("data must either be a numpy ndarray or a path")

    def __iter__(self):
        """Return self so that Membrane.children can be looped easily."""
        return self

    def __str__(self):
        """Say your name, rather than your address."""
        return self.name

    def __repr__(self):
        """Say your name, rather than your address."""
        return self.name

    def __array__(self):
        """Make the Trajectory accessible to numpy."""
        return self.traj

    def _parse_tcl_output(self, path, parent, quantity, leaflet):
        """
        Read in the tcl output data, update the parent Membrane's grid_dims,\
        and generate the traj array.

        Parameters
        ----------
        path : Path or str
            The path to the nougat.tcl results folder.
        parent  :  Membrane
            The Membrane object that this field belongs to.
        quantity : str
            A valid nougat.tcl output quantity, E.G. "height", "order", etc.
        leaflet : str
            A valid nougat.tcl leaflet name, I.E. "zone", "ztwo", or "zzero".

        Returns
        -------
        ndarray
            A 3D array containing nougat.tcl output data.
        """
        # import traj values
        input_file_path = path.joinpath("tcl_output", quantity, leaflet + ".dat")
        unrolled_data = np.genfromtxt(input_file_path, missing_values='nan', filling_values=np.nan)

        config_file_path = path.joinpath("tcl_output", "nougat.log")
        system_dict = read_log(config_file_path)
        N1, N2, d1, d2 = system_dict['bin_info'].values()

        # determine Nframes
        Nframes = int(np.shape(unrolled_data)[0] / N1)

        # error checks
        err_msg = "This ndarray doesn't have the same dimensions as its parent Membrane."
        if parent.grid_dims["N1"] is not None:
            assert parent.grid_dims["N1"] == N1, err_msg
            assert parent.grid_dims["d1"] == d1, err_msg
        else:
            parent.grid_dims["N1"] = N1
            parent.grid_dims["d1"] = d1
        if parent.grid_dims["N2"] is not None:
            assert parent.grid_dims["N2"] == N2, err_msg
            assert parent.grid_dims["d2"] == d2, err_msg
        else:
            parent.grid_dims["N2"] = N2
            parent.grid_dims["d2"] = d2
        if parent.grid_dims["Nframes"] is not None:
            assert parent.grid_dims["Nframes"] == Nframes, err_msg
        else:
            parent.grid_dims["Nframes"] = Nframes

        # create a new array that has each frame in a different array level
        field_data = np.zeros((Nframes, N1, N2))
        for frm in range(Nframes):
            field_data[frm, :, :] = unrolled_data[frm * N1: (frm + 1) * N1, 2:]

        field_data = mostly_empty(field_data)
        return field_data

    def save_to_file(self, path, file_type, name=None):
        """
        Save the trajectory or average to file.

        Parameters
        ----------
        path : Path or str
            The path to the directory where you would like to save the file.
        file_type : str
            "trajectory" or "average".
        name : str
            The file name you wish to use, without a suffix. Default is the \
            name of the field (self.name).

        Returns
        -------
        None.

        """
        if name is None:
            name = self.name
        if file_type == "trajectory":
            np.save(path.joinpath(name + ".npy"), self.traj._traj_to_3darray())
        elif file_type == "average":
            np.savetxt(path.joinpath(name + ".dat"), np.round(self.traj.avg(), decimals=5), delimiter=",")

    # BASIC MATH MAGIC METHODS BELOW #
    # These make it so that you can do math on the Field object, rather than\
    # having to specify the object's .traj attribute every time.

    def __add__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, Field):
            return self.traj + other.traj
        elif isinstance(other, (np.ndarray, int, float)):
            return self.traj + other
        else:
            return NotImplemented

    def __radd__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.traj + other
        else:
            return NotImplemented

    def __sub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, Field):
            return self.traj - other.traj
        elif isinstance(other, (np.ndarray, int, float)):
            return self.traj - other
        else:
            return NotImplemented

    def __rsub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, (np.ndarray, int, float)):
            return other - self.traj
        else:
            return NotImplemented

    def __mul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, Field):
            return self.traj * other.traj
        elif isinstance(other, (np.ndarray, int, float)):
            return self.traj * other
        else:
            return NotImplemented

    def __rmul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.traj * other
        else:
            return NotImplemented

    def __div__(self, other):
        """Use numpy to divide things."""
        if isinstance(other, Field):
            return self.traj / other.traj
        elif isinstance(other, (np.ndarray, int, float)):
            return self.traj / other
        else:
            return NotImplemented

    def __pow__(self, exponent):
        """Use numpy's power() on the array stored in this Field."""
        return np.power(self.traj, exponent)


class Trajectory:
    """
    A 1D array of Frames.

    Attributes
    ----------
    frames  :  np.ndarray
        1D array of the Frame objects that constitute the trajectory, in time-\
        order.
    polar  :  bool
        If True, allow averaging over theta dimension.
    """

    def __init__(self, frames, polar):
        """
        Construct a Trajectory.

        Parameters
        ----------
        frames : np.ndarray
            Contains data from one or multiple frames in an ndarray. Will be\
            converted to Frame objects.
        polar : bool
            If True, allow averaging over theta dimension.

        Returns
        -------
        None.

        """
        assert isinstance(frames, np.ndarray), "frames must be a numpy ndarray"
        if len(np.shape(frames)) == 2:
            self.frames = np.empty(1, dtype=Frame)
            self.frames[0] = Frame(0, frames)
        elif len(np.shape(frames)) == 3:
            self.frames = np.empty(np.shape(frames)[0], dtype=Frame)
            for f in range(np.shape(frames)[0]):
                self.frames[f] = Frame(f, frames[f, :, :])
        self.polar = polar

    def __len__(self):
        """Trajectory length = number of frames in trajectory."""
        return np.shape(self.frames)[0]

    def __array__(self):
        """Return the frames to numpy."""
        return self.frames

    def __iter__(self):
        """Iterate through all the Frames in the Trajectory."""
        for f in self.frames:
            yield f

    def __str__(self):
        """Print out all the Frames in the Trajectory."""
        for frame in self:
            print(frame)
        return ""

    def __repr__(self):
        """Return the .frames attribute."""
        return self.frames

    def _traj_to_3darray(self):
        """Turn Trajectory into 3D numpy array."""
        frames = self.frames
        frmlist = []
        for frame in frames:
            frmlist.append(frame.bins)
        frmlist = tuple(frmlist)
        return np.stack(frmlist, axis=0)

    def avg(self):
        """Calculate the trajectory average (over time)."""
        data_array = self._traj_to_3darray()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            avg = np.nanmean(data_array, axis=0)
            return avg

    def stdev(self):
        """Calculate standard deviation of time average."""
        data_array = self._traj_to_3darray()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stdev = np.nanstd(data_array, axis=0)
            return stdev

    def avg_over_theta(self):
        """Calculate the average over theta."""
        assert self.polar, "Trajectory must be polar in order to average over theta."
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            avg_over_theta = np.nanmean(self.avg(), axis=1)
            return avg_over_theta

    def stdev_over_theta(self):
        """Calculate standard deviation of time averages, averaged over theta."""
        assert self.polar, "Trajectory must be polar in order to average over theta."
        stdev2 = np.square(self.stdev())
        num_samples = np.sum(~np.isnan(stdev2), axis=1)
        num_samples[num_samples == 0] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sums_over_theta = np.nansum(stdev2, axis=1)
            sqrt_over_theta = np.sqrt(sums_over_theta)
            avg_stdev_over_theta = sqrt_over_theta / num_samples
            return avg_stdev_over_theta

    def __getitem__(self, item):
        """Make Trajectory object subscriptable."""
        if isinstance(item, int):
            return self.frames[item]
        else:
            return NotImplemented

    # BASIC MATH MAGIC METHODS BELOW #
    # These make it so that you can do math on the Trajectory object, rather\
    # than having to specify the object's .frames attribute every time.

    def __add__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, Trajectory):
            return self.frames + other.frames
        elif isinstance(other, (np.ndarray, int, float)):
            return self.frames + other
        else:
            return NotImplemented

    def __radd__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.frames + other
        else:
            return NotImplemented

    def __sub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, Trajectory):
            return self.frames - other.frames
        elif isinstance(other, (np.ndarray, int, float)):
            return self.frames - other
        else:
            return NotImplemented

    def __rsub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, Trajectory):
            return other.frames - self.frames
        elif isinstance(other, (np.ndarray, int, float)):
            return other - self.frames
        else:
            return NotImplemented

    def __mul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, Trajectory):
            return self.frames * other.frames
        elif isinstance(other, (np.ndarray, int, float)):
            return self.frames * other
        else:
            return NotImplemented

    def __rmul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.frames * other
        else:
            return NotImplemented

    def __div__(self, other):
        """Use numpy to divide things."""
        if isinstance(other, Trajectory):
            return self.frames / other.frames
        elif isinstance(other, (np.ndarray, int, float)):
            return self.frames / other
        else:
            return NotImplemented

    def __pow__(self, exponent):
        """Use numpy's power() on the array stored in this Field."""
        return np.power(self.frames, exponent)


class Frame:
    """A Frame is a 2D array of binned data.

    Attributes
    ----------
    index  :  int
        The Trajectory index at which this frame belongs.
    bins  :  np.ndarray
        The 2D array that contains either scalar or vector values.
    """

    def __init__(self, index, bins):
        """
        Create a Frame object.

        Parameters
        ----------
        index : int
            The index value of this frame in its parent Trajectory.
        bins : ndarray
            Binned data from one frame of a trajectory. Can contain scalar or \
            3-vector values.

        Returns
        -------
        None.

        """
        self.index = index
        assert isinstance(bins, np.ndarray), "bins must be a numpy ndarray."
        self.bins = bins

    def __iter__(self):
        """Iterate through all the Bins in the Frame."""
        for i in range(np.shape(self.bins)[0]):
            for j in range(np.shape(self.bins)[1]):
                yield self.bins[i, j]

    def __getitem__(self, item):
        """Make Frame object subscriptable."""
        if isinstance(item, tuple):
            assert len(item) == 2
            return self.bins[item[0], item[1]]
        else:
            return NotImplemented

    def __str__(self):
        """Print out the bins numpy style."""
        return str(self.bins)

    def __array__(self):
        """Make the bins accessible to numpy easily."""
        return self.bins

    # BASIC MATH MAGIC METHODS BELOW #
    # These make it so that you can do math on the frame object, rather than\
    # having to specify the object's .bins attribute every time.

    def __add__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, Frame):
            return self.bins + other.bins
        elif isinstance(other, (np.ndarray, int, float)):
            return self.bins + other
        else:
            return NotImplemented

    def __radd__(self, other):
        """Use numpy to add things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.bins + other
        else:
            return NotImplemented

    def __sub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, Frame):
            return self.bins - other.bins
        elif isinstance(other, (np.ndarray, int, float)):
            return self.bins - other
        else:
            return NotImplemented

    def __rsub__(self, other):
        """Use numpy to subtract things."""
        if isinstance(other, Frame):
            return other.bins - self.bins
        elif isinstance(other, (np.ndarray, int, float)):
            return other - self.bins
        else:
            return NotImplemented

    def __mul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, Frame):
            return self.bins * other.bins
        elif isinstance(other, (np.ndarray, int, float)):
            return self.bins * other
        else:
            return NotImplemented

    def __rmul__(self, other):
        """Use numpy to multiply things together."""
        if isinstance(other, (np.ndarray, int, float)):
            return self.bins * other
        else:
            return NotImplemented

    def __div__(self, other):
        """Use numpy to divide things."""
        if isinstance(other, Frame):
            return self.bins / other.bins
        elif isinstance(other, (np.ndarray, int, float)):
            return self.bins / other
        else:
            return NotImplemented

    def __pow__(self, exponent):
        """Use numpy's power() on the array stored in this Field."""
        return np.power(self.bins, exponent)


class Field_set:
    """A Field_set contains four Fields: the outer and inner leaflets, plus\
    the symmetric and anti-symmetric variables <Field>_plus and <Field>_minus.\
    This object serves as a constructor for the _plus and _minus variables.

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
            The name of this Field set, and the prefix that will be given to\
            the _plus and _minus fields.
        parent  :  Membrane
            The Membrane object to which this Field_set belongs.

        Returns
        -------
        None.
        """
        self.outer = outer
        self.inner = inner
        self.name = name
        self.plus = Field((outer + inner) / 2., self.name + "_plus", parent)
        self.minus = Field((outer - inner) / 2., self.name + "_minus", parent)

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
