"""
Created on Fri Jul 21 11:18:40 2023.

@author: js2746
"""

import pytest
import numpy as np
from pathlib import Path
from collections import namedtuple


Comparison = namedtuple("Comparison", "test_data ref_data")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIXTURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@pytest.fixture
def cwd():
    """
    Get current working directory.

    """
    return Path.cwd()


@pytest.fixture(scope='function', params=["cart", "polar"])
def coordsys(request):
    """
    Supply the proper coordinate system name to the test function requesting it.

    """
    return request.param


@pytest.fixture(scope='function', params=["zone", "ztwo", "zplus", "zzero"])
def surface4(request):
    """
    Supply all four surfaces to the test function requesting it.

    """
    return request.param


@pytest.fixture(scope='function', params=["zone", "ztwo"])
def surface2(request):
    """
    Supply z1 and z2 to the test function requesting it.

    """
    return request.param


@pytest.fixture(scope='function', params=["meancurvature", "height", "normal_vectors", "gausscurvature"])
def quantity(request):
    """
    Supply the quantity being measured to the test function requesting it.

    """
    return request.param


@pytest.fixture(scope='function', params=["E-protein", "flat"])
def system(request):
    """
    Supply the system being tested to the test function requesting it.

    """
    return request.param


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def make_py_paths(wd, system, coord, surf, quant, file_format):
    """
    Concatenate strings together to make the path to the correct test files.

    Parameters
    ----------
    wd: string
        Path to current working directory.
    system: string
        Which test system the test pertains to.
    coord: string
        Coordinate system; 'cart' or 'polar'
    surf: sring
        The membrane surface in question (z1, z2, z0, or z+)
    quant: string
        The quantity being measured(height, thickness, curvature, etc.)
    file_format: string
        "dat" or "npy"

    Returns
    -------
    test_input: string
        Path to the test data npy file.
    expected: string
        Path to the saved 'correct' data npy file.

    """
    if file_format == "npy":
        file_type = "trajectory"
    elif file_format == "dat":
        file_type = "average"
    if coord == "cart":
        coordsys_path = "_cart_5_5_0_-1_1"
    elif coord == "polar":
        coordsys_path = "_polar_3_12_0_-1_1"
    if system == "E-protein":
        directory = "E-protein_trajectory"
    elif system == "flat":
        directory = "flat_surface_test"
    if (quant == "meancurvature"):
        expected = wd.joinpath(directory, system + coordsys_path, file_type, "curvature", "mean", surf + "." + file_format)
        test_input = wd.joinpath(directory, "test" + coordsys_path, file_type, "curvature", "mean", surf + "." + file_format)
    elif (quant == "gausscurvature"):
        expected = wd.joinpath(directory, system + coordsys_path, file_type, "curvature", "gaussian", surf + "." + file_format)
        test_input = wd.joinpath(directory, "test" + coordsys_path, file_type, "curvature", "gaussian", surf + "." + file_format)
    elif quant == "normal_vectors":
        expected = wd.joinpath(directory, system + coordsys_path, file_type, "curvature", quant, surf + "." + file_format)
        test_input = wd.joinpath(directory, "test" + coordsys_path, file_type, "curvature", quant, surf + "." + file_format)
    elif (quant == "height") or (quant == "thickness"):
        expected = wd.joinpath(directory, system + coordsys_path, file_type, quant, surf + "." + file_format)
        test_input = wd.joinpath(directory, "test" + coordsys_path, file_type, quant, surf + "." + file_format)
    return Comparison(test_input, expected)


def make_tcl_paths(wd, system, coord, surf):
    """
    Concatenate strings together to make the path to the correct test files.

    Parameters
    ----------
    wd: string
        Path to current working directory.
    system: string
        Which test system the test pertains to.
    coord: string
        Coordinate system; 'cart' or 'polar'
    surf: sring
        The membrane surface in question (z1, z2, z0, or z+)

    Returns
    -------
    test_input: string
        Path to the test data npy file.
    expected: string
        Path to the saved 'correct' data npy file.

    """
    if coord == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/tcl_output"
    elif coord == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/tcl_output"
    if system == "E-protein":
        directory = "E-protein_trajectory"
    elif system == "flat":
        directory = "flat_surface_test"
    expected = wd.joinpath(directory, system + coordsys_path, "height", surf + ".dat")
    test_input = wd.joinpath(directory, system + coordsys_path, "height", surf + ".dat")
    return Comparison(test_input.resolve(), expected.resolve())


def arrays_equal(paths, filetype, tolerance):
    """
    Determine whether two arrays have identical elements.

    Parameters
    ----------
    paths: Comparison named tuple
        Named tuple of two path objects.
    filetype: string
        The filetype of f1 and f2
    tolerance: float
        The amount of tolerance you have for differences between files. Should \
            be a small number!

    Returns
    -------
    Bool
        Whether or not the two arrays contain identical elements.

    """
    if filetype == "npy":
        f1 = np.load(paths[0])
        f2 = np.load(paths[1])
    elif filetype == "dat":
        f1 = np.genfromtxt(paths[0], missing_values="nan", filling_values=np.nan)
        f2 = np.genfromtxt(paths[1], missing_values="nan", filling_values=np.nan)

    if tolerance == 0:
        return np.array_equal(f1, f2, equal_nan=True)
    elif tolerance > 0:
        return np.allclose(f1, f2, rtol=0, atol=tolerance, equal_nan=True)
    else:
        raise ValueError("Tolerance must be positive or 0.")


# %%%%%%%%%%%%%%%%%%%%%%%%%% TESTS ARE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test if TCL outputs match

def test_if_tcl_heights_match(cwd, coordsys, surface2, system):
    paths = make_tcl_paths(cwd, system, coordsys, surface2)
    assert arrays_equal(paths, 'dat', 1e-11)

# Still needed: density, order, tilt tests


# Test if python trajectory outputs match

def test_if_heights_and_curvatures_match(cwd, coordsys, surface4, quantity, system):
    paths = make_py_paths(cwd, system, coordsys, surface4, quantity, "npy")
    assert arrays_equal(paths, 'npy', 1e-11)


def test_if_thicknesses_match(cwd, coordsys, surface2, system):
    paths = make_py_paths(cwd, system, coordsys, surface2, "thickness", "npy")
    assert arrays_equal(paths, 'npy', 1e-11)


def test_whether_flat(cwd, coordsys):
    if coordsys == "cart":
        settings = "_5_5_0_-1_1"
    else:
        settings = "_3_12_0_-1_1"
    Hone = np.load(cwd.joinpath("flat_surface_test", "test_" + coordsys + settings, "trajectory", "curvature", "mean", "zone.npy"))
    Htwo = np.load(cwd.joinpath("flat_surface_test", "test_" + coordsys + settings, "trajectory", "curvature", "mean", "ztwo.npy"))
    Hplus = Hone + Htwo / 2.0
    avgHplus = np.nanmean(Hplus)
    assert avgHplus <= 0.000000000001 and avgHplus >= -0.000000000001


def test_whether_flat_gaussian(cwd, coordsys):
    if coordsys == "cart":
        settings = "_5_5_0_-1_1"
    else:
        settings = "_3_12_0_-1_1"
    Kone = np.load(cwd.joinpath("flat_surface_test", "test_" + coordsys + settings, "trajectory", "curvature", "gaussian", "zone.npy"))
    Ktwo = np.load(cwd.joinpath("flat_surface_test", "test_" + coordsys + settings, "trajectory", "curvature", "gaussian", "ztwo.npy"))
    Kplus = Kone + Ktwo / 2.0
    avgKplus = np.nanmean(Kplus)
    assert avgKplus <= 0.000000000001 and avgKplus >= -0.000000000001


@pytest.mark.xfail(strict=True)
def test_if_leaflets_are_distinct(cwd, coordsys, quantity, system):
    zone_test, _ = make_py_paths(cwd, system, coordsys, "zone", quantity, "npy")
    ztwo_test, _ = make_py_paths(cwd, system, coordsys, "ztwo", quantity, "npy")
    assert arrays_equal((zone_test, ztwo_test), 'npy', 0)


@pytest.mark.xfail(strict=True)
def test_if_leaflet_thicknesses_are_distinct(cwd, coordsys, system):
    zone_test, _ = make_py_paths(cwd, system, coordsys, "zone", "thickness", "npy")
    ztwo_test, _ = make_py_paths(cwd, system, coordsys, "ztwo", "thickness", "npy")
    assert arrays_equal((zone_test, ztwo_test), 'npy', 0)

# Still needed: order, tilt


# Test if python time-averages match

def test_if_avg_heights_and_curvatures_match(cwd, coordsys, surface4, quantity, system):
    if quantity != "normal_vectors":
        paths = make_py_paths(cwd, system, coordsys, surface4, quantity, "dat")
        assert arrays_equal(paths, 'dat', 1e-11)
    else:
        return

# Still needed: thickness, order, tilt, normal_vectors


# BELOW: old density tests that need to be re-implemented
"""
def test_if_densities_match(cwd, coordsys, surface2):
    if coordsys == "cart":
        coordsys_path = "_cart_5_5_0_-1_1"
    elif coordsys == "polar":
        coordsys_path = "_polar_3_12_0_-1_1"
    exp = cwd.joinpath("E-protein_trajectory/E-protein" + coordsys_path, "trajectory", "density", "DTPC", surface2 + ".npy")
    test = cwd.joinpath("E-protein_trajectory/test" + coordsys_path, "trajectory", "density", "DTPC", surface2 + ".npy")
    paths = Comparison(test, exp)
    assert arrays_equal(paths, 'npy', 1e-11)


def test_if_avg_densities_match(cwd, coordsys, surface2, system):
    paths = make_avg_paths(cwd, system, coordsys, surface2, "avgdensity")
    assert arrays_equal(paths, 'dat', 1e-11)
"""
