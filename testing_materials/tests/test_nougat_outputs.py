"""
Created on Fri Jul 21 11:18:40 2023.

@author: js2746
"""

import pytest
import numpy as np
import os


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIXTURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@pytest.fixture
def cwd():
    """
    Get current working directory.

    """
    return os.getcwd()


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


@pytest.fixture(scope='function', params=["curvature", "height", "Kcurvature"])
def avg_quantities(request):
    """
    Supply the average quantities being compared to the test function requesting it.

    """
    return request.param


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def make_npy_paths(wd, system, coord, surf, quant):
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
        The membrane surface in question (z1, z2, z0, or z +)
    quant: string
        The quantity being measured(height, thickness, curvature, etc.)

    Returns
    -------
    test_input: string
        Path to the test data npy file.
    expected: string
        Path to the saved 'correct' data npy file.

    """
    if coord == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/npy/"
    elif coord == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/npy/"
    if system == "E-protein":
        directory = "/E-protein_trajectory/"
    elif system == "flat":
        directory = "/flat_surface_test/"
    expected = wd + directory + system + coordsys_path + system + "." + surf + ".C1A.C1B." + coord + "." + quant + ".npy"
    test_input = wd + directory + "test" + coordsys_path + "test." + surf + ".C1A.C1B." + coord + "." + quant + ".npy"
    return test_input, expected


def make_avg_paths(wd, system, coord, surf, quant):
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
        The membrane surface in question (z1, z2, z0, or z +)
    quant: string
        The quantity being measured(height, thickness, curvature, etc.)

    Returns
    -------
    test_input: string
        Path to the test data npy file.
    expected: string
        Path to the saved 'correct' data npy file.

    """
    if coord == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/dat/"
    elif coord == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/dat/"
    if system == "E-protein":
        directory = "/E-protein_trajectory/"
    elif system == "flat":
        directory = "/flat_surface_test/"
    if quant != "avgdensity":
        expected = wd + directory + system + coordsys_path + system + "." + surf + ".C1A.C1B." + coord + ".avg" + quant + ".dat"
        test_input = wd + directory + "test" + coordsys_path + "test." + surf + ".C1A.C1B." + coord + ".avg" + quant + ".dat"
    elif quant == "avgdensity":
        expected = wd + directory + system + coordsys_path + system + ".DTPC." + surf + "." + coord + "." + quant + ".dat"
        test_input = wd + directory + "test" + coordsys_path + "test." + "DTPC." + surf + "." + coord + "." + quant + ".dat"
    return test_input, expected


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
        The membrane surface in question (z1, z2, z0, or z +)

    Returns
    -------
    test_input: string
        Path to the test data npy file.
    expected: string
        Path to the saved 'correct' data npy file.

    """
    if coord == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/tcl_output/"
    elif coord == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/tcl_output/"
    if system == "E-protein":
        directory = "/E-protein_trajectory/"
    elif system == "flat":
        directory = "/flat_surface_test/"
    expected = wd + directory + system + coordsys_path + system + "." + surf + ".C1A.C1B." + coord + ".height.dat"
    test_input = wd + directory + "test" + coordsys_path + "test." + surf + ".C1A.C1B." + coord + ".height.dat"
    return test_input, expected


def arrays_equal(path1, path2, filetype, tolerance):
    """
    Determine whether two arrays have identical elements.

    Parameters
    ----------
    path1: string
        path to first npy file.
    path2: string
        path to second npy file.
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
        f1 = np.load(path1)
        f2 = np.load(path2)
    elif filetype == "dat":
        f1 = np.genfromtxt(path1, missing_values="nan", filling_values=np.nan)
        f2 = np.genfromtxt(path2, missing_values="nan", filling_values=np.nan)

    if tolerance == 0:
        return np.array_equal(f1, f2, equal_nan=True)
    elif tolerance > 0:
        return np.allclose(f1, f2, rtol=0, atol=tolerance, equal_nan=True)
    else:
        raise ValueError("Tolerance must be positive or 0.")


# %%%%%%%%%%%%%%%%%%%%%%%%%% TESTS ARE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test if TCL outputs match

def test_if_tcl_heights_match(cwd, coordsys, surface2, system):
    test_input, expected = make_tcl_paths(cwd, system, coordsys, surface2)
    assert arrays_equal(test_input, expected, 'dat', 1e-11)

# Still needed: density, order, tilt tests


# Test if npy outputs match

def test_if_heights_and_curvatures_match(cwd, coordsys, surface4, quantity, system):
    test_input, expected = make_npy_paths(cwd, system, coordsys, surface4, quantity)
    assert arrays_equal(test_input, expected, 'npy', 1e-11)


def test_if_thicknesses_match(cwd, coordsys, surface2, system):
    test_input, expected = make_npy_paths(cwd, system, coordsys, surface2, "thickness")
    assert arrays_equal(test_input, expected, 'npy', 1e-11)


def test_if_densities_match(cwd, coordsys, surface2):
    if coordsys == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/npy/"
    elif coordsys == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/npy/"
    exp = cwd + "/E-protein_trajectory/E-protein" + coordsys_path + "E-protein.DTPC." + surface2 + "." + coordsys + ".density.npy"
    test = cwd + "/E-protein_trajectory/test" + coordsys_path + "test.DTPC." + surface2 + "." + coordsys + ".density.npy"
    assert arrays_equal(test, exp, 'npy', 1e-11)


@pytest.mark.xfail(strict=True)
def test_if_leaflets_are_distinct(cwd, coordsys, quantity, system):
    zone_test, _ = make_npy_paths(cwd, system, coordsys, "zone", quantity)
    ztwo_test, _ = make_npy_paths(cwd, system, coordsys, "ztwo", quantity)
    assert arrays_equal(zone_test, ztwo_test, 'npy', 0)


@pytest.mark.xfail(strict=True)
def test_if_leaflet_thicknesses_are_distinct(cwd, coordsys, system):
    zone_test, _ = make_npy_paths(cwd, system, coordsys, "zone", "thickness")
    ztwo_test, _ = make_npy_paths(cwd, system, coordsys, "ztwo", "thickness")
    assert arrays_equal(zone_test, ztwo_test, 'npy', 0)


def test_whether_flat_cartesian(cwd):
    Hone = np.load(cwd + "/flat_surface_test/test_cart_5_5_0_-1_1/npy/test.zone.C1A.C1B.cart.meancurvature.npy")
    Htwo = np.load(cwd + "/flat_surface_test/test_cart_5_5_0_-1_1/npy/test.ztwo.C1A.C1B.cart.meancurvature.npy")
    Hplus = Hone + Htwo / 2.0
    avgHplus = np.nanmean(Hplus)
    assert avgHplus <= 0.000000000001 and avgHplus >= -0.000000000001


def test_whether_flat_polar(cwd):
    Hone = np.load(cwd + "/flat_surface_test/test_polar_3_12_0_-1_1/npy/test.zone.C1A.C1B.polar.meancurvature.npy")
    Htwo = np.load(cwd + "/flat_surface_test/test_polar_3_12_0_-1_1/npy/test.ztwo.C1A.C1B.polar.meancurvature.npy")
    Hplus = Hone + Htwo / 2.0
    avgHplus = np.nanmean(Hplus)
    assert avgHplus <= 0.000000000001 and avgHplus >= -0.000000000001


# Still needed: order, tilt

# Test if time-averages match

def test_if_avg_heights_and_curvatures_match(cwd, coordsys, surface4, system, avg_quantities):
    test_input, expected = make_avg_paths(cwd, system, coordsys, surface4, avg_quantities)
    assert arrays_equal(test_input, expected, 'dat', 1e-11)


def test_if_avg_densities_match(cwd, coordsys, surface2, system):
    test_input, expected = make_avg_paths(cwd, system, coordsys, surface2, "avgdensity")
    assert arrays_equal(test_input, expected, 'dat', 1e-11)

# Still needed: thickness, order, tilt
