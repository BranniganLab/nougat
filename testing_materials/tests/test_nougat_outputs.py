"""
Created on Fri Jul 21 11:18:40 2023.

@author: js2746
"""

import pytest
import numpy as np
import os


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


def make_paths(wd, coord, surf, quant):
    """
    Concatenate strings together to make the path to the correct test files.

    Parameters
    ----------
    wd : string
        Path to current working directory.
    coord : string
        Coordinate system; 'cart' or 'polar'
    surf : sring
        The membrane surface in question (z1, z2, z0, or z+)
    quant : string
        The quantity being measured (height, thickness, curvature, etc.)

    Returns
    -------
    test_input : string
        Path to the test data npy file.
    expected : string
        Path to the saved 'correct' data npy file.

    """
    if coord == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/npy/"
    elif coord == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/npy/"
    expected = wd + "/E-protein_trajectory/E-protein" + coordsys_path + "E-protein." + surf + ".C1A.C1B." + coord + "." + quant + ".npy"
    test_input = wd + "/E-protein_trajectory/test" + coordsys_path + "test." + surf + ".C1A.C1B." + coord + "." + quant + ".npy"
    return test_input, expected


def arrays_equal(path1, path2):
    """
    Determine whether two arrays have identical elements.

    Parameters
    ----------
    path1 : string
        path to first npy file.
    path2 : string
        path to second npy file.

    Returns
    -------
    Bool
        Whether or not the two arrays contain identical elements.

    """
    f1 = np.load(path1)
    f2 = np.load(path2)
    return np.array_equal(f1, f2, equal_nan=True)


def test_if_height_and_curvature_files_match(cwd, coordsys, surface4, quantity):
    test_input, expected = make_paths(cwd, coordsys, surface4, quantity)
    assert arrays_equal(test_input, expected)


def test_if_thickness_files_match(cwd, coordsys, surface2):
    test_input, expected = make_paths(cwd, coordsys, surface2, "thickness")
    assert arrays_equal(test_input, expected)


@pytest.mark.xfail(strict=True)
def test_if_thickness_files_dont_match(cwd, coordsys, surface2):
    test_input, expected = make_paths(cwd, coordsys, surface2, "thickness")
    assert (arrays_equal(test_input, expected) is False)


@pytest.mark.xfail(strict=True)
def test_if_height_and_curvature_files_dont_match(cwd, coordsys, surface4, quantity):
    test_input, expected = make_paths(cwd, coordsys, surface4, quantity)
    assert (arrays_equal(test_input, expected) is False)


@pytest.mark.xfail(strict=True)
def test_if_leaflets_are_distinct(cwd, coordsys, quantity):
    zone_test, _ = make_paths(cwd, coordsys, "zone", quantity)
    ztwo_test, _ = make_paths(cwd, coordsys, "ztwo", quantity)
    assert arrays_equal(zone_test, ztwo_test)


@pytest.mark.xfail(strict=True)
def test_if_leaflet_thicknesses_are_distinct(cwd, coordsys):
    zone_test, _ = make_paths(cwd, coordsys, "zone", "thickness")
    ztwo_test, _ = make_paths(cwd, coordsys, "ztwo", "thickness")
    assert arrays_equal(zone_test, ztwo_test)


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


def test_if_densities_match(cwd, coordsys, surface2):
    if coordsys == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/npy/"
    elif coordsys == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/npy/"
    exp = cwd + "/E-protein_trajectory/E-protein" + coordsys_path + "E-protein.DTPC." + surface2 + "." + coordsys + ".density.npy"
    test = cwd + "/E-protein_trajectory/test" + coordsys_path + "test.DTPC." + surface2 + "." + coordsys + ".density.npy"
    assert arrays_equal(test, exp)
