"""
Created on Fri Jul 21 11:18:40 2023.

@author: js2746
"""

import pytest
import numpy as np
import os


@pytest.fixture
def cwd():
    return os.getcwd()


@pytest.fixture(scope='function', params=["cart", "polar"])
def coordsys(request):
    return request.param


@pytest.fixture(scope='function', params=["zone", "ztwo", "zplus", "zzero"])
def surface(request):
    return request.param


@pytest.fixture(scope='function', params=["meancurvature", "height", "normal_vectors", "gausscurvature"])
def quantity(request):
    return request.param


def test_if_height_and_curvature_files_match(cwd, coordsys, surface, quantity):
    if coordsys == "cart":
        coordsys_path = "_cart_5_5_0_-1_1/npy/"
    elif coordsys == "polar":
        coordsys_path = "_polar_3_12_0_-1_1/npy/"
    expected = cwd + "/E-protein_trajectory/E-protein" + coordsys_path + "E-protein." + surface + ".C1A.C1B." + coordsys + "." + quantity + ".npy"
    test_input = cwd + "/E-protein_trajectory/test" + coordsys_path + "test." + surface + ".C1A.C1B." + coordsys + "." + quantity + ".npy"
    f1 = np.load(test_input)
    f2 = np.load(expected)
    assert f1.all() == f2.all()


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
