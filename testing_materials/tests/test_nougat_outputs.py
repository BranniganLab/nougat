"""
Created on Fri Jul 21 11:18:40 2023.

@author: js2746
"""

import pytest
import numpy as np
import os


@pytest.fixture
def path():
    return os.getcwd()


def test_whether_flat_cartesian(path):
    Hone = np.load(path + "/flat_surface_test/test_cart_5_5_0_-1_1/npy/test.zone.C1A.C1B.cart.meancurvature.npy")
    Htwo = np.load(path + "/flat_surface_test/test_cart_5_5_0_-1_1/npy/test.ztwo.C1A.C1B.cart.meancurvature.npy")
    Hplus = Hone + Htwo / 2.0
    avgHplus = np.nanmean(Hplus)
    assert avgHplus <= 0.000000000001 and avgHplus >= -0.000000000001


def test_whether_flat_polar(path):
    Hone = np.load(path + "/flat_surface_test/test_polar_3_12_0_-1_1/npy/test.zone.C1A.C1B.polar.meancurvature.npy")
    Htwo = np.load(path + "/flat_surface_test/test_polar_3_12_0_-1_1/npy/test.ztwo.C1A.C1B.polar.meancurvature.npy")
    Hplus = Hone + Htwo / 2.0
    avgHplus = np.nanmean(Hplus)
    assert avgHplus <= 0.000000000001 and avgHplus >= -0.000000000001


def do_files_match(path1, path2):
    f1 = np.load(path1)
    f2 = np.load(path2)
    assert f1.all() == f2.all()


def test_E_protein_z1_height_cartesian(path):
    original = path + "/E-protein_trajectory/E-protein_cart_5_5_0_-1_1/npy/E-protein.zone.C1A.C1B.cart.height.npy"
    test_data = path + "/E-protein_trajectory/test_cart_5_5_0_-1_1/npy/test.zone.C1A.C1B.cart.height.npy"
    do_files_match(original, test_data)


def test_E_protein_z1_height_polar(path):
    original = path + "/E-protein_trajectory/E-protein_polar_3_12_0_-1_1/npy/E-protein.zone.C1A.C1B.polar.height.npy"
    test_data = path + "/E-protein_trajectory/test_polar_3_12_0_-1_1/npy/test.zone.C1A.C1B.polar.height.npy"
    do_files_match(original, test_data)


def test_E_protein_z1_curvature_cartesian(path):
    original = path + "/E-protein_trajectory/E-protein_cart_5_5_0_-1_1/npy/E-protein.zone.C1A.C1B.cart.meancurvature.npy"
    test_data = path + "/E-protein_trajectory/test_cart_5_5_0_-1_1/npy/test.zone.C1A.C1B.cart.meancurvature.npy"
    do_files_match(original, test_data)


def test_E_protein_z1_curvature_polar(path):
    original = path + "/E-protein_trajectory/E-protein_polar_3_12_0_-1_1/npy/E-protein.zone.C1A.C1B.polar.meancurvature.npy"
    test_data = path + "/E-protein_trajectory/test_polar_3_12_0_-1_1/npy/test.zone.C1A.C1B.polar.meancurvature.npy"
    do_files_match(original, test_data)
