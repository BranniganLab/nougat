"""
Created on Fri Jul 21 11:18:40 2023.

@author: js2746
"""

import pytest
import numpy as np
import os
from itertools import product


@pytest.fixture
def cwd():
    return os.getcwd()


@pytest.fixture
def test_paths(cwd):
    coordsys = ["cart", "polar"]
    surface = [".zone.", ".ztwo.", ".zplus.", ".zzero."]
    quant = [".meancurvature.npy", ".height.npy"]
    paths = []

    for prod in product(coordsys, surface, quant):
        if prod[0] == "cart":
            coordsys_path = "_cart_5_5_0_-1_1/npy/"
        elif prod[0] == "polar":
            coordsys_path = "_polar_3_12_0_-1_1/npy/"
        expected = cwd + "/E-protein_trajectory/E-protein" + coordsys_path + "E-protein" + prod[1] + "C1A.C1B." + prod[0] + prod[2]
        test_data = cwd + "/E-protein_trajectory/test" + coordsys_path + "test" + prod[1] + "C1A.C1B." + prod[0] + prod[2]
        paths.append(tuple([expected, test_data]))
    return paths


def test_if_files_match(test_paths):
    for file_pair in test_paths:
        f1 = np.load(file_pair[0])
        f2 = np.load(file_pair[1])
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
