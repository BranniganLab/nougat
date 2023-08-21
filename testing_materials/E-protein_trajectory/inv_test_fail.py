#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 09:30:09 2023.

@author: js2746
"""
import numpy as np

exp = np.load("E-protein_polar_3_12_0_-1_1/npy/E-protein.ztwo.C1A.C1B.polar.height.npy")
test = np.load("test_polar_3_12_0_-1_1/npy/test.ztwo.C1A.C1B.polar.height.npy")
data = exp - test
comp = (exp == test) | np.isnan(exp) & np.isnan(test)
