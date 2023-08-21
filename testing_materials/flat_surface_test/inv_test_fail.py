#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 09:30:09 2023.

@author: js2746
"""
import numpy as np

exp = np.load("flat_polar_3_12_0_-1_1/npy/flat.ztwo.C1A.C1B.polar.meancurvature.npy")
test = np.load("test_polar_3_12_0_-1_1/npy/test.ztwo.C1A.C1B.polar.meancurvature.npy")
data = exp - test
print(exp == test)
