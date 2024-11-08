#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:07:30 2024

@author: js2746
"""

import pytest
import sys
import os
import numpy as np
sys.path.append(os.path.abspath('../../python/'))
from nougat import *
from utils import *


def test_make_todo_list_None():
    todo_list = make_todo_list(None)
    testm = Membrane(True, todo_list)
    assert len(testm.to_analyze) == 4
    # change this number when you add more features!


def test_make_todo_list():
    todo_list = make_todo_list('htoc')
    testm = Membrane(True, todo_list)
    assert len(testm.to_analyze) == 4
    assert "height" in testm.to_analyze
    assert "thickness" in testm.to_analyze
    assert "order" in testm.to_analyze
    assert "curvature" in testm.to_analyze


def test_create_Field_numpy():
    a = np.zeros((4, 4))
    b = np.ones((4, 4))
    b[0, 0] = 7
    testm = Membrane(True, None)
    testf = testm.create_Field(a + b, "test")
    assert np.array_equal(testf.traj[0].bins, (a + b))
