# -*- coding: utf-8 -*-
# Copyright (c) 2018 Aaron Parsons
# Licensed under the GPLv3

import numpy as np
from aipy import twodgauss

def test_moments():
    """Test twodgauss.moments"""
    a0 = 2.0
    x0 = 5.0
    y0 = 7.0
    wx = 1.0
    wy = 1.0
    x = np.arange(11, dtype=np.float64)
    y = np.arange(11, dtype=np.float64)
    xx, yy = np.meshgrid(x, y)
    data = a0 * np.exp(-((xx - x0) / wx)**2 / 2.0 - ((yy - y0) / wy)**2 / 2.0)
    params = twodgauss.moments(data)
    assert np.isclose(a0, params[1], atol=1e-3)
    assert np.isclose(x0, params[2], atol=1e-3)
    assert np.isclose(y0, params[3], atol=1e-3)
    return

def test_twodgaussian():
    """Test twodgauss.twodgaussian"""
    a0 = 2.0
    x0 = 5.0
    y0 = 7.0
    wx = 1.0
    wy = 1.0
    x = np.arange(11, dtype=np.float64)
    y = np.arange(11, dtype=np.float64)
    xx, yy = np.meshgrid(x, y)
    data0 = a0 * np.exp(-((xx - x0) / wx)**2 / 2.0 - ((yy - y0) / wy)**2 / 2.0)
    data1 = twodgauss.twodgaussian([0.0, a0, x0, y0, wx, wy], shape=yy.shape)
    assert np.allclose(data0, data1, atol=1e-6)
    return
