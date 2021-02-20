# -*- coding: utf-8 -*-
# Copyright (c) 2011 Aaron Parsons
# Licensed under the GPLv3

import numpy as np
import pytest
import aipy._dsp as _dsp


def test_testgrid1D_c():
    buf = np.zeros(32, dtype=np.complex64)
    ind = np.array([5, 10.1, 14.9], dtype=np.float32)
    dat = np.array([1, 1, 1], dtype=np.complex64)
    _dsp.grid1D_c(buf, ind, dat)
    x = np.arange(32)
    ans = (
        np.exp(-((x - 5) ** 2) / (2 * 0.5 ** 2)) / np.sqrt(2 * np.pi * 0.5 ** 2)
        + np.exp(-((x - 10.1) ** 2) / (2 * 0.5 ** 2)) / np.sqrt(2 * np.pi * 0.5 ** 2)
        + np.exp(-((x - 14.9) ** 2) / (2 * 0.5 ** 2)) / np.sqrt(2 * np.pi * 0.5 ** 2)
    )
    assert np.allclose(np.max(np.abs(buf - ans)), 0, atol=1e-6)


def test_testgrid2D_c():
    buf = np.zeros((32, 32), dtype=np.complex64)
    ind = np.array([[5, 5], [10.1, 10.1], [14.5, 15.5]], dtype=np.float32)
    ind1 = ind[:, 0].copy()  # copy necessary b/c not contigous memory otherwise
    ind2 = ind[:, 1].copy()
    dat = np.array([1, 1, 1], dtype=np.complex64)
    _dsp.grid2D_c(buf, ind1, ind2, dat)
    assert np.allclose(buf[5, 5], 0.63661977236758149, atol=1e-6)
    assert np.allclose(
        buf[10, 10], 0.63661977236758149 * np.exp(-2 * (0.1 ** 2 + 0.1 ** 2)), atol=1e-6
    )
    assert np.allclose(
        buf[14, 15], 0.63661977236758149 * np.exp(-2 * (0.5 ** 2 + 0.5 ** 2)), atol=1e-6
    )


def test_testdegrid2D_c():
    buf = np.ones((32, 32), dtype=np.complex64)
    ind = np.array([[5, 5], [10.1, 10.1], [14.5, 15.5]], dtype=np.float32)
    ind1 = ind[:, 0].copy()  # copy necessary b/c not contigous memory otherwise
    ind2 = ind[:, 1].copy()
    dat = np.zeros(ind1.shape, dtype=np.complex64)
    _dsp.degrid2D_c(buf, ind1, ind2, dat)
    assert np.allclose(dat, 1.0)
