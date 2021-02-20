# -*- coding: utf-8 -*-
# Copyright (c) 2011 Aaron Parsons
# Licensed under the GLPv3

import pytest
import random
import ephem
import aipy
import numpy as np

DIM = 128


def test_clean2dc():
    res = np.zeros((DIM, DIM), dtype=np.complex64)
    ker = np.zeros((DIM, DIM), dtype=np.complex64)
    mdl = np.zeros((DIM, DIM), dtype=np.complex64)
    area = np.zeros((DIM, DIM), dtype=np.int64)
    with pytest.raises(ValueError):
        aipy._deconv.clean(res, ker, mdl, area.astype(np.float64))
    ker[0, 0] = 1.0
    res[0, 0] = 1.0
    res[5, 5] = 1.0
    area[:4, :4] = 1
    rv = aipy._deconv.clean(res, ker, mdl, area, tol=1e-8)
    assert np.allclose(res[0, 0], 0, atol=1e-3)
    assert res[5, 5] == 1

    return


def test_clean2dr():
    res = np.zeros((DIM, DIM), dtype=np.float32)
    ker = np.zeros((DIM, DIM), dtype=np.float32)
    mdl = np.zeros((DIM, DIM), dtype=np.float32)
    area = np.zeros((DIM, DIM), dtype=np.int64)
    with pytest.raises(ValueError):
        aipy._deconv.clean(res, ker, mdl, area.astype(np.float32))
    ker[0, 0] = 1.0
    res[0, 0] = 1.0
    res[5, 5] = 1.0
    area[:4, :4] = 1
    rv = aipy._deconv.clean(res, ker, mdl, area, tol=1e-8)
    assert np.allclose(res[0, 0], 0, atol=1e-3)

    return


def test_clean1dc():
    res = np.zeros((DIM,), dtype=np.complex64)
    ker = np.zeros((DIM,), dtype=np.complex64)
    mdl = np.zeros((DIM,), dtype=np.complex64)
    area = np.zeros((DIM,), dtype=np.int64)
    with pytest.raises(ValueError):
        aipy._deconv.clean(res, ker, mdl, area.astype(np.float32))
    ker[0] = 1.0
    res[0] = 1.0
    res[5] = 1.0
    area[:4] = 1
    rv = aipy._deconv.clean(res, ker, mdl, area, tol=1e-8)
    assert np.allclose(res[0], 0, atol=1e-3)
    assert res[5] == 1

    return


def test_clean1dr():
    res = np.zeros((DIM,), dtype=np.float32)
    ker = np.zeros((DIM,), dtype=np.float32)
    mdl = np.zeros((DIM,), dtype=np.float32)
    area = np.zeros((DIM,), dtype=np.int64)
    with pytest.raises(ValueError):
        aipy._deconv.clean(res, ker, mdl, area.astype(np.float32))
    ker[0] = 1.0
    res[0] = 1.0
    res[5] = 1.0
    area[:4] = 1
    rv = aipy._deconv.clean(res, ker, mdl, area, tol=1e-8)
    assert np.allclose(res[0], 0, atol=1e-3)
    assert res[5] == 1

    return


def test_clean2d_stop_if_div():
    DIM1, DIM2 = 1000, 250
    dim = np.random.normal(size=(DIM1, DIM2))
    dbm = np.random.normal(size=(DIM1, DIM2))
    mdl = np.zeros(dim.shape, dtype=dim.dtype)
    area = np.ones(dim.shape, dtype=np.int64)
    init_val = dim[0, 0]
    rv = aipy._deconv.clean(
        dim, dbm, mdl, area, gain=0.1, tol=1e-2, stop_if_div=0, maxiter=100
    )
    assert np.allclose(dim[0, 0], init_val, atol=1e-3)

    return
