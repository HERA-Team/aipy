# -*- coding: utf-8 -*-
# Copyright (c) 2011 Aaron Parsons
# Licensed under the GPLv3

import pytest
import aipy.healpix as ah
import numpy as np

def test_init():
    """Test initializing a aipy._alm object"""
    alm = ah.Alm(10, 5)
    assert alm.lmax() == 10
    assert alm.mmax() == 5
    with pytest.raises(Exception):
        ah.Alm(5, 10)

    return

def test_get_data():
    alm = ah.Alm(1,1)
    data = alm.get_data()
    assert data.size == 3

    return

def test_set_data():
    alm = ah.Alm(1,1)
    data = np.array([1 + 1j, 2 + 2j, 3 + 3j], dtype=np.complex128)
    alm.set_data(data)
    assert np.allclose(alm.get_data(), data)

    return

def test_set_to_zero():
    alm = ah.Alm(3,3)
    assert np.allclose(alm.get_data(), 0.0)
    return

def test_get_set():
    alm = ah.Alm(1, 1)
    data = np.array([1 + 1j, 2 + 2j, 3 + 3j], dtype=np.complex128)
    alm.set_data(data)
    assert alm[0, 0] == data[0]
    assert alm[1, 0] == data[1]
    assert alm[1, 1] == data[2]

    alm[1,1] = 4+4j
    assert alm[1,1] == 4+4j

    return

def test_to_map():
    alm = ah.Alm(1, 1)
    alm[0, 0] = 1
    d = alm.to_map(32, 'RING')
    assert np.allclose(d, np.ones(12 * 32**2, dtype=np.float32) * 0.2820948)

    return

def test_from_map():
    alm = ah.Alm(1,1)
    d = np.ones(12 * 32**2, dtype=np.float32)
    alm.from_map(d, 3)
    assert np.allclose(alm[0, 0], 3.5449077018110313)
    assert np.allclose(alm[1, 0], 0.0)
    assert np.allclose(alm[1, 1], 0.0)

    return
