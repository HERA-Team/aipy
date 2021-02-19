# -*- coding: utf-8 -*-
# Copyright (c) 2010 Aaron Parsons
# Licensed under the GPLv3

import pytest
import numpy as np

from aipy import miriad
from aipy import _miriad

@pytest.fixture(scope="function")
def test_file_r(tmp_path):
    filename1 = str(tmp_path / "test1.uv")
    filename2 = str(tmp_path / "test2.uv")
    uv = miriad.UV(filename1, status="new", corrmode="r")
    uv["history"] = "Made this file from scratch.\n"
    uv.add_var("nchan", "i")
    uv.add_var("pol", "i")
    uv["nchan"] = 4
    uv["pol"] = -5
    uvw = np.array([1, 2, 3], dtype=np.float64)
    preamble = (uvw, 12345.6789, (0, 1))
    data = np.ma.array([1j, 2, 3j, 4], mask=[0, 0, 1, 0], dtype=np.complex64)
    uv.write(preamble, data)
    uv["pol"] = -6
    uv.write(preamble, data)
    del uv

    yield filename1, filename2, data

    return

@pytest.fixture(scope="function")
def test_file_j(tmp_path):
    filename = str(tmp_path / "test1.uv")
    uv = miriad.UV(filename, status="new", corrmode="j")
    uv["history"] = "Made this file from scratch.\n"
    uv.add_var("nchan", "i")
    uv.add_var("pol", "i")
    uv["nchan"] = 4
    uv["pol"] = -5
    uvw = np.array([1, 2, 3], dtype=np.float64)
    preamble = (uvw, 12345.6789, (0, 1))
    data = np.ma.array([1j, 2, 3j, 4], mask=[0, 0, 1, 0], dtype=np.complex64)
    uv.write(preamble, data)
    uv["pol"] = -6
    uv.write(preamble, data)
    del uv

    yield filename, data

    return

def test_maxchan():
    assert type(_miriad.MAXCHAN) == int
    return

def test_immediate_corr(test_file_r):
    """Test immediate corr of a Miriad UV file"""
    filename1, filename2, data = test_file_r
    uv = miriad.UV(filename2, status="new")
    assert uv.vartable["corr"] == "r"
    return

def test_vartable_r(test_file_r):
    """Test accesing vartable data in a Miriad UV file"""
    filename1, filename2, data = test_file_r
    uv1 = miriad.UV(filename1)
    assert uv1.vartable["corr"] == "r"
    assert uv1.vartable["nchan"] == "i"
    assert uv1.vartable["pol"] == "i"
    return

def test_data_r(test_file_r):
    """Test writing data from a Miriad UV file"""
    filename1, filename2, data = test_file_r
    uv = miriad.UV(filename1)
    assert uv["history"] == "Made this file from scratch.\n"
    (uvw, t, bl), d = uv.read()
    assert uv["nchan"] == 4
    assert uv["pol"] == -5
    assert bl == (0, 1)
    assert np.isclose(t, 12345.6789)
    assert np.allclose(uvw, np.array([1, 2, 3], dtype=np.float64))
    assert np.allclose(d, data)

    (uvw, t, bl), d = uv.read()
    assert uv["nchan"] == 4
    assert uv["pol"] == -6
    assert bl == (0, 1)
    assert np.isclose(t, 12345.6789)
    assert np.allclose(uvw, np.array([1, 2, 3], dtype=np.float64))
    assert np.allclose(d, data)
    return

def test_vartable_j(test_file_j):
    """Test accesing vartable data in a Miriad UV file"""
    filename, data = test_file_j
    uv = miriad.UV(filename, corrmode="j")
    assert uv.vartable["corr"] == "j"
    assert uv.vartable["nchan"] == "i"
    assert uv.vartable["pol"] == "i"
    return

def test_data_j(test_file_j):
    """Test writing data from a Miriad UV file"""
    filename, data = test_file_j
    uv = miriad.UV(filename)
    assert uv["history"] == "Made this file from scratch.\n"

    (uvw, t, bl), d = uv.read()
    assert uv["nchan"] == 4
    assert uv["pol"] == -5
    assert bl == (0, 1)
    assert np.isclose(t, 12345.6789)
    assert np.allclose(uvw, np.array([1, 2, 3], dtype=np.float64))
    assert np.allclose(np.abs(d - data), 0.0, atol=1e-4)

    (uvw, t, bl), d = uv.read()
    assert uv["nchan"] == 4
    assert uv["pol"] == -6
    assert bl == (0, 1)
    assert np.isclose(t, 12345.6789)
    assert np.allclose(uvw, np.array([1, 2, 3], dtype=np.float64))
    assert np.allclose(np.abs(d - data), 0.0, atol=1e-4)
    return
