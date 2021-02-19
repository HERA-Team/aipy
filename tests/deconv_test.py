# -*- coding: utf-8 -*-
# Copyright (c) 2008 Aaron Parsons
# Licensed under the GPLv3

import pytest
import numpy as np
import aipy

@pytest.fixture(scope="function")
def init_deconv():
    SIZE = 100
    NOISE = 0.001
    img = np.zeros((SIZE, SIZE), np.float64)
    img[10, 10] = 10
    img[20:25, 20:25] = 1
    img[30:40, 30:40] = 0.1

    bm = aipy.img.gaussian_beam(2, shape=img.shape)
    bm[0, 0] = 0.05

    data = np.abs(np.fft.ifft2(np.fft.fft2(img) * np.fft.fft2(bm)))
    ns = np.random.normal(scale=NOISE, size=img.shape)
    data = np.abs(data + ns)

    yield data, bm

    # clean up when done
    del data, bm

    return

def test_clean(init_deconv):
    """Test that the standard clean deconvolution runs"""
    data, bm = init_deconv
    cln, info = aipy.deconv.clean(data, bm, verbose=False)
    assert True

    return

def test_lsq(init_deconv):
    """Test that least squared deconvolution runs"""
    data, bm = init_deconv
    cln, info = aipy.deconv.lsq(data, bm, verbose=False)
    assert True

    return

def test_mem(init_deconv):
    """Test the maximum entropy deconvolution runs"""
    data, bm = init_deconv
    cln, info = aipy.deconv.maxent(data, bm, np.var(data**2) * 0.5, verbose=False)
    assert True

    return

def test_anneal(init_deconv):
    """Test that simulated annealing deconvolution runs"""
    data, bm = init_deconv
    cln, info = aipy.deconv.anneal(data, bm, verbose=False)
    assert True

    return
