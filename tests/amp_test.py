# -*- coding: utf-8 -*-
# Copyright (c) 2011 Aaron Parsons
# Licensed under the GPLv3

import pytest
import aipy.amp as amp
import numpy as np


@pytest.fixture(scope="function")
def generate_AntennaArray():
    freqs = np.arange(0.1, 0.2, 0.01)
    bm = amp.Beam(freqs)
    ant0 = amp.Antenna(0, 0, 0, bm)
    aa = amp.AntennaArray(("0", "0"), [ant0])

    yield freqs, aa

    # clean up when done
    del freqs, aa

    return


@pytest.fixture(scope="function")
def generate_Beam():
    freqs = np.arange(0.1, 0.2, 0.01)
    bm = amp.Beam(freqs)

    yield freqs, bm

    # clean up when done
    del freqs, bm

    return


@pytest.fixture(scope="function")
def generate_Beam2DGaussian():
    freqs = np.arange(0.1, 0.2, 0.01)
    bm = amp.Beam2DGaussian(freqs, 0.05, 0.025)

    yield freqs, bm

    # clean up when done
    del freqs, bm

    return


@pytest.fixture(scope="function")
def generate_BeamAlm():
    freqs = np.arange(0.1, 0.2, 0.01)
    coeffs = {0: np.array([1.0, 0, 0], dtype=np.complex128)}
    bm = amp.BeamAlm(freqs, lmax=1, mmax=1, deg=0, coeffs=coeffs)

    yield freqs, coeffs, bm

    # clean up when done
    del freqs, coeffs, bm

    return


@pytest.fixture(scope="function")
def generate_Antenna():
    freqs = np.arange(0.1, 0.2, 0.01)
    bm = amp.Beam2DGaussian(freqs, 0.05, 0.025)
    ant = amp.Antenna(0, 0, 0, beam=bm)

    yield freqs, ant

    # clean up when done
    del freqs, ant

    return


def test_radiobody_attributes(generate_AntennaArray):
    """Test aipy.amp.RadioFixedBody attributes"""
    freqs, aa = generate_AntennaArray

    source = amp.RadioFixedBody(
        "0:00", "0:00", jys=100, index=-2, mfreq=0.1, name="src1"
    )
    assert source._jys == 100
    assert source.index == -2

    source.compute(aa)
    assert np.allclose(source.get_jys(), 100 * (freqs / 0.1) ** (-2))

    return


def test_response_TestBeam(generate_Beam):
    """Test retrieving aipy.amp.Beam response"""
    freqs, bm = generate_Beam
    xyz = (0, 0, 1)
    assert np.allclose(bm.response(xyz), np.ones_like(freqs))

    x = np.array([0, 0.1, 0.2])
    y = np.array([0.1, 0.2, 0])
    z = np.array([0.2, 0, 0.1])
    xyz = (x, y, z)
    assert np.allclose(bm.response(xyz), np.ones((freqs.size, 3)))

    bm.select_chans([0, 1, 2])
    assert np.allclose(bm.response(xyz), np.ones((3, 3)))

    return


def test_response_TestBeam2DGaussian(generate_Beam2DGaussian):
    """Test retrieving a 2D Gaussian beam response"""
    freqs, bm = generate_Beam2DGaussian
    xyz = (0, 0, 1)
    assert np.allclose(bm.response(xyz), np.ones_like(freqs))

    x = np.array([0, 0.05, 0])
    y = np.array([0, 0, 0.05])
    z = np.array([1, 1, 1])
    xyz = (x, y, z)
    resp = bm.response(xyz)
    assert resp.shape == (freqs.size, 3)

    ans = np.sqrt(np.array([1.0, np.exp(-0.5), np.exp(-2)]))
    ans.shape = (1, 3)
    bm.select_chans([0])
    resp = bm.response(xyz)
    assert np.allclose(np.round(resp - ans, 3), 0.0)

    return


def test_init_update_BeamAlm(generate_BeamAlm):
    freqs, coeffs, bm = generate_BeamAlm
    assert np.allclose(bm.hmap[0].map, 1.0 / np.sqrt(4 * np.pi))

    return


def test_response_BeamAlm(generate_BeamAlm):
    freqs, coeffs, bm = generate_BeamAlm
    x = np.array([0, 0.05, 0])
    y = np.array([0, 0, 0.05])
    z = np.array([1, 1, 1])
    xyz = (x, y, z)
    resp = bm.response(xyz)
    assert resp.shape == (freqs.size, 3)

    ans = np.ones(3, dtype=np.complex128) / np.sqrt(4 * np.pi)
    ans.shape = (1, 3)
    bm.select_chans([0])
    resp = bm.response(xyz)
    assert np.allclose(np.round(resp - ans, 3), 0.0)

    return


def test_passband_Antenna(generate_Antenna):
    """Test the Antenna passband"""
    freqs, ant = generate_Antenna
    pb = ant.passband()
    assert np.allclose(pb, np.ones_like(freqs))

    ant.select_chans([0, 1, 2])
    pb = ant.passband()
    assert pb.shape == (3,)

    return
