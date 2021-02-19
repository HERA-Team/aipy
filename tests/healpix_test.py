# -*- coding: utf-8 -*-
# Copyright (c) 2018 Aaron Parsons
# Licensed under the GPLv3

import pytest
import aipy.healpix as ah
import numpy as np

@pytest.fixture(scope="function")
def healpix_base():
    hpb = ah.HealpixBase(nside=32, scheme="RING")
    yield hpb

    # cleanup
    del hpb
    return

@pytest.fixture(scope="function")
def healpix_map():
    hpm = ah.HealpixMap(32, "RING")
    yield hpm

    # cleanup
    del hpm
    return

def test_order(healpix_base):
    """Test HEALpix order functional attribute"""
    hpb = healpix_base
    assert hpb.order() == 5
    return

def test_nside(healpix_base):
    """Test HEALpix nside functional attribute"""
    hpb = healpix_base
    assert hpb.nside() == 32
    return

def test_npix(healpix_base):
    """Test HEALpix npix functional attribute"""
    hpb = healpix_base
    assert hpb.npix() == 12 * hpb.nside()**2
    return

def test_scheme(healpix_base):
    """Test HEALpix scheme functional attribute"""
    hpb = healpix_base
    assert hpb.scheme() == "RING"
    return

def test_npix2nside(healpix_base):
    """Test HEALpix npix2nside functional attribute"""
    hpb = healpix_base
    assert hpb.npix2nside(12 * 2**12) == 2**6
    return

def test_set_nside_scheme(healpix_base):
    hpb = healpix_base
    hpb.set_nside_scheme(64, "NEST")
    assert hpb.nside() == 64
    assert hpb.scheme() == "NEST"

    hpb.set_nside_scheme(32, "RING")
    assert hpb.nside() == 32
    assert hpb.scheme() == "RING"
    return

def test_nest_ring_conv(healpix_base):
    hpb = healpix_base
    ipx = np.array([0])
    px = hpb.nest_ring_conv(ipx, "NEST")
    assert px == 1023

    hpb.set_nside_scheme(32, "NEST")
    ipx = np.array([0])
    px = hpb.nest_ring_conv(ipx, "RING")
    assert px == 5968
    return

def test_ang2px(healpix_base):
    hpb = healpix_base
    th = np.linspace(0, np.pi, 10)
    ph = np.linspace(-np.pi, np.pi, 10)
    px = hpb.crd2px(th,ph)
    assert len(px) == len(th)
    np.testing.assert_allclose(
        px, np.array([2, 398, 1375, 3114, 5049, 7239, 9173, 10912, 11889, 12286])
    )
    return

def test_vec2px(healpix_base):
    hpb = healpix_base
    x = np.linspace(-0.5, 0.5, 10)
    y = np.linspace(-0.5, 0.5, 10)
    z = 1 - np.sqrt(x**2 + y**2)
    px = hpb.crd2px(x, y, z)
    assert len(px) == len(x)
    np.testing.assert_allclose(
        px, np.array([3728, 2192, 1069, 247, 19, 13,225, 1023, 2128, 3664])
    )
    return

def test_px2vec(healpix_base):
    hpb = healpix_base
    px = np.arange(hpb.npix())
    x, y, z = hpb.px2crd(px)
    assert len(px) == len(x)

    px_recovered = hpb.crd2px(x, y, z)
    np.testing.assert_equal(px, px_recovered)
    return

def test_px2ang(healpix_base):
    hpb = healpix_base
    px = np.arange(hpb.npix())
    th, phi = hpb.px2crd(px, ncrd=2)
    assert len(px) == len(th)

    px_recovered = hpb.crd2px(th, phi)
    np.testing.assert_equal(px, px_recovered)
    return

def test_ang2px_interp(healpix_base):
    hpb = healpix_base
    th = np.linspace(0, np.pi, 3)
    ph = np.linspace(-np.pi, np.pi, 3)
    px, wgt = hpb.crd2px(th, ph, interpolate=True)
    assert len(px) == len(th)
    assert wgt.shape == (len(th), 4)
    np.testing.assert_allclose(
        px, np.array(
            [
                [ 3, 0, 1, 2],
                [ 6207, 6080, 6208, 6209],
                [12285, 12286, 12287, 12284],
            ]
        ),
    )
    np.testing.assert_allclose(
        wgt,
        np.array(
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.5 , 0.5 , 0.,  0.],
                [0.25, 0.25, 0.25, 0.25],
            ]
        ),
    )
    return

def test_vec2px_interp(healpix_base):
    hpb = healpix_base
    x = np.linspace(-0.5, 0.5, 3)
    y = np.linspace(-0.5, 0.5, 3)
    z = 1 - np.sqrt(x**2 + y**2)
    px, wgt = hpb.crd2px(x, y, z, interpolate=True)
    assert len(px) == len(x)
    assert wgt.shape == (len(x), 4)
    np.testing.assert_allclose(
        px,
        np.array(
            [
                [3728, 3729, 3855, 3856],
                [1, 2, 3, 0],
                [3664, 3665, 3791, 3792],
            ]
        ),
    )
    np.testing.assert_allclose(
        wgt,
        np.array(
            [
                [0.367711, 0., 0.316145, 0.316145],
                [0.25, 0.25, 0.25, 0.25],
                [0.367711, 0., 0.316145, 0.316145],
            ]
        ),
        rtol=1e-5,
    )
    return

def test_set(healpix_map):
    hpm = healpix_map
    hpm[0] = 1
    assert hpm[0] == 1

    hpm[0, 0] = 1
    assert hpm[0,0] == 1

    hpm[0, 1, 0] = 1
    assert hpm[0, 1, 0] == 1
    return

def test_set_px(healpix_map):
    hpm = healpix_map
    cval = np.arange(15)
    assert len(hpm[cval]) == len(cval)

    hpm[cval] += cval
    np.testing.assert_allclose(hpm[cval], cval)

    hpm[cval] += cval
    np.testing.assert_allclose(hpm[cval], 2 * cval)

    hpm[cval] = cval
    np.testing.assert_allclose(hpm[cval], cval)
    return

def test_set_ang(healpix_map):
    hpm = healpix_map
    th = np.linspace(0, np.pi, 10)
    ph = np.linspace(-np.pi, np.pi, 10)
    assert len(hpm[th, ph]) == len(th)

    hpm[th, ph] += th
    np.testing.assert_allclose(hpm[th, ph], th)

    hpm[th, ph] += th
    np.testing.assert_allclose(hpm[th, ph], 2 * th)

    hpm[th, ph] = th
    np.testing.assert_allclose(hpm[th, ph], th)
    return

def test_set_vec(healpix_map):
    hpm = healpix_map
    x = np.linspace(-0.5, 0.5, 10)
    y = np.linspace(-0.5, 0.5, 10)
    z = 1 - np.sqrt(x**2 + y**2)
    assert len(hpm[x, y, z]) == len(x)

    hpm[x, y, z] += x
    np.testing.assert_allclose(hpm[x, y, z], x)

    hpm[x, y, z] += x
    np.testing.assert_allclose(hpm[x, y, z], 2 * x)

    hpm[x, y, z] = x
    np.testing.assert_allclose(hpm[x, y, z], x)
    return
