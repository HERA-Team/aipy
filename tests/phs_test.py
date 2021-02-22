# -*- coding: utf-8 -*-
# Copyright (c) 2008 Aaron Parsons
# Licensed under the GPLv3

import random

import aipy
import pytest
import ephem
import numpy as np


@pytest.fixture(scope="function")
def source_catalog():
    src1 = aipy.phs.RadioFixedBody("1:00", "1:00", name="src1")
    src2 = aipy.phs.RadioFixedBody("2:00", "2:00", name="src2")
    src3 = aipy.phs.RadioFixedBody("3:00", "3:00", name="src3")
    srcs = [src1, src2, src3]
    cat = aipy.phs.SrcCatalog(srcs)
    yield srcs, cat
    return


@pytest.fixture(scope="function")
def test_beam():
    freqs = np.arange(0, 1, 0.1)
    bm = aipy.phs.Beam(freqs)
    yield freqs, bm
    return


@pytest.fixture(scope="function")
def test_antenna():
    freqs = np.arange(0, 1, 0.1)
    bm = aipy.phs.Beam(freqs)
    ant = aipy.phs.Antenna(1, 2, 3, bm, phsoff=[0, 1])
    yield freqs, bm, ant
    return


@pytest.fixture(scope="function")
def test_array_location():
    yield aipy.phs.ArrayLocation(("0", "0"))
    return


@pytest.fixture(scope="function")
def test_antenna_array():
    bm = aipy.phs.Beam(np.arange(0, 1, 0.1))
    a1 = aipy.phs.Antenna(0, 0, 0, bm, [0, 0])
    a2 = aipy.phs.Antenna(1, 0, 0, bm, [0, 1])
    a3 = aipy.phs.Antenna(0, 1, 0, bm, [0, 2])
    a4 = aipy.phs.Antenna(0, 0, 1, bm, [0, 3])
    ants = [a1, a2, a3, a4]
    aa = aipy.phs.AntennaArray(("0", "0"), ants)
    yield ants, aa
    return


def test_pointing_error():
    pnterr = aipy.phs.PointingError("Error String")
    assert str(pnterr) == "Error String"
    with pytest.raises(aipy.phs.PointingError):
        raise pnterr
    return


def test_ephem_zero():
    """Test converting between ephem dates and JD - zero point"""
    ephemzero = ephem.date("1899/12/31 12:00")
    jdzero = 2415020.0
    assert np.isclose(aipy.phs.juldate2ephem(jdzero), ephemzero)
    assert np.isclose(aipy.phs.ephem2juldate(ephemzero), jdzero)
    return


def test_ephem_random():
    """Test converting between ephem dates and JD - various"""
    for i in range(10):
        d1 = random.random() * ephem.now()
        d2 = aipy.phs.juldate2ephem(aipy.phs.ephem2juldate(d1))
        assert np.isclose(d1, d2)

    return


def test_radiobody_attributes():
    """Test aipy.phs.RadioFixedBody attributes"""
    epoch = ephem.B1950
    source = aipy.phs.RadioFixedBody(
        "0:00",
        "0:00",
        mfreq=0.200,
        name="src1",
        epoch=epoch,
        ionref=(0.0, 0.0),
        srcshape=(0.003, 0.005, 0.6),
    )
    assert source._ra == ephem.hours("0:00")
    assert source._dec == ephem.degrees("0:00")
    assert np.isclose(source.mfreq, 0.200)
    assert source.src_name == "src1"
    assert np.isclose(source._epoch, epoch, atol=1e-3)
    assert np.allclose(source.ionref, np.array([0.0, 0.0]))
    assert source.srcshape == [0.003, 0.005, 0.6]

    with pytest.raises(RuntimeError):
        source.ra
    with pytest.raises(RuntimeError):
        source.dec
    with pytest.raises(AttributeError):
        source.map

    obs = aipy.phs.ArrayLocation(("0:00", "0:00"))
    obs.set_ephemtime(epoch)
    source.compute(obs)
    assert len(source.get_crds("eq", ncrd=2)) == 2
    assert len(source.get_crds("top", ncrd=2)) == 2
    assert len(source.get_crds("eq", ncrd=3)) == 3
    assert len(source.get_crds("top", ncrd=3)) == 3
    assert source.map.shape == (3, 3)
    return


def test_ephem_interface():
    """Test the aipy.phs.RadioFixedBody ephem interface"""
    obs = aipy.phs.ArrayLocation(("0:00", "0:00"))
    for epoch in [ephem.B1900, ephem.B1950, ephem.J2000]:
        for ra in np.arange(np.pi / 8, 2 * np.pi, np.pi / 8):
            for dec in np.arange(-3 * np.pi / 8, np.pi / 2, np.pi / 8):
                source = aipy.phs.RadioFixedBody(ra, dec, epoch=epoch)
                obs.set_ephemtime(epoch)
                source.compute(obs)
                assert np.isclose(source.a_ra, source._ra, atol=1e-9)
                assert np.isclose(source.a_dec, source._dec, atol=1e-9)
    return


def test_ephem_compute():
    """Test the aipy.phs.RadioFixedBody ephem calculations"""
    epoch = ephem.J2000
    obs = aipy.phs.ArrayLocation(("0:00", "0:00"))
    obs.set_ephemtime(epoch)
    diagonal = np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])
    for ra in np.arange(np.pi / 8, 2 * np.pi, np.pi / 8):
        for dec in np.arange(-3 * np.pi / 8, np.pi / 2, np.pi / 8):
            source = aipy.phs.RadioFixedBody(ra, dec, epoch=epoch)
            source.compute(obs)
            mat = aipy.coord.top2eq_m(obs.sidereal_time() - source.ra, source.dec)
            err = np.abs(diagonal - np.dot(mat, source.map)).sum()
            assert np.isclose(err, 0.0, atol=1e-10)
    return


def test_get_crds(self):
    """Test the aipy.phs.RadioFixedBody calculated coordinates"""
    epoch = ephem.J2000
    obs = aipy.phs.ArrayLocation(("0:00", "0:00"))
    obs.set_ephemtime(epoch)
    for ra in np.arange(np.pi / 8, 2 * np.pi, np.pi / 8):
        for dec in np.arange(-3 * np.pi / 8, np.pi / 2, np.pi / 8):
            source = aipy.phs.RadioFixedBody(ra, dec, epoch=epoch)
            source.compute(obs)
            ra, dec = source.get_crds("eq", ncrd=2)
            assert np.isclose(source.ra, ra)
            assert np.isclose(source.dec, dec)

            az, alt = source.get_crds("top", ncrd=2)
            assert np.isclose(source.az, az)
            assert np.isclose(source.alt, alt)

            eq = source.get_crds("eq", ncrd=3)
            ra, dec = aipy.coord.eq2radec(eq)
            assert np.isclose(source.ra, ra, atol=1e-10)
            assert np.isclose(source.dec, dec, atol=1e-10)

            top = source.get_crds("top", ncrd=3)
            az, alt = aipy.coord.top2azalt(top)
            assert np.isclose(source.az, az, atol=1e-10)
            assert np.isclose(source.alt, alt, atol=1e-10)
    return


def test_add_srcs(source_catalog):
    """Test adding sources to a aipy.phs.SrcCatalog() catalog"""
    srcs, cat = source_catalog
    cat2 = aipy.phs.SrcCatalog()
    src1b = aipy.phs.RadioFixedBody("0:00", "0:00", name="src1")
    src4 = aipy.phs.RadioFixedBody("4:00", "4:00", name="src4")
    cat2.add_srcs(srcs)
    assert len(cat2) == 3
    cat2.add_srcs(src1b)
    assert len(cat2) == 3
    cat2.add_srcs(src1b, src4)
    assert len(cat2) == 4

    srclist = [src1b] + srcs[1:] + [src4]
    for name, src in zip([s.src_name for s in srclist], srclist):
        assert src == cat2[name]
    return


def test_get_srcs(source_catalog):
    """Test retrieving sources from a aipy.phs.SrcCatalog() catalog"""
    srcs, cat = source_catalog
    assert cat.get_srcs("src1", "src2", "src3") == srcs
    assert cat.get_srcs(["src1", "src2"]) == srcs[:2]
    with pytest.raises(KeyError):
        cat.get_srcs("bad")
    return


def test_compute(source_catalog):
    """Test the ephem interfaces for a aipy.phs.SrcCatalog() catalog"""
    srcs, cat = source_catalog
    obs = ephem.Observer()
    cat.compute(obs)
    for src in cat.values():
        assert src.ra is not None
        assert src.dec is not None
    return


def test_get_crds(source_catalog):
    """Test coordinates calculated from a aipy.phs.SrcCatalog() catalog"""
    srcs, cat = source_catalog
    obs = ephem.Observer()
    cat.compute(obs)
    crd1 = cat.get_crds("eq", srcs=["src1"])
    assert crd1.shape == (3, 1)
    assert np.allclose(crd1[:, 0], srcs[0].get_crds("eq"))

    crd2 = cat.get_crds("top", srcs=["src1", "src2"])
    assert crd2.shape == (3, 2)
    return


def test_get(source_catalog):
    """Test retrieving source attributes from a aipy.phs.SrcCatalog() catalog"""
    srcs, cat = source_catalog
    mfreq = cat.get("mfreq", srcs=["src1"])
    assert mfreq.shape == (1,)

    mfreq = cat.get("mfreq", srcs=["src1", "src2"])
    assert mfreq.shape == (2,)

    ionrefs = cat.get("ionref", srcs=["src1", "src2"])
    assert ionrefs.shape == (2, 2)

    srcshapes = cat.get("srcshape", srcs=["src1", "src2"])
    assert srcshapes.shape == (3, 2)
    return


def test_beam_attributes(test_beam):
    """Test accessing aipy.phs.Beam attributes"""
    freqs, bm = test_beam
    assert np.allclose(bm.freqs, freqs)
    assert np.allclose(bm.chans, np.arange(freqs.size))
    assert np.allclose(bm.afreqs, freqs)
    return


def test_beam_select_chans(test_beam):
    """Test selecting various aipy.phs.Beam channels"""
    freqs, bm = test_beam
    chans = np.array([1, 2, 3])
    bm.select_chans(chans)
    assert np.allclose(bm.chans, chans)
    assert np.allclose(bm.afreqs, freqs.take(chans))
    return


def test_antenna_attributes(test_antenna):
    """Test accessing aipy.phs.Antenna attributes"""
    freqs, bm, ant = test_antenna
    assert ant.beam == bm
    pos = np.array([1, 2, 3], dtype=np.float64)
    assert np.allclose(pos, ant.pos)
    x, y, z = ant
    assert np.isclose(x, pos[0])
    assert np.isclose(y, pos[1])
    assert np.isclose(z, pos[2])
    assert np.allclose(ant + ant, pos * 2)
    assert np.allclose(ant - ant, 0.0)
    assert np.allclose(-ant, -pos)
    return


def test_antenna_select_chans(test_antenna):
    """Test selecting various aipy.phs.Antenna channels"""
    freqs, bm, ant = test_antenna
    chans = np.array([1, 2, 3])
    ant.select_chans(chans)
    assert np.allclose(bm.chans, chans)
    assert np.allclose(ant.phsoff, 1)
    return


def test_array_location_attributes(test_array_location):
    al = test_array_location
    assert al.pressure == 0.0
    assert al.lat == 0.0
    assert al.long == 0.0
    assert al.elev == 0.0
    assert np.allclose(al._eq2zen, aipy.coord.eq2top_m(0.0, 0.0))
    return


def test_array_location_set_jultime(test_array_location):
    al = test_array_location
    jd = 2454554.9798841
    al.set_jultime(jd)
    assert al.date == al.epoch
    assert al.date == aipy.phs.juldate2ephem(jd)
    assert np.isclose(al.sidereal_time(), 0, atol=1e-7)

    eq2now_rnd = np.round(al._eq2now, 7)
    assert np.allclose(
        eq2now_rnd,
        np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    return


def test_array_location_get_jultime(test_array_location):
    al = test_array_location
    al.set_jultime(2454555)
    assert al.get_jultime() == 2454555
    return


def test_antenna_array_attributes(test_antenna_array):
    ants, aa = test_antenna_array
    assert len(aa) == 4
    for ai, aj in zip(aa, ants):
        assert ai == aj
    for i in range(4):
        assert len(aa[:i]) == i
    return


def test_antenna_array_select_chans(test_antenna_array):
    ants, aa = test_antenna_array
    chans = np.array([1, 2, 3])
    aa.select_chans(chans)
    for ant in aa.ants:
        assert np.allclose(ant.beam.chans, chans)
    return


def test_antenna_array_ij2bl(test_antenna_array):
    ants, aa = test_antenna_array
    assert aa.ij2bl(0, 1) == 258
    assert aa.ij2bl(0, 2) == 259
    assert aa.ij2bl(1, 2) == 515
    return


def test_antenna_array_bl2ij(test_antenna_array):
    ants, aa = test_antenna_array
    assert aa.bl2ij(258) == (0, 1)
    assert aa.bl2ij(259) == (0, 2)
    assert aa.bl2ij(515) == (1, 2)
    return


def test_antenna_array_get_baseline(test_antenna_array):
    ants, aa = test_antenna_array
    for j, ant in enumerate(aa):
        if j in [0, 3]:
            aa.set_jultime(2454554.9)
            assert np.allclose(aa.get_baseline(0, j, "r"), ant.pos)
            assert np.allclose(aa.get_baseline(0, j, "e"), ant.pos)
        else:
            assert np.allclose(aa.get_baseline(0, j, "r"), ant.pos)
            aa.set_jultime(2454554.9798841)
            bl_rnd = np.round(aa.get_baseline(0, j, "e"), 7)
            assert np.allclose(bl_rnd, ant.pos)

            aa.set_jultime(2454554.9)
            bl_rnd = np.round(aa.get_baseline(0, j, "e"), 7)
            assert not np.allclose(bl_rnd, ant.pos)
    src = aipy.phs.RadioFixedBody("12:00", "0:00")
    src.compute(aa)
    with pytest.raises(aipy.phs.PointingError):
        aa.get_baseline(0, 1, src)
    for t in np.random.random((10,)):
        aa.set_jultime(2454555.0 + t)
        src = aipy.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, epoch=aa.epoch)
        src.compute(aa)
        zbl_rnd = np.round(aa.get_baseline(0, 1, "z"), 3)
        sbl_rnd = np.round(aa.get_baseline(0, 1, src), 3)
        assert np.allclose(zbl_rnd, sbl_rnd)
    return


def test_antenna_array_get_phs_offset(test_antenna_array):
    ants, aa = test_antenna_array
    assert np.allclose(aa.get_phs_offset(0, 0), 0)
    assert np.allclose(aa.get_phs_offset(0, 1), 1)
    assert np.allclose(aa.get_phs_offset(0, 2), 2)
    assert np.allclose(aa.get_phs_offset(0, 3), 3)
    return


def test_antenna_array_gen_uvw(test_antenna_array):
    ants, aa = test_antenna_array
    aa.select_chans()
    afreqs = aa[0].beam.afreqs
    u, v, w = aa.gen_uvw(0, 1, "z")
    assert u.shape == (1, afreqs.size)
    assert np.allclose(u, 0 * afreqs)
    assert np.allclose(v, 0 * afreqs)
    assert np.allclose(w, 1 * afreqs)

    u, v, w = aa.gen_uvw(0, 2, "z")
    assert np.allclose(u, 1 * afreqs)
    assert np.allclose(v, 0 * afreqs)
    assert np.allclose(w, 0 * afreqs)

    u, v, w = aa.gen_uvw(0, 3, "z")
    assert np.allclose(u, 0 * afreqs)
    assert np.allclose(v, 1 * afreqs)
    assert np.allclose(w, 0 * afreqs)
    return


def test_antenna_array_gen_phs(test_antenna_array):
    ants, aa = test_antenna_array
    aa.select_chans([1, 2, 3])
    afreqs = aa[0].beam.afreqs
    for t in np.random.random((40,)):
        aa.set_jultime(2454555.0 + t)
        src = aipy.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, epoch=aa.epoch)
        src.compute(aa)
        if t > 0.5:
            seq = src.get_crds("eq", ncrd=3)
            if t > 0.75:
                seq = np.array([seq, seq]).transpose()
        else:
            seq = src
        phs = np.round(aa.gen_phs(seq, 0, 1, mfreq=0.1), 6)
        ans = np.round(np.exp(-1j * 2 * np.pi * afreqs), 6)
        if t > 0.75:
            assert phs.shape == (2, 3)
        else:
            assert phs.shape == (3,)
        assert np.allclose(phs, ans)
        phs = np.round(aa.gen_phs(src, 0, 2, mfreq=0.1), 3)
        assert np.allclose(phs, 1 + 0j)
        phs = np.round(aa.gen_phs(src, 0, 3, mfreq=0.1), 3)
        assert np.allclose(phs, 1 + 0j)

    phs1 = aa.gen_phs(src, 0, 2, mfreq=0.1, ionref=(0.001, 0.001))
    phs2 = aa.gen_phs(src, 0, 2, mfreq=0.1, srcshape=(0.01, 0.01, 0), resolve_src=True)
    assert np.all(phs1 != 1 + 0j)
    assert np.all(phs2 != 1 + 0j)
    return


def test_antenna_array_resolve_src(test_antenna_array):
    ants, aa = test_antenna_array
    # we get a runtime warning when x ~ 0
    with pytest.warns(RuntimeWarning) as record:
        amp = aa.resolve_src(100.0, 100.0, srcshape=(0, 0, 0))
    assert len(record) == 1
    assert record[0].message.args[0] == "invalid value encountered in double_scalars"

    # check that value is correct
    assert amp == 1

    amp1 = aa.resolve_src(100.0, 50.0, srcshape=(0.01, 0, np.pi / 2))
    amp2 = aa.resolve_src(100.0, 50.0, srcshape=(0, 0.01, 0))
    assert np.isclose(amp1, amp2, atol=1e-15)

    amp1 = aa.resolve_src(100.0, 100.0, srcshape=(0.02, 0, 0))
    amp2 = aa.resolve_src(100.0, 100.0, srcshape=(0.01414, 0.01414, 0))
    assert np.isclose(amp1, amp2, atol=1e-3)

    amp = aa.resolve_src(100.0, 0.0, srcshape=(0.001, 0, 0))
    x = 2 * np.pi * 0.1
    assert np.isclose(amp, 2 * aipy.phs.j1(x) / x)
    return


def test_refract(test_antenna_array):
    ants, aa = test_antenna_array
    aa.select_chans([1, 2, 3])
    afreqs = aa[0].beam.afreqs
    zeros = np.zeros((1, 3), dtype=np.float64)
    ones = np.ones((1, 3), dtype=np.float64)
    # Test non-vectors, dra->u association
    dw = aa.refract(ones, zeros, mfreq=0.1, ionref=(0.001, 0))
    ans = 0.001 / (afreqs / 0.1) ** 2
    ans.shape = (1, ans.size)
    assert len(dw.shape) == 2
    assert np.allclose(np.round(dw, 10), np.round(ans, 10))

    # Test non-vectors, no dra-> v association
    dw = aa.refract(zeros, ones, mfreq=0.1, ionref=(0.001, 0))
    assert np.allclose(dw, 0)

    # Test non-vectors, no ddec->u association
    dw = aa.refract(ones, zeros, mfreq=0.1, ionref=(0, 0.001))
    assert np.allclose(dw, 0)

    # Test non-vectors, ddec->v association, v scaling
    dw = aa.refract(zeros, 2 * ones, mfreq=0.1, ionref=(0, 0.001))
    assert np.allclose(np.round(dw, 10), np.round(2 * ans, 10))

    # Test vectors, mfreq scaling
    ones = np.ones((2, 3), dtype=np.float64)
    ionref = (np.array([0, 0.001]), np.array([0.001, 0]))
    mfreq = np.array([0.1, 0.2])
    ans = np.array([0.002 / (afreqs / 0.1) ** 2, 0.002 / (afreqs / 0.2) ** 2])
    dw = aa.refract(2 * ones, 2 * ones, mfreq=mfreq, ionref=ionref)
    assert np.allclose(np.round(dw, 10), np.round(ans, 10))
    return


def test_antenna_array_phs2src(test_antenna_array):
    ants, aa = test_antenna_array
    aa.select_chans([1, 2, 3])
    aa.set_jultime(2454555.0)
    src = aipy.phs.RadioFixedBody("0:00", "20:00")
    src.compute(aa)
    assert np.allclose(aa.phs2src(1.0, src, 0, 1), aa.gen_phs(src, 0, 1))
    return


def test_unphs2src(test_antenna_array):
    ants, aa = test_antenna_array
    aa.select_chans([1, 2, 3])
    aa.set_jultime(2454555.0)
    src = aipy.phs.RadioFixedBody("0:00", "20:00")
    src.compute(aa)
    assert np.allclose(aa.unphs2src(aa.gen_phs(src, 0, 1), src, 0, 1), 1.0)
    return
