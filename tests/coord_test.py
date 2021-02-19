# -*- coding: utf-8 -*-
# Copyright (c) 2009 Aaron Parsons
# Licensed under the GPLv3

import pytest
import ephem
import random
import aipy
import numpy as np

def test_convert_typecheck():
    """Test coordinate conversion types"""
    crd2, crd3 = (0, 0), (1, 0, 0)
    for sys in aipy.coord.sys_dict:
        ans2 = aipy.coord.convert(crd2, sys, sys)
        assert crd2[0] == ans2[0]
        assert crd2[1] == ans2[1]

        ans3 = aipy.coord.convert(crd3, sys, sys)
        assert crd2[0] == ans3[0]
        assert crd2[1] == ans3[1]

        with pytest.raises(KeyError):
            aipy.coord.convert(crd2, "bad", "bad")

    return

def test_convert_precession():
    """Test coordinate precessions for accuracy"""
    crdpairs = [
        [('0:00','0:00'), ('00:02:33.77','00:16:42.1')],
        [('6:00','0:00'), ('06:02:33.75','-00:00:05.6')],
        [('0:00','89:00'), ('00:03:03.75','+89:16:41.7')],
    ]
    for b1950, j2000 in crdpairs:
        c1 = aipy.coord.convert(
            b1950, 'eq', 'eq', iepoch=ephem.B1950, oepoch=ephem.J2000
        )
        c2 = aipy.coord.convert(
            j2000, 'eq', 'eq', iepoch=ephem.J2000, oepoch=ephem.B1950
        )
        c1_ck = ephem.Equatorial(j2000[0], j2000[1], epoch=ephem.J2000).get()
        c2_ck = ephem.Equatorial(b1950[0], b1950[1], epoch=ephem.B1950).get()
        assert np.allclose(c1[0], c1_ck[0], atol=1e-4)
        assert np.allclose(c1[1], c1_ck[1], atol=1e-4)
        assert np.allclose(c2[0], c2_ck[0], atol=1e-4)
        assert np.allclose(c2[1], c2_ck[1], atol=1e-4)

    return

def test_convert_crdsys():
    """Test coordinate conversions for accuracy"""
    eq_ec_ga = [
        [
            ('19:59:28.3566', '40:44:02.096'),
            ('317.7054323', '59.3254895'),
            ('76.1898379', '5.7554756'),
        ],
        [
            ('12:30:49.4233', '12:23:28.043'),
            ('182.0592608', '14.4166861'),
            ('283.7777978', '74.4911308'),
        ],
        [
            ('13:25:27.6152', '-43:01:08.805'),
            ('217.1433477', '-31.3319020'),
            ('309.5158743', '19.4173247'),
        ]
    ]
    for eq, ec, ga in eq_ec_ga:
        eq_ec = aipy.coord.convert(eq,'eq','ec')
        eq_ga = aipy.coord.convert(eq,'eq','ga')
        ec_eq = aipy.coord.convert(ec,'ec','eq')
        ec_ga = aipy.coord.convert(ec,'ec','ga')
        ga_eq = aipy.coord.convert(ga,'ga','eq')
        ga_ec = aipy.coord.convert(ga,'ga','ec')
        eq_ck = ephem.Equatorial(eq[0], eq[1], epoch=ephem.J2000).get()
        ec_ck = ephem.Ecliptic(ec[0], ec[1], epoch=ephem.J2000).get()
        ga_ck = ephem.Galactic(ga[0], ga[1], epoch=ephem.J2000).get()

        assert np.allclose(eq_ec[0], ec_ck[0], atol=1e-4)
        assert np.allclose(eq_ec[1], ec_ck[1], atol=1e-4)
        assert np.allclose(eq_ga[0], ga_ck[0], atol=1e-4)
        assert np.allclose(eq_ga[1], ga_ck[1], atol=1e-4)
        assert np.allclose(ec_eq[0], eq_ck[0], atol=1e-4)
        assert np.allclose(ec_eq[1], eq_ck[1], atol=1e-4)
        assert np.allclose(ec_ga[0], ga_ck[0], atol=1e-4)
        assert np.allclose(ec_ga[1], ga_ck[1], atol=1e-4)
        assert np.allclose(ga_eq[0], eq_ck[0], atol=1e-4)
        assert np.allclose(ga_eq[1], eq_ck[1], atol=1e-4)
        assert np.allclose(ga_ec[0], ec_ck[0], atol=1e-4)
        assert np.allclose(ga_ec[1], ec_ck[1], atol=1e-4)

    return

def test_convert_m_shape():
    """Test conversion matrices for shape preservation"""
    for c1 in aipy.coord.sys_dict:
        for c2 in aipy.coord.sys_dict:
            mat = aipy.coord.convert_m(c1,c2)
            assert len(mat.shape) == 2
            assert mat.shape[0] == 3
            assert mat.shape[1] == 3

    return

def test_convert_m_diagonal():
    """Test conversion matrices for normalcy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    for c1 in aipy.coord.sys_dict:
        for c2 in aipy.coord.sys_dict:
            m1 = aipy.coord.convert_m(c1, c2)
            m2 = aipy.coord.convert_m(c2, c1)
            mat = np.round(np.dot(m1, m2), 10)
            assert np.allclose(mat, diag)

    return

def test_rot_m_shape():
    """Test rotation matrices for shape preservation"""
    mat = aipy.coord.rot_m(0, np.array([1., 0, 0]))
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 3

    mat = aipy.coord.rot_m(np.array([0., 0.]), np.array([[1., 0, 0], [1., 0, 0]]))
    assert len(mat.shape) == 3
    assert mat.shape[0] == 2
    assert mat.shape[1] == 3
    assert mat.shape[2] == 3

    return

def test_rot_m_accuracy():
    """Test rotation matricies for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.rot_m(2 * np.pi, np.array([1., 0, 0])), 10)
    assert np.allclose(mat, diag)

    mat = np.round(aipy.coord.rot_m(np.pi/2, np.array([0.,0,1])), 10)
    assert np.allclose(mat, np.array([-e2, e1, e3]))

    mat = np.round(aipy.coord.rot_m(np.pi / 2, np.array([0, 1, 0])), 10)
    assert np.allclose(mat, np.array([e3, e2, -e1]))

    mat = np.round(aipy.coord.rot_m(np.pi / 2, np.array([1, 0, 0])), 10)
    assert np.allclose(mat, np.array([e1, -e3, e2]))

    return

def test_xyz2thphi_shape():
    """Test the x,y,z to theta,phi conversion for shape preservation"""
    mat = aipy.coord.xyz2thphi((0, 0, 1))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 2

    mat = aipy.coord.xyz2thphi(
        (np.array([0., 0, 0]), np.array([0., 0, 0]), np.array([1., 1, 1]))
    )

    assert len(mat.shape) == 2
    assert mat.shape[0] == 2
    assert mat.shape[1] == 3

    return

def test_xyz2thphi_accuracy():
    """Test the x,y,z to theta,phi conversion for accuracy"""
    diag = np.array([[1,0,0],[0,1,0],[0,0,1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.xyz2thphi(e3), 10)
    assert np.allclose(mat, np.array([0.,0]))

    mat = np.round(aipy.coord.xyz2thphi(e1), 10)
    assert np.allclose(mat, np.round(np.array([np.pi / 2, 0]), 10))

    mat = np.round(aipy.coord.xyz2thphi(e2), 10)
    assert np.allclose(mat, np.round(np.array([np.pi / 2, np.pi / 2]), 10))

    return

def test_thphi2xyz_shape():
    """Test the theta,phi to x,y,z conversion for shape preservation"""
    mat = aipy.coord.thphi2xyz((0, 0))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 3

    mat = aipy.coord.thphi2xyz((np.array([0., 0]), np.array([0., 0])))
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 2

    return

def test_thphi2xyz_accuracy():
    """Test the theta,phi to x,y,z conversion for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.thphi2xyz((0, 0)), 10)
    assert np.allclose(mat, e3)

    mat = np.round(aipy.coord.thphi2xyz((np.pi / 2, 0)), 10)
    assert np.allclose(mat, e1)

    mat = np.round(aipy.coord.thphi2xyz((np.pi / 2, np.pi / 2)), 10)
    assert np.allclose(mat, e2)

    return

def test_eq2radec_shape():
    """Test the equatorial to ra,dec conversion for shape preservation"""
    mat = aipy.coord.eq2radec((0, 0, 1))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 2

    mat = aipy.coord.eq2radec(
        (np.array([0., 0, 0]), np.array([0., 0, 0]), np.array([1., 1, 1]))
    )
    assert len(mat.shape) == 2
    assert mat.shape[0] == 2
    assert mat.shape[1] == 3

    return

def test_eq2radec_accuracy():
    """Test the equatorial to ra,dec conversion for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.eq2radec(e3), 10)
    assert np.allclose(mat, np.round(np.array([0., np.pi / 2]), 10))

    mat = np.round(aipy.coord.eq2radec(e1), 10)
    assert np.allclose(mat, np.round(np.array([0, 0]), 10))

    mat = np.round(aipy.coord.eq2radec(e2), 10)
    assert np.allclose(mat, np.round(np.array([np.pi / 2,0]), 10))

    return

def test_radec2eq_shape():
    """Test the ra,dec to equatorial conversion for shape preservation"""
    mat = aipy.coord.radec2eq((0, 0))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 3


    mat = aipy.coord.radec2eq((np.array([0., 0]), np.array([0., 0])))
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 2

    return

def test_radec2eq_accuracy():
    """Test the ra,dec to equatorial conversion for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.radec2eq((0, 0)), 10)
    assert np.allclose(mat, e1)

    mat = np.round(aipy.coord.radec2eq((np.pi / 2, 0)), 10)
    assert np.allclose(mat, e2)

    mat = np.round(aipy.coord.radec2eq((np.pi / 2, np.pi / 2)), 10)
    assert np.allclose(mat, e3)

    return

def test_latlong2xyz_shape():
    """Test the lat,long to x,y,z conversion for shape preservation"""
    mat = aipy.coord.latlong2xyz((0, 0))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 3

    mat = aipy.coord.latlong2xyz((np.array([0., 0]), np.array([0., 0])))
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 2

    return

def test_latlong2xyz_accuracy():
    """Test the lat,long to x,y,z conversion for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.latlong2xyz((0, 0)), 10)
    assert np.allclose(mat, e1)

    mat = np.round(aipy.coord.latlong2xyz((np.pi / 2, 0)), 10)
    assert np.allclose(mat, e3)

    mat = np.round(aipy.coord.latlong2xyz((0, np.pi / 2)), 10)
    assert np.allclose(mat, e2)

    return

def test_top2azalt_shape():
    """Test the x,y,z to az,alt conversion for shape preservation"""
    mat = aipy.coord.top2azalt((0, 0, 1))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 2

    mat = aipy.coord.top2azalt(
        (np.array([0., 0, 0]), np.array([0., 0, 0]), np.array([1., 1, 1]))
    )
    assert len(mat.shape) == 2
    assert mat.shape[0] == 2
    assert mat.shape[1] == 3

    return

def test_top2azalt_accuracy():
    """Test the x,y,z to az,alt conversion for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = aipy.coord.top2azalt(e3)
    assert np.allclose(mat[1], np.pi / 2, atol=1e-10)

    mat = np.round(aipy.coord.top2azalt(e1), 10)
    assert np.allclose(mat, np.round(np.array([np.pi / 2, 0]), 10))

    mat = np.round(aipy.coord.top2azalt(e2), 10)
    assert np.allclose(mat, np.round(np.array([0,0]), 10))

    return

def test_azalt2top_shape():
    """Test the az,alt to x,y,z conversion for shape preservation"""
    mat = aipy.coord.azalt2top((0, 0))
    assert len(mat.shape) == 1
    assert mat.shape[0] == 3

    mat = aipy.coord.azalt2top((np.array([0., 0]), np.array([0., 0])))
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 2

    return

def test_azalt2top_accuracy():
    """Test the az,alt to x,y,z conversion for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.azalt2top((0, 0)), 10)
    assert np.allclose(mat, e2)

    mat = np.round(aipy.coord.azalt2top((np.pi / 2, 0)), 10)
    assert np.allclose(mat, e1)

    mat = np.round(aipy.coord.azalt2top((0, np.pi / 2)), 10)
    assert np.allclose(mat, e3)

    return

def test_eq2top_m_shape():
    """Test the equatorial/x,y,z rotation matrix for shape preservation"""
    mat = aipy.coord.eq2top_m(0, 0)
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 3

    mat = aipy.coord.eq2top_m(np.array([0., 0]), np.array([0., 0]))
    assert len(mat.shape) == 3
    assert mat.shape[0] == 2
    assert mat.shape[1] == 3
    assert mat.shape[2] == 3

    return

def test_eq2top_m_accuracy():
    """Test the equatorial/x,y,z rotation matrix for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.eq2top_m(0, 0), 10)
    assert np.allclose(mat, np.array([e2, e3, e1]))

    mat = np.round(aipy.coord.eq2top_m(-np.pi / 2, 0.), 10)
    assert np.allclose(mat, np.array([-e1, e3, e2]))

    mat = np.round(aipy.coord.eq2top_m(0, np.pi / 2), 10)
    assert np.allclose(mat, np.array([e2, -e1, e3]))

    return

def test_top2eq_m_shape():
    """Test the x,y,z/equatorial rotation matrix for shape preservation"""
    mat = aipy.coord.top2eq_m(0, 0)
    assert len(mat.shape) == 2
    assert mat.shape[0] == 3
    assert mat.shape[1] == 3

    mat = aipy.coord.top2eq_m(np.array([0., 0]), np.array([0.,0]))
    assert len(mat.shape) == 3
    assert mat.shape[0] == 2
    assert mat.shape[1] == 3
    assert mat.shape[2] == 3

    return

def test_top2eq_m_accuracy():
    """Test the x,y,z/equatorial rotation matrix for accuracy"""
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    e1, e2, e3 = diag
    mat = np.round(aipy.coord.top2eq_m(0, 0), 10)
    assert np.allclose(mat, np.array([e3, e1, e2]))

    mat = np.round(aipy.coord.top2eq_m(-np.pi / 2, 0.), 10)
    assert np.allclose(mat, np.array([-e1, e3, e2]))

    mat = np.round(aipy.coord.top2eq_m(0, np.pi / 2), 10)
    assert np.allclose(mat, np.array([-e2, e1, e3]))

    return
