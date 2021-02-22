# -*- coding: utf-8 -*-
# Copyright (c) 2009 Aaron Parsons
# Licensed under the GPLv3

import re
import pytest
import numpy as np

import aipy
from aipy.miriad import ij2bl


@pytest.fixture(scope="function")
def chan_cases():
    nchan = 256
    cases = {
        "all": [np.arange(nchan)],
        "0,1,2,3": [np.array([0]), np.array([1]), np.array([2]), np.array([3])],
        "0_20": [np.arange(21)],
        "0_20_4": [np.arange(0, 21, 4)],
        "5,7,40_50": [np.array([5]), np.array([7]), np.arange(40, 51)],
    }
    yield nchan, cases
    return


@pytest.fixture(scope="function")
def test_strs():
    yield [
        "src1=jys",
        "src2=jys/100",
        "src3=jys/100/1",
        "src4=jys//1",
        "(src5/src6)=jys",
        "(src7/src8)=jys/(100/200)",
        "(src9/src10)=jys/(100/200)/(1/2)",
        "(src11/src12)=jys/(100/200)/1",
        "(src13/src14/src15)=(jys/index)//1",
    ]
    return


def test_parse_ants():
    """Test aipy.scripting.parse_ants()"""
    nants = 4
    cases = {
        "all": [],
        "auto": [("auto", 1)],
        "cross": [("auto", 0)],
        "0_1": [(ij2bl(0, 1), 1)],
        "0_1,1_2": [(ij2bl(0, 1), 1), (ij2bl(1, 2), 1)],
        "0x_1x": [(ij2bl(0, 1), 1, "xx")],
        "(0x,0y)_1x": [(ij2bl(0, 1), 1, "xx"), (ij2bl(0, 1), 1, "yx")],
        "(0,1)_2": [(ij2bl(0, 2), 1), (ij2bl(1, 2), 1)],
        "0_(1,2)": [(ij2bl(0, 1), 1), (ij2bl(0, 2), 1)],
        "(0,1)_(2,3)": [
            (ij2bl(0, 2), 1),
            (ij2bl(0, 3), 1),
            (ij2bl(1, 2), 1),
            (ij2bl(1, 3), 1),
        ],
        "0_(1,-2)": [(ij2bl(0, 1), 1), (ij2bl(0, 2), 0)],
        "(-0,1)_(2,-3)": [
            (ij2bl(0, 2), 0),
            (ij2bl(0, 3), 0),
            (ij2bl(1, 2), 1),
            (ij2bl(1, 3), 0),
        ],
        "0,1,all": [],
    }
    for i in range(nants):
        cases[str(i)] = list(map(lambda x: (ij2bl(x, i), 1), range(nants)))
        cases["-" + str(i)] = list(map(lambda x: (ij2bl(x, i), 0), range(nants)))
    # inelegantly paste on the new pol parsing flag on the above tests
    # XXX really should add some new tests for the new pol parsing
    for k in cases:
        cases[k] = [(v + (-1,))[:3] for v in cases[k]]
    for ant_str in cases:
        assert aipy.scripting.parse_ants(ant_str, nants) == cases[ant_str]

    with pytest.raises(ValueError):
        aipy.scripting.parse_ants("(0_1)_2", nants)

    return


def test_parse_srcs():
    """Test aipy.scripting.parse_srcs()"""
    cat = "misc,helm"
    src_opt = "cyg,cas,12_40,12:00_40:00,12:00:00.0_40:00:00.0,bug_man"
    srclist, cutoff, catalogs = aipy.scripting.parse_srcs(src_opt, cat)
    assert srclist[0] == "cyg"
    assert srclist[1] == "cas"
    assert srclist[2].src_name == "12_40"
    assert srclist[3].src_name == "12:00_40:00"
    assert srclist[4].src_name == "12:00:00.0_40:00:00.0"
    assert srclist[5] == "bug_man"
    return


def test_parse_chans_concat(chan_cases):
    """Test aipy.scripting.parse_chans() - concatenate"""
    nchan, cases = chan_cases
    for case in cases:
        chans = aipy.scripting.parse_chans(case, nchan, concat=True)
        assert np.allclose(chans, np.concatenate(cases[case]))

    with pytest.raises(AssertionError):
        aipy.scripting.parse_chans("0_1_2_3", nchan)

    return


def test_parse_chans_noconcat(chan_cases):
    """Test aipy.scripting.parse_chans() - without concatenate"""
    nchan, cases = chan_cases
    for case in cases:
        chans = aipy.scripting.parse_chans(case, nchan, concat=False)
        for ch1, ch2 in zip(chans, cases[case]):
            assert np.allclose(ch1, ch2)
    return


def test_name():
    """Test aipy.scripting.name matching"""
    regex = re.compile(aipy.scripting.name)
    assert regex.match("src1").groups() == ("src1",)
    assert regex.match("1").groups() == ("1",)
    assert regex.match("1,1").groups() == ("1",)
    assert regex.match("1=1").groups() == ("1",)
    assert regex.match("1/1").groups() == ("1",)
    assert regex.match("1(1").groups() == ("1",)
    assert regex.match("1)1").groups() == ("1",)
    assert regex.match(")1") is None
    return


def test_grp():
    """Test aipy.scripting.grp matching"""
    regex = re.compile(aipy.scripting.grp)
    assert regex.match("src1").groups()[0] == "src1"
    assert regex.match("src1/src2").groups()[0] == "src1"
    assert regex.match("(src1/src2)").groups()[0] == "(src1/src2)"
    assert regex.match("(1/2/3)").groups()[0] == "(1/2/3)"
    assert regex.match("(src1/)") is None
    assert regex.match("(src1/src2") is None
    return


def test_prm_rgx(test_strs):
    """Test aipy.scripting.prm_rgx matching"""
    regex = aipy.scripting.prm_rgx
    for t in test_strs:
        assert regex.match(t).groups()[0] == t
    return


def test_parse_prms(test_strs):
    """Test aipy.scripting.parse_prms()"""
    t = ",".join(test_strs)
    prms = aipy.scripting.parse_prms(t)
    assert prms["src1"]["jys"] == (None, None)
    assert prms["src2"]["jys"] == (100.0, None)
    assert prms["src3"]["jys"] == (100.0, 1.0)
    assert prms["src4"]["jys"] == (None, 1.0)
    assert prms["src5"]["jys"] == (None, None)
    assert prms["src6"]["jys"] == (None, None)
    assert prms["src7"]["jys"] == (100.0, None)
    assert prms["src8"]["jys"] == (200.0, None)
    assert prms["src9"]["jys"] == (100.0, 1.0)
    assert prms["src10"]["jys"] == (200.0, 2.0)
    assert prms["src11"]["jys"] == (100.0, 1.0)
    assert prms["src12"]["jys"] == (200.0, 1.0)
    assert prms["src13"]["jys"] == (None, 1.0)
    assert prms["src13"]["index"] == (None, 1.0)
    assert prms["src14"]["jys"] == (None, 1.0)
    assert prms["src14"]["index"] == (None, 1.0)
    assert prms["src15"]["jys"] == (None, 1.0)
    assert prms["src15"]["index"] == (None, 1.0)

    with pytest.raises(AssertionError):
        aipy.scripting.parse_prms("(a/b)=(c/d)/(1/2)/(3/4)")

    t = "(1/2/3)=jys,(2/3)=index"
    prms = aipy.scripting.parse_prms(t)
    assert len(prms["1"]) == 1
    assert len(prms["2"]) == 2
    assert len(prms["3"]) == 2
    prms = aipy.scripting.parse_prms("a=(b/c)")
    assert len(prms["a"]) == 2
