# -*- coding: utf-8 -*-
# Copyright (c) 2018 Aaron Parsons
# Licensed under the GPLv3

import aipy.dsp as dsp
import numpy as np

def test_kaiser2():
    win = dsp.gen_window(1024, 'kaiser2')
    assert np.allclose(win[0], 0.01147993)
    assert np.allclose(win[1], 0.01192681)
    assert np.allclose(win[2], 0.01238142)
