#! /usr/bin/env python
"""
Apply the bandpass function in a UV file to the raw data, and then write
to a new UV file which will not have a bandpass file.  Has the (recommended)
option of linearizing for quantization gain effects.  For now, these are
applied only to the auto-correlations.  Cross-correlation quantization gain is
less sensitive to level, so linearization will come later.
Author: Aaron Parsons
Date: 8/14/07
Revisions:
    9/25/07 arp Uninverted bandpass and untransposed to sync with what
                Miriad expects from the bandpass variable.  Removes variables
                associated with bandpass (correctly).
    12/11/07 arp    Updated to new miriad file interface
    05/12/08 arp    Sped up reading w/ "raw" mode
"""

__version__ = '0.0.1'

import aipy as a, numpy as np, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('apply_bp.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-l', '--linearization', dest='linearization', default='comb',
    help='Apply the specified quantization linearization function to raw correlator values before applying bandpass.  Options are null, digi, full, and comb.  Default is comb')
o.add_option('-s', '--scale', dest='scale', type='float', default=12250000.,
    help='An additional numerical scaling to apply to the data.  Default: 12250000.')
opts, args = o.parse_args(sys.argv[1:])

# Digital Gain, Power Output (Channel 1024+512)
# NOD-5109 Noise Generator (Coarse=20dB, Fine=0dB attenuation)
# -3 dB, SHP-100, -3 dB, VLF-225
# port A pocket_correlator, shift=x3ff, clk=600MHz
# Data is [digital gain, power output]

qdata1 = np.array([
     [12157.543882402342, 69.423910140999993],
     [10941.789494162109, 66.790805816700001],
     [9847.6105447458976, 63.422136306799999],
     [8862.8494902713082, 59.848737716700001],
     [7976.5645412441772, 56.260694503800003],
     [7178.9080871197593, 52.3417243958],
     [6461.0172784077831, 48.231815338099999],
     [5814.9155505670051, 43.995954513500003],
     [5233.4239955103048, 39.623107910199998],
     [4710.0815959592746, 35.128669738799999],
     [4239.0734363633474, 30.719709396399999],
     [3815.1660927270127, 26.564506530799999],
     [3433.6494834543114, 22.680536270099999],
     [3090.2845351088804, 18.9764328003],
     [2781.2560815979923, 15.7225093842],
     [2503.1304734381933, 12.9534130096],
     [2252.8174260943742, 10.588439941400001],
     [2027.535683484937, 8.6283283233599999],
     [1824.7821151364433, 7.0549869537400003],
     [1642.303903622799, 5.7384414672900004],
     [1478.0735132605191, 4.6558799743700003],
     [1330.2661619344672, 3.8215141296400001],
     [1197.2395457410205, 3.1258945465100001],
     [1077.5155911669185, 2.5557441711400002],
     [969.76403205022666, 2.1025066375699999],
     [872.78762884520404, 1.72937583923],
     [785.50886596068369, 1.4391326904299999],
     [706.95797936461531, 1.19428634644],
     [636.26218142815378, 0.99624633789100003],
     [572.63596328533845, 0.83918571472199999],
     [515.37236695680463, 0.70255088806199995],
])

# Analog Gain (Coarse, Fine Attenuation), Power Output (Channel 1024+512)
# D-5109 Noise Generator (Coarse, Fine attenuation)
# -3 dB, SHP-100, -3 dB, VLF-225
# port A pocket_correlator, shift=x3ff, gain=[3500], clk=600MHz
# Data is [Attenuation (Coarse, Fine), Power Output]

qdata2 = np.array([
    [40,  4,  0.43  ],
    [40,  3,  0.46  ],
    [40,  2,  0.50  ],
    [40,  1,  0.55  ],
    [40,  0,  0.60  ],
    [30, 10,  0.60  ],
    [30,  9,  0.66  ],
    [30,  8,  0.74  ],
    [30,  7,  0.84  ],
    [30,  6,  1.0   ],
    [30,  5,  1.1   ],
    [30,  4,  1.3   ],
    [30,  3,  1.6   ],
    [30,  2,  1.9   ],
    [30,  1,  2.3   ],
    [30,  0,  2.7   ],
    [20, 10,  2.7   ],
    [20,  9,  3.4   ],
    [20,  8,  4.2   ],
    [20,  7,  5.2   ],
    [20,  6,  6.5   ],
    [20,  5,  8.1   ],
    [20,  4, 10.1   ],
    [20,  3, 12.6   ],
    [20,  2, 15.7   ],
    [20,  1, 19.3   ],
    [20,  0, 23.3   ],
    [10, 10, 23.7   ],
    [10,  9, 28.3   ],
    [10,  8, 32.9   ],
    [10,  7, 37.9   ],
    [10,  6, 42.6   ],
    [10,  5, 47.4   ],
    [10,  4, 51.9   ],
    [10,  3, 56.2   ],
    [10,  2, 60.2   ],
    [10,  1, 63.8   ],
    [10,  0, 66.9   ],
])

# Define some polynomials for interpreting between the datapoints above
# and use them to apply the correction function.  Polynomials are normalized
# to map a raw value of 10 to itself.

digi_cpoly = np.polyfit(qdata1[:,1], qdata1[:,0]**2, 6)
full_cpoly = np.polyfit(qdata2[:,2], 10**(-(qdata2[:,0]+qdata2[:,1])/10), 6)
digi_cpoly *= 10 / np.polyval(digi_cpoly, 10)
full_cpoly *= 10 / np.polyval(full_cpoly, 10)
null_cpoly = np.array([0, 0, 0, 0, 0, 1, 0])
# Unscientifically average the 3 linearization polynomials.  In practice, this
# seems to work best.
comb_cpoly = (null_cpoly + full_cpoly + digi_cpoly) / 3
cpolys = {'null':null_cpoly, 'digi':digi_cpoly, 
    'full':full_cpoly, 'comb':comb_cpoly}

cpoly = cpolys[opts.linearization]

print 'Using %s quantization correction' % opts.linearization

# These are all the items which should be removed once bandpass applied.
ignore_vars = ['bandpass', 'freqs', 'ngains', 'nspect0', 
    'nchan0', 'ntau', 'nfeeds', 'nsols', 'header', 'vartable']

# Process all files passed from the command line.
for filename in args:
    print filename,'->',filename+'b'
    if os.path.exists(filename+'b'):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'b', status='new')
    uvo.init_from_uv(uvi, exclude=ignore_vars)
    nchan = uvi['nchan']
    nants = uvi['nants']
    try:
        # If there is a bandpass item, we're going apply it to data
        bp = uvi['bandpass'].real   # Sync'd with pocket_corr.py in corr pkg.
        bp.shape = (nants, nchan)
        print
    except:
        print 'No bandpass found'
        bp = np.ones((nants, nchan))
        print '.'
    def f(uv, preamble, data, flags):
        uvw, t, (i,j) = preamble
        if i == j: data = np.polyval(cpoly, data)
        data *= bp[i,:] * bp[j,:] * opts.scale
        return preamble, data, flags
    uvo.pipe(uvi, mfunc=f, raw=True,
        append2hist='APPLY_BP: ver=%s, corr type=%s, scale=%f\n' % \
            (__version__, opts.linearization, opts.scale))
