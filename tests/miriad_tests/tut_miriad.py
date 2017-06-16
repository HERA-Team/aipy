# --------------------------------------------------------------------

import aipy
uv = aipy.miriad.UV('test.uv')
print uv.items()
# ['vartable', 'obstype', 'history']
print uv['history']
# C2M (Python): Version=0.1.1.Fixed bandpass inversion & ordering, and pol 
# label.APPLY_BP: version=0.0.1, corr type = combXRFI: version 0.0.2XTALK2: 
# version 0.0.1 Miniaturized...
print uv.vars()
# ['latitud', 'npol', 'nspect', 'obsdec', 'vsource', 'ischan', 'operator', 
# 'nants', 'baseline', 'sfreq', 'inttime', 'source', 'epoch', 'version', 
# 'ra', 'restfreq', 'nschan', 'sdf', 'corr', 'freq', 'longitu', 'nchan', 
# 'tscale', 'antpos', 'telescop', 'pol', 'coord', 'veldop', 'lst', 'time', 
# 'dec', 'obsra']
print uv['nchan']
# 64
print uv['antpos']
# [  -8.48  205.47  187.1  -262.7   455.28  319.53 -352.95 -219.07    9.82
#  -251.71 -232.59  318.7 ]

# --------------------------------------------------------------------

preamble, data = uv.read()
print preamble
# (array([ 0.,  0.,  0.]), 2454302.8700115741, (0, 0))
print data
# [(3.55898427963+0j) (5.16037225723+0j) (7.65382957458+0j)
#  (11.5349502563+0j) (17.6214637756+0j) (26.8085384369+0j)
#  (40.0749702454+0j) (56.860118866+0j) (74.8811569214+0j) (89.6064910889+0j)
#  (98.601524353+0j) (101.491455078+0j) (100.617973328+0j) (98.0315933228+0j)
#  (95.0735092163+0j) (92.583152771+0j) (90.0556259155+0j) (88.2838745117+0j)
#  (86.3324737549+0j) (84.3934631348+0j) (82.3522338867+0j) --
#  (77.4334640503+0j) (74.7851333618+0j) (71.8084716797+0j)
#  (68.7729568481+0j) (65.6971817017+0j) (62.5315704346+0j) (59.719078064+0j)
#  (56.9530410767+0j) (54.4193191528+0j) (52.1953392029+0j)
#  (50.2718162537+0j) (48.5867958069+0j) (47.1000137329+0j)
#  (45.9260749817+0j) (44.9503746033+0j) (44.1512298584+0j) (43.609172821+0j)
#  (43.2684516907+0j) (43.1135787964+0j) (42.8874664307+0j)
#  (42.9587059021+0j) (43.0020713806+0j) (43.1228713989+0j)
#  (43.1600418091+0j) (43.1321640015+0j) (43.1135787964+0j)
#  (43.0020713806+0j) (42.726398468+0j) (42.5312576294+0j) (42.1409759521+0j)
#  (41.6794548035+0j) (41.0073051453+0j) (40.369228363+0j) (39.5948638916+0j)
#  (38.8019142151+0j) (38.0523262024+0j) (37.1168937683+0j)
#  (36.1814575195+0j) (35.2924880981+0j) (34.3105926514+0j)
#  (33.4278144836+0j) (32.3839683533+0j)]
print uv['pol'], aipy.miriad.pol2str[uv['pol']]
# -6 yy
preamble, data = uv.read()
print preamble
# (array([-538.6,  298.79613781, -674.73816035]), 2454302.8700115741, (1, 3))
preamble, data, flags = uv.read(raw=True)

# --------------------------------------------------------------------

uv.rewind()
uv.select('antennae', 0, 1, include=True)
for preamble, data in uv.all():
    uvw, t, (i,j) = preamble
    print i, j, t
# 0 1 2454302.87001
# 0 1 2454302.87009
# 0 1 2454302.87017
# 0 1 2454302.87025
# [snip]
# 0 1 2454302.91135
# 0 1 2454302.91144
# 0 1 2454302.91152
# 0 1 2454302.9116

# --------------------------------------------------------------------

import aipy
uvi = aipy.miriad.UV('test.uv')
uvo = aipy.miriad.UV('new1.uv', status='new')
uvo.init_from_uv(uvi)
def conjugate_01(uv, preamble, data):
    uvw, t, (i,j) = preamble
    if i == 0 and j == 1: return preamble, data.conjugate()
    else: return preamble, data

uvo.pipe(uvi, mfunc=conjugate_01, append2hist="Conjugated (0,1)\n")
del(uvo)

# --------------------------------------------------------------------

uvi = aipy.miriad.UV('new1.uv')
uvo = aipy.miriad.UV('new2.uv', status='new')
uvo.init_from_uv(uvi, override={'pol':-7}, exclude=['ra','lst'])
uvo.pipe(uvi)
del(uvo)

# --------------------------------------------------------------------

uvi = aipy.miriad.UV('new2.uv')
uvo = aipy.miriad.UV('new3.uv', status='new')
def change_pol(uv, p, d):
    uvw, t, (i,j) = p
    if i == j: uvo['pol'] = -5
    else: uvo['pol'] = -6
    return p, d

uvo.init_from_uv(uvi, override={'pol':-7})
uvo.pipe(uvi, mfunc=change_pol)
del(uvo)

# --------------------------------------------------------------------

import aipy, numpy
uv = aipy.miriad.UV('newest.uv', 'new')
uv['history'] = 'Made this file from scratch.\n'
uv.add_var('nchan', 'i')
uv.add_var('pol', 'i')
uv['nchan'] = 4
uv['pol'] = -5
uvw = numpy.array([1,2,3], dtype=numpy.double)
preamble = (uvw, 12345.6789, (0,1))
data = numpy.ma.array([1j,2,3j,4], mask=[0,0,1,0], dtype=numpy.complex64)
uv.write(preamble,data)
uv['pol'] = -6
uv.write(preamble,data)
del(uv)
uv = aipy.miriad.UV('newest.uv')
for p, d in uv.all():
    print p, uv['pol']
    print d

# (array([ 1.,  2.,  3.]), 12345.678900000001, (0, 1)) -5
# [1j (2+0j) -- (4+0j)]
# (array([ 1.,  2.,  3.]), 12345.678900000001, (0, 1)) -6
# [1j (2+0j) -- (4+0j)]

# --------------------------------------------------------------------

