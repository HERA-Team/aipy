# --------------------------------------------------------------------

import aipy
o = aipy.ant.ArrayLocation(location=('40:00:00','123:00:00'))
print (o.lat, o.long), (float(o.lat), float(o.long))
# (40:00:00.00, 123:00:00.00) (0.69813170079773179, 2.1467549799530254)
o.set_ephemtime('1980/8/28 12:00') # Set a time in UTC (number or string)
print aipy.ant.ephem2juldate(o.date)
# 2444480.0
o.set_jultime(2444480.0) # Same thing using Julian dates
s1 = aipy.ant.RadioFixedBody('20:00','40:00',name="CygnusA")
# ra,dec of source can be strings or radians.  Coords in Epoch 2000.
s2 = aipy.ant.RadioSpecial("Sun")
s1.compute(o)
print (s1.ra, s1.dec), (float(s1.ra), float(s1.dec))
# (19:59:22.72, 39:57:02.81) (5.2332763671875, 0.6972726583480835)
s2.compute(o)
print (s2.ra, s2.dec), (s2.az, s2.alt)
# (1:22:19.41, 8:39:55.77) (298:56:16.64, -18:15:35.20)
print s1.map
# [[ 0.84075657 -0.54141333  0.        ]
#  [ 0.08157227  0.12667295  0.98858481]
#  [-0.535233   -0.83115917  0.15066543]]

print ' --------------------------------------------------------------------'

cat = aipy.ant.SrcCatalog([s1, s2])
cat.add_src(aipy.ant.RadioSpecial("Moon"))
print cat.keys()
# ['Sun', 'CygnusA', 'Moon']
o.set_ephemtime("2007/11/20 15:15")
cat.compute(o)
print (cat['Moon'].ra, cat['Moon'].dec), cat['CygnusA'].alt
# (0:02:09.15, 2:33:27.7) 11:45:10.4

print ' --------------------------------------------------------------------'

import aipy, numpy
freqs = numpy.array([.150,.160,.170])
beam = aipy.ant.Beam(freqs)
ants = []
ants.append(aipy.ant.Antenna(0,0,0,beam,delay=0,offset=0))
ants.append(aipy.ant.Antenna(0,100,0,beam,delay=1))
ants.append(aipy.ant.Antenna(100,0,0,beam,offset=.5))
aa = aipy.ant.AntennaArray(ants=ants,location=("18:20:39","-66:45:10"))
print aa.get_baseline(0,2,'r'), aa.get_delay(1,2), aa.get_offset(0,2)
# [ 100.    0.    0.] -1.0 0.5
aa.set_jultime(2454447.37472)
srcs = []
srcs.append(aipy.ant.RadioSpecial("Sun"))
srcs.append(aipy.ant.RadioSpecial("Venus"))
cat = aipy.ant.SrcCatalog(srcs)
cat.compute(aa) # REMEMBER to call this before projecting!
print aa.get_baseline(0,1,src=cat['Sun'])
# [ 34.66645092 -36.79755511 -86.27964487]
#print aa.get_baseline(0,1,src=cat['Venus'])
try: print aa.get_baseline(0,1,src=cat['Venus'])
except aipy.ant.PointingError,e: print aipy.ant.PointingError, e
# <class 'aipy.ant.PointingError'>: 'Venus is below horizon'

print ' --------------------------------------------------------------------'

print aa.gen_phs(cat['Sun'], 0, 1)
# [ 0.26051837-0.96546889j -0.61418129-0.78916496j -0.99988051+0.01545836j]
data = aa.phs2src(numpy.array([.5,1,0j]),cat['Sun'],1,2)
print data
# [ 0.44540993-0.22717833j  0.02050412-0.99978977j  0.        -0.j        ]
uvw = aa.gen_uvw(1,2,cat['Sun'])
print uvw
# [[  8.86987019   7.55959383  17.72506541]
#  [  9.46119487   8.06356676  18.90673644]
#  [ 10.05251955   8.56753968  20.08840747]]
print aa.unphs2src(data,cat['Sun'],1,2)
# [ 0.5+0.j  1.0+0.j  0.0+0.j]
print aa.rmsrc(numpy.array([.5,1,0j]),cat.values(),1,2)
# [ 0.36514582+0.25266785j  0.39471115-0.13255143j  0.12804513-0.41987284j]

print ' --------------------------------------------------------------------'

import aipy
uvi = aipy.miriad.UV('test.uv')
aa = aipy.loc.get_aa('pwa303',uvi['sdf'],uvi['sfreq'],uvi['nchan'])
print len(aa.ants), aa.get_baseline(0,1)
# 4 [ 138.14       -222.37112004   -0.70162169]
cat = aipy.src.get_catalog(srcs=['Sun','cen'])
uvo = aipy.miriad.UV('test.uv.mod','new')
uvo.init_from_uv(uvi)
def f(uv,preamble,data):
    uvw,t,(i,j) = preamble
    aa.set_jultime(t)
    cat.compute(aa)
    data = aa.rmsrc(data, cat['Sun'], i, j)
    data = aa.phs2src(data, cat['cen'], i, j)
    return preamble, data

uvo.pipe(uvi, mfunc=f)
del(uvo)

