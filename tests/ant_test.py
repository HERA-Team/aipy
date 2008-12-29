#! /usr/bin/env python
import aipy as a, numpy as n, ephem as e

JULDATE = 2454000
epoch = a.ant.juldate2ephem(JULDATE)
#epoch = e.J2000

src1 = a.ant.RadioFixedBody('12:00','45:00', epoch=epoch, name='src1')
src2 = a.ant.RadioFixedBody('12:00','90:00', epoch=epoch, name='src2')
sun = a.ant.RadioSpecial('Sun')
cat = a.ant.SrcCatalog([src1,src2])

bm = a.ant.Beam(n.array([.1,.2,.3]), active_chans=n.array([0,1]))
ant1 = a.ant.Antenna(0, 0, 0, bm)
ant2 = a.ant.Antenna(100, 0, 0, bm)
aa = a.ant.AntennaArray(('0:00','0:00'), [ant1,ant2])

aa.set_jultime(JULDATE)
sun.compute(aa)
cat.compute(aa)

print '---------------------------------------------------------------'
print '(0,1) eq raw:', aa.get_baseline(0,1,src='r')
print '(0,1) eq now:', aa.get_baseline(0,1,src='e')
print '(0,1) top:', aa.get_baseline(0,1,src='z')
print '---------------------------------------------------------------'
print 'Sun at RA=%s, DEC=%s' % (sun.ra, sun.dec)
print '       AZ=%s, ALT=%s' % (sun.az, sun.alt)
print '(0,1) Sun:', aa.get_baseline(0,1,src=sun)
print '---------------------------------------------------------------'
print 'src1 at RA=%s, DEC=%s' % (cat['src1'].ra, cat['src1'].dec)
print '        AZ=%s, ALT=%s' % (cat['src1'].az, cat['src1'].alt)
print '(0,1) src1:', aa.get_baseline(0,1,src=cat['src1'])
print '---------------------------------------------------------------'
print 'src2 at RA=%s, DEC=%s' % (cat['src2'].ra, cat['src2'].dec)
print '        AZ=%s, ALT=%s' % (cat['src2'].az, cat['src2'].alt)
print '(0,1) src2:'
x,y,z = aa.get_baseline(0,1,src=cat.get_crds('eq', srcs=['src2']))
print 'x=%s\ny=%s\nz=%s' % (x,y,z)
print '---------------------------------------------------------------'
print '(0,1) %s:' % cat.keys()
x,y,z = aa.get_baseline(0,1,src=cat.get_crds('eq'))
print 'x=%s\ny=%s\nz=%s' % (x,y,z)
print '---------------------------------------------------------------'
print 'uvw Sun for freqs %s:' % (bm.afreqs)
u,v,w = aa.gen_uvw(0,1, src=sun)
print 'u=%s\nv=%s\nw=%s' % (u,v,w)
print '---------------------------------------------------------------'
print 'uvw %s:'% cat.keys()
u,v,w = aa.gen_uvw(0,1, src=cat.get_crds('eq'))
print 'u=%s\nv=%s\nw=%s' % (u,v,w)
print '---------------------------------------------------------------'
print aa.resolve_src(cat.get_crds('eq', srcs=['src2']), 0,1, 
    srcshape=(.02,.01,0.))
print aa.resolve_src(cat['src2'], 0, 1, srcshape=(.02,.01,n.pi/2))
print aa.resolve_src(cat['src2'], 0, 1, srcshape=(.01,.02,0.))
print aa.resolve_src(cat.get_crds('eq'), 0, 1, srcshape=(.01,.02,0.))

