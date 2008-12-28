"""
An astronomy library for precessing coordinates between epochs and converting 
between 'Ecliptic', 'Equatorial', and 'Galactic' coordinate systems.
Vectors are 3 dimensional with magnitude 1, representing a point on the
unit sphere.  Also includes generic 3 vector rotation code.

Precession and coordinate transformation code was adapted to Python from code 
written as a part of Healpix, copyright (C) 2005 Max-Planck-Society, 
Martin Reinecke, and released under the terms of the GNU General Public 
License.

Author: Aaron Parsons
Date: 11/14/2007
Revisions: None
"""

import numpy as n

degr2rad = n.pi / 180.

def precess(v, epoch, iepoch=2000.):
    """Precess a vector from 'iepoch' to the 'epoch' specified."""
    # Z-axis rotation by OBL_LONG
    Tm = ((epoch+iepoch)*0.5 - 1900.) *0.01
    gp_long  = (epoch-iepoch) * (50.2564+0.0222*Tm) / 3600.
    obl_long = 180. - (173. + (57.06+54.77*Tm) / 60.) + gp_long*0.5
    dco, dso = n.cos(obl_long*degr2rad), n.sin(obl_long*degr2rad)
    v[...,0], v[...,1] = v[...,0]*dco-v[...,1]*dso, v[...,0]*dso+v[...,1]*dco
    # X-axis rotation by dE:
    dE = (epoch-iepoch) * (0.4711-0.0007*Tm) / 3600.
    dce, dse = n.cos(dE*degr2rad), n.sin(dE*degr2rad)
    v[...,1], v[...,2] = v[...,1]*dce-v[...,2]*dse, v[...,1]*dse+v[...,2]*dce
    # Z-axis rotation by GP_LONG - OBL_LONG:
    dL = gp_long - obl_long
    dcl, dsl = n.cos(dL*degr2rad), n.sin(dL*degr2rad)
    v[...,0], v[...,1] = v[...,0]*dcl-v[...,1]*dsl, v[...,0]*dsl+v[...,1]*dcl
    return v

def get_epsilon(epoch):
    T = (epoch - 1900.) * 0.01
    epsilon = 23.452294 - 0.0130125*T - 1.63889e-6*T**2 + 5.02778e-7*T**3
    return epsilon*degr2rad

def e_to_q(v, epoch=2000.):
    """Routine to convert from ecliptic (celestial) coordinates to equatorial
    coordinates at the given epoch."""
    epsilon = get_epsilon(epoch)
    dc, ds = n.cos(epsilon), n.sin(epsilon)
    v[...,1], v[...,2] = dc*v[...,1]-ds*v[...,2], dc*v[...,2]+ds*v[...,1]
    return v

def q_to_e(v, epoch=2000.):
    """Convert equatorial coordinates to ecliptic (celestial) coordinates
    at the given epoch."""
    epsilon = -get_epsilon(epoch)
    dc, ds = n.cos(epsilon), n.sin(epsilon)
    v[...,1], v[...,2] = dc*v[...,1]-ds*v[...,2], dc*v[...,2]+ds*v[...,1]
    return v

def g_to_e(v, epoch=2000.):
    """Convert galactic coordinates to ecliptic (celestial) 
    coordinates at the given epoch.  First the conversion to ecliptic 
    2000 is done, then the results are precessed (if necessary)."""
    T = n.array([[-0.054882486,  0.494116468, -0.867661702,],
                 [-0.993821033, -0.110993846, -0.000346354,],
                 [-0.096476249,  0.862281440,  0.497154957,]])
    v = n.dot(T, v)
    if epoch != 2000.: v = precess(v, 2000., iepoch=epoch)
    return v

def e_to_g(v, epoch=2000.):
    """Convert ecliptic (celestial) coordinates at the given
    epoch to galactic coordinates.  The ecliptic coordinates are first
    precessed to 2000. (if necessary), then converted."""
    T = n.array([[-0.054882486, -0.993821033, -0.096476249,],
                 [ 0.494116468, -0.110993846,  0.862281440,],
                 [-0.867661702, -0.000346354,  0.497154957,]])
    if epoch != 2000.: v = precess(v, 2000., iepoch=epoch)
    return n.dot(T, v)

def convert(v, isys, osys, iepoch=2000., oepoch=2000.):
    """Convert vector 'v' from coordinate system 'isys' to 'osys', including
    precession.  Valid coordinate systems are 'Ecliptic', 'Equatorial', and
    'Galactic'."""
    if isys == 'Ecliptic': pass
    elif isys == 'Equatorial': v = q_to_e(v, iepoch)
    elif isys == 'Galactic': v = g_to_e(v, iepoch)
    else: raise ValueError("Unknown input coordinate system")
    if iepoch != oepoch: v = precess(v, oepoch, iepoch=iepoch)
    if osys == 'Ecliptic': pass
    elif osys == 'Equatorial': v = e_to_q(v, oepoch)
    elif osys == 'Galactic': v = e_to_g(v, oepoch)
    else: raise ValueError("Unknown output coordinate system")
    return v

def coordsys2matrix(isys, osys, iepoch=2000., oepoch=2000.):
    """Return the rotation matrix corresponding to a coordinate 
    system/precssion transformation (see 'convert')."""
    m = n.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    return convert(m, isys, osys, iepoch=iepoch, oepoch=oepoch)

def rotmatrix(ang, vec):
    """Return 3x3 rotation matrix defined by rotation by 'ang' around the
    axis 'vec' (according to the right-hand rule).  Both can be vectors,
    returning a vector of rotation matrices.  Rotation matrix will have a 
    scaling of |vec| (i.e. normalize |vec|=1 for a pure rotation)."""
    c = n.cos(ang); s = n.sin(ang); C = 1-c
    x,y,z = vec[...,0], vec[...,1], vec[...,2]
    xs,ys,zs = x*s, y*s, z*s
    xC,yC,zC = x*C, y*C, z*C
    xyC,yzC,zxC = x*yC, y*zC, z*xC
    rm = n.array([[x*xC+c, xyC-zs, zxC+ys],
                  [xyC+zs, y*yC+c, yzC-xs],
                  [zxC-ys, yzC+xs, z*zC+c]])
    axes = range(rm.ndim)
    return rm.transpose(axes[-1:] + axes[:-1])

