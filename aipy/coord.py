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

import numpy as n, aipy

degr2rad = n.pi / 180.

def precess(v, epoch, iepoch=2000.):
    """Precess a vector from 'iepoch' to the 'epoch' specified.  This does
    not give the same answer as ephemeris (!)."""
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
                 [-0.096476249,  0.862281440,  0.497154957,]], dtype=n.double)
    v = n.dot(T, v)
    if epoch != 2000.: v = precess(v, 2000., iepoch=epoch)
    return v

def e_to_g(v, epoch=2000.):
    """Convert ecliptic (celestial) coordinates at the given
    epoch to galactic coordinates.  The ecliptic coordinates are first
    precessed to 2000. (if necessary), then converted."""
    T = n.array([[-0.054882486, -0.993821033, -0.096476249,],
                 [ 0.494116468, -0.110993846,  0.862281440,],
                 [-0.867661702, -0.000346354,  0.497154957,]], dtype=n.double)
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
    m = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
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
                  [zxC-ys, yzC+xs, z*zC+c]], dtype=n.double)
    axes = range(rm.ndim)
    return rm.transpose(axes[-1:] + axes[:-1])

def xyz_to_th_phi(xyz):
    """Convert xyz vectors (with x,y,z along the first axis if 
    multidimensional) into theta (angle from z axis), phi (angle clockwise
    around z, with phi=0 at x axis)."""
    x,y,z = xyz
    phi = n.arctan2(y, x)
    th = n.arctan2(n.sqrt(x**2+y**2),z)
    return n.array([th,phi], dtype=n.double)

def th_phi_to_xyz(th_phi):
    th,phi = th_phi
    z = n.cos(th)
    r = n.sin(th)
    x,y = r*n.cos(phi), r*n.sin(phi)
    return n.array([x,y,z], dtype=n.double)

def xyz_to_ra_dec(xyz):
    th_phi = xyz_to_th_phi(xyz)
    th,phi = th_phi
    dec = n.pi/2 - th
    ra = n.where(phi < 0, phi + 2*n.pi, phi)
    th_phi[0],th_phi[1] = ra,dec
    return th_phi

def ra_dec_to_xyz(ra_dec):
    phi,th = ra_dec
    return th_phi_to_xyz((n.pi/2 - th, phi))

def az_alt_to_xyz(az_alt):
    az,alt = az_alt
    return th_phi_to_xyz((n.pi/2 - alt, -az))

def xyz_to_az_alt(xyz):
    th_phi = xyz_to_th_phi(xyz)
    th,phi = th_phi
    alt = n.pi/2 - th
    phi = -phi
    az = n.where(phi < 0, phi + 2*n.pi, phi)
    th_phi[0],th_phi[1] = az,alt
    return th_phi

def eq_to_topo_matrix(observer):
    aa = aipy.loc.get_loc('pwa303')
    m = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
    ra, dec = xyz_to_ra_dec(m)
    for i in range(ra.shape[0]):
        s = aipy.ant.RadioFixedBody(ra[i], dec[i])
        s.compute(observer)
        ra[i], dec[i] = s.az, s.alt
    return az_alt_to_xyz((ra,dec))

def topo_to_eq_matrix(observer):
    m = eq_to_topo_matrix(observer)
    return n.linalg.inv(m)

