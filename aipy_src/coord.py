"""
An astronomy library for precessing coordinates between epochs and converting 
between topocentric (z = up, x = east), ecliptic (heliocentric), equatorial
(celestial), and galactic coordinate systems.  Vectors are 3 dimensional 
with magnitude 1, representing a point on the unit sphere.  Includes generic 
3 vector rotation code.
"""

import numpy as np
import ephem as e

sys_dict = {
    'eq': e.Equatorial,
    'ec': e.Ecliptic,
    'ga': e.Galactic
}

def convert(crd, isys, osys, iepoch=e.J2000, oepoch=e.J2000):
    """Convert 'crd' from coordinate system isys to osys, including
    epoch precession.  Valid coordinate systems are 'ec' (Ecliptic), 'eq' 
    (Equatorial), and 'ga' (Galactic).  Epochs may be date strings, or
    numerical ephemeris times."""
    if len(crd) == 3: crd = eq2radec(crd)
    c1 = sys_dict[isys[:2].lower()](crd[0], crd[1], epoch=iepoch)
    return sys_dict[osys[:2].lower()](c1, epoch=oepoch).get()

def convert_m(isys, osys, iepoch=e.J2000, oepoch=e.J2000):
    """Return the 3x3 matrix corresponding to a coordinate 
    system/precssion transformation (see 'convert'). 
    NOTE: To obtain a transformation matrix to dot with isys to get osys, 
    must reverse isys and osys arguments in convert_m. """
    m = np.array([[1,0,0],[0,1,0],[0,0,1]], dtype=np.double)
    for i in range(3):
        c = convert(m[:,i], isys, osys, iepoch=iepoch, oepoch=oepoch)
        m[:,i] = radec2eq(c)
    return m

def rot_m(ang, vec):
    """Return 3x3 matrix defined by rotation by 'ang' around the
    axis 'vec', according to the right-hand rule.  Both can be vectors,
    returning a vector of rotation matrices.  Rotation matrix will have a 
    scaling of |vec| (i.e. normalize |vec|=1 for a pure rotation)."""
    c = np.cos(ang); s = np.sin(ang); C = 1-c
    x,y,z = vec[...,0], vec[...,1], vec[...,2]
    xs,ys,zs = x*s, y*s, z*s
    xC,yC,zC = x*C, y*C, z*C
    xyC,yzC,zxC = x*yC, y*zC, z*xC
    rm = np.array([[x*xC+c, xyC-zs, zxC+ys],
                   [xyC+zs, y*yC+c, yzC-xs],
                   [zxC-ys, yzC+xs, z*zC+c]], dtype=np.double)
    if rm.ndim > 2:
        axes = range(rm.ndim)
        return rm.transpose(axes[-1:] + axes[:-1])
    else:
        return rm

def xyz2thphi(xyz):
    """Convert xyz vectors (x,y,z along first axis) into angles theta
    (from z axis), phi (counter-clockwise around z, 0 at x axis)."""
    x,y,z = xyz
    phi = np.arctan2(y, x)
    th = np.arctan2(np.sqrt(x**2+y**2),z)
    if np.ma.isMA(x):
        try:
            return np.ma.array([th.filled(0),phi.filled(0)], 
                mask=[x.mask,x.mask], dtype=np.double)
        except(np.core.ma.MAError):
            return np.ma.array([th,phi], dtype=np.double)
    return np.array([th,phi], dtype=np.double)

def thphi2xyz(th_phi):
    """Convert angles theta (from z axis), phi (counter-clockwise around z, 
    0 at x axis) into xyz vectors (x,y,z along first axis)."""
    th,phi = th_phi
    z = np.cos(th)
    r = np.sin(th)
    x,y = r*np.cos(phi), r*np.sin(phi)
    if np.ma.isMA(th):
        try:
            return np.ma.array([x.filled(),y.filled(),z.filled()], 
                mask=[th.mask, th.mask, th.mask], dtype=np.double)
        except(np.core.ma.MAError):
            return np.ma.array([x,y,z], dtype=np.double)
    return np.array([x,y,z], dtype=np.double)

def eq2radec(xyz):
    """Convert equatorial xyz vectors (x,y,z along first axis) into angles ra 
    (counter-clockwise around z = north, 0 at x axis), dec (from equator)."""
    th_phi = xyz2thphi(xyz)
    th,phi = th_phi
    dec = np.pi/2 - th
    ra = np.ma.where(phi < 0, phi + 2*np.pi, phi)
    th_phi[0],th_phi[1] = ra,dec
    return th_phi

def radec2eq(ra_dec):
    """Convert angles ra (counter-clockwise around z = north, 0 at x axis), dec 
    (from equator) into equatorial xyz vectors (x,y,z along first axis)."""
    phi,th = ra_dec
    return thphi2xyz((np.pi/2 - th, phi))

def latlong2xyz(lat_long):
    """Convert angles lat (from equator), long (counter-clockwise around
    z = north, 0 at x axis) into xyz vectors (x,y,z along first axis)."""
    lat,long = lat_long
    return radec2eq((long,lat))

def top2azalt(xyz):
    """Convert topocentric xyz vectors (x,y,z along first axis) into angles az 
    (clockwise around z = up, 0 at x axis = north), alt (from horizon)."""
    th_phi = xyz2thphi(xyz)
    th,phi = th_phi
    alt = np.pi/2 - th
    az = np.pi/2 - phi
    az = np.ma.where(az < 0, az + 2*np.pi, az)
    th_phi[0],th_phi[1] = az,alt
    return th_phi

def azalt2top(az_alt):
    """Convert angles az (clockwise around z = up, 0 at x axis = north), alt 
    (from horizon) into topocentric xyz vectors (x,y,z along first axis)."""
    az,alt = az_alt
    return thphi2xyz((np.pi/2 - alt, np.pi/2 - az))

def eq2top_m(ha, dec):
    """Return the 3x3 matrix converting equatorial coordinates to topocentric
    at the given hour angle (ha) and declination (dec)."""
    sin_H, cos_H = np.sin(ha), np.cos(ha)
    sin_d, cos_d = np.sin(dec), np.cos(dec)
    zero = np.zeros_like(ha)
    map =  np.array([[    sin_H    ,       cos_H  ,       zero  ],
                     [ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
                     [  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])
    if len(map.shape) == 3: map = map.transpose([2, 0, 1])
    return map

def top2eq_m(ha, dec):
    """Return the 3x3 matrix converting topocentric coordinates to equatorial
    at the given hour angle (ha) and declination (dec)."""
    m = eq2top_m(ha, dec)
    if len(m.shape) == 3:
        for i in range(m.shape[0]): m[i] = np.linalg.inv(m[i])
        return m
    else: return np.linalg.inv(eq2top_m(ha, dec))

