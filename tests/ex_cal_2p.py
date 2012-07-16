import aipy as a

def get_aa(freqs):
    # Define the location of your instrument
    lat, lon = '45:00', '90:00'
    # Create a model of the primary beam.  BeamFlat is a minimal model 
    # that has unity gain in all directions.  
    beam = a.fit.Beam(freqs)
    # Make a list of antennas with requisite nanosecond locations, 
    # primary beams, and any other calibration parameters you wish to provide.
    ants = [
        a.pol.Antenna(  0,   0,  0, beam, num=0, pol='x'),
        a.pol.Antenna(  0, 100,  0, beam, num=0, pol='y'),
        a.pol.Antenna(100,   0,  0, beam, num=1, pol='x'),
        a.pol.Antenna(100, 100,  0, beam, num=1, pol='y'),
        a.pol.Antenna(  0,   0,  0, beam, num=2, pol='x'),
        a.pol.Antenna(  0, 100,  0, beam, num=2, pol='y'),
        a.pol.Antenna(100,   0,  0, beam, num=3, pol='x'),
        a.pol.Antenna(100, 100,  0, beam, num=3, pol='y'),
    ]
    # Create an AntennaArray at the specified location with the listed antennas
    aa = a.pol.AntennaArray((lat,lon), ants)
    return aa

def get_catalog(srcs=None, cutoff=None, catalogs=[]):
    # Pass off the request for sources to the AIPY source catalog.  If desired, 
    # you can substitute your own sources or source calibrations.
    return a.src.get_catalog(srcs=srcs, cutoff=cutoff, catalogs=catalogs)

