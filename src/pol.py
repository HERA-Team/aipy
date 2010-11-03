import fit,numpy as n
def invdict(d):
    return dict([[v,k] for k,v in d.iteritems()])
def main():
    import aipy as a
    aa = a.cal.get_aa('pgb050_v002',0.1,0.1,1)
    prms = aa.get_params(ant_prms={'*':['x','y','z','dly','phs','phsoff']})
    #print prms
    fit.print_params(aa.get_params(ant_prms={'*':['x','y','z','dly','amp']}))

    
if __name__=='__main__':
    main()
class Antenna(fit.Antenna):
    """Representation of physical location and beam pattern of individual 
    antenna in array."""
    def __init__(self, x, y, z, beam, phsoff=[0.,0.], bp_r=n.array([1]),
            bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0),pol='x',
            num=-1,**kwargs):
        """x,y z = antenna coordinates in equatorial (ns) coordinates
        beam = Beam object (implements response() function)
        phsoff = polynomial phase vs. frequency.  Phs term that is linear
                 with freq is often called 'delay'.
        bp_r = polynomial (in freq) modeling real component of passband
        bp_i = polynomial (in freq) modeling imaginary component of passband
        amp = overall multiplicative scaling of gain
        pointing = antenna pointing (az,alt).  Default is zenith.
        pol = 'x'(default) or 'y'
        """
        fit.Antenna.__init__(self, x, y, z, beam, phsoff=phsoff,
            bp_r=bp_r,bp_i=bp_i,
            amp=amp, pointing=pointing, **kwargs)
        self.pol = pol
        self.num = num
class AntennaArray(fit.AntennaArray):        
    def get_ant_list(self):
        try: #return [str(ant.num)+ant.pol for ant in self.ants]
            ants = {}
            for i,ant in enumerate(self.ants):
                ants[str(ant.num)+ant.pol]=i
            return ants
        except(NameError): return [str(i) for i in range(len(self.ants))]
    def get_params(self, ant_prms={'*':'*'}):
        """Return all fitable parameters in a dictionary."""
        prms = {}
        for k in ant_prms:
            ants = self.get_ant_list()            
            if k.startswith('*'):
                ants = self.get_ant_list()
   #             ants = map(str, range(len(self)))
            else: ants = {k:ants[k]}
           # print ants
            prm_list = ant_prms[k]
            if type(prm_list) is str: prm_list = [prm_list]
            for a,i in ants.iteritems():
     #           try: prms[a] = self.ants[a].get_params(prm_list)
                try:prms[a]=self.ants[i].get_params(prm_list)
                except(ValueError): pass
        return prms
