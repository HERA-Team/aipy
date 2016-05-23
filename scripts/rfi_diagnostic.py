#!/usr/bin/python
"""
Diagnose RFI and RFI removal strategies by plotting either maximum amplitude per channel or fraction of flagged data per channel.
"""
import numpy as np
import aipy as a
from matplotlib import pylab as p
import os,sys,optparse,math

#Set up the options, make sure to define them well.

o = optparse.OptionParser()
o.set_usage('rfi_diagnostic.py [options] *.uvcb')
o.set_description(__doc__)

a.scripting.add_standard_options(o,ant=True,pol=True,chan=False,dec=False,cmap=False,max=False,drng=False)
o.add_option('-m','--max_hold',action='store_true',dest='max_hold',help='Display the maximum amplitude in each frequency bin.')
o.add_option('-d','--dwell_time',action='store_true',dest='dwell',help='Display the fraction of time lost in each frequency bin.')
o.add_option('-x','--x_axis',dest='x_axis',default='freq',help='Indicates the units of the x-axis. freq [default] = frequency (MHz), z = redshift,chan = channel number')
o.add_option('-l','--log',dest='log',action='store_true',help='Display max-hold values in a logarithmic scale. This is recommended.')
o.add_option('-s','--str',dest='strategy',default='normal',help='Chooses the RFI removal strategy. simple = the strategy of xrfi_simple.py. normal [default] = the strategy of xrfi.py. naive = another strategy of xrfi_simple.py, in which those data a certain number of standard deviations above the median are flagged. none = RFI removal has already been performed. Only count flagged data. This is moved to default if filename ends in r or R.')
o.add_option('--best',dest='best',action='store_true',help='Summarizes all baselines. Dispays the median, lowest, and highest average dwell-times (or max-holds) for all selected baselines.')
o.add_option('-g','--gom',dest='gom',default=None,help='Flag the gainometer channels, because we dont care about it.')
o.add_option('--sum',dest='sum',action='store_true',help='Average over all baselines, plot one line.')

opts, args = o.parse_args(sys.argv[1:])

#Run through the uv files and plot what you want to...

plot_y = {}
infected = {}
total = {}
hi_bl,med_bl,lo_bl = '','',''
freqs = []

for uvfile in args:
	print 'Reading',uvfile
	uv = a.miriad.UV(uvfile)
	a.scripting.uv_selector(uv,opts.ant,opts.pol)	
	if uvfile == args[0]:
		freqs = 1000.*np.linspace(uv['sfreq'],uv['sfreq']+uv['nchan']*uv['sdf'],uv['nchan']) #Units are sensible Megahertz.
	if uvfile[-1] == 'r' or uvfile[-1] == 'R': opts.strategy = 'none'
	data,mask = {},{}
	#This is the max-hold part of the process.
	for (uvw,t,(i,j)),d,f in uv.all(raw=True):
		bl = '%d,%d' % (i,j)
		if i==opts.gom or j==opts.gom: continue #Don't look at gainometer channels.
		f = np.where(np.abs(d) == 0,1,f)
		if opts.max_hold:
			if opts.log: d = np.log10(d)
			if not plot_y.has_key(bl): plot_y[bl] = np.zeros(uv['nchan'])
			for I in range(uv['nchan']):
				if d[I] >= plot_y[bl][I]: plot_y[bl][I]=d[I]
		#This is the part that figures out fractional time spent in RFI-mode. Use xrfi filter (derivatives, not polynomial fit) and just compute (flagged_samples/total_samples) for each channel.
		if opts.dwell:
			#Following xrfi.py, record all datums... make a mask later.
			if not data.has_key(bl):
				data[bl]={}
				mask[bl] = {}
			#create data array, mask.
			data[bl][t] = d
			mask[bl][t] = f	
		if not opts.max_hold and not opts.dwell:
			print 'Please indicate whether you want max-hold (-m) or dwell-time (-d) plots'
			sys.exit(0)
	
	if opts.dwell:
		for bl in data: #This is a blatant rip-off of xrfi_simple.py, but that's kind of what we want...	
			data_times = data[bl].keys()
			data_times.sort()	
			rfi_counter = np.array([mask[bl][t] for t in data_times]) 
			
			if opts.strategy == 'simple':
				d = np.array([data[bl][t] for t in data_times])
				#find sigma for frequency.
				ddf = d[:,1:-1]-0.5*(d[:,:-2]+d[:,2:])
				ddf2 = np.abs(ddf)**2
				sigf = np.sqrt(np.median(ddf2,axis=1))
				sigf.shape = (sigf.size,1)
				rfi_counter[:,0] += 1;rfi_counter[:,-1] += 1
				rfi_counter[:,1:-1] += np.where(ddf2/sigf**2 > 4.**2,1,0)
				#same thing, but for time.
				ddt = d[1:-1,:]-0.5*(d[:-2,:]+d[2:,:])
				ddt2 = np.abs(ddt)**2
				sigt = np.sqrt(np.median(ddt2,axis=1))
				sigt.shape = (sigt.size,1)
				rfi_counter[0,:] += 1;rfi_counter[-1,:] += 1
				rfi_counter[1:-1,:] += np.where(ddt2/sigt**2 > 4.**2,1,0)			
	
				rfi_counter = np.where(rfi_counter > 0,1,0)
	
			if opts.strategy == 'naive':
				ad = np.abs(d)
				med = np.median(ad)
				sig = np.sqrt(np.median(np.abs(ad-med)**2))
				rfi_counter += np.where(ad > med +2.*sig,1,0)
				
				rfi_counter = np.where(rfi_counter > 0,1,0) #Don't want to count twice!!!
			
			if opts.strategy == 'normal':
				d = np.ma.array([data[bl][t] for t in data_times],mask=[mask[bl][t] for t in data_times])
				hi,low = a.rfi.gen_rfi_thresh(d)
				m = np.where(np.abs(d) > hi,1,0)
				rfi_counter += m
				ch_cnt = np.array([mask[bl][t] for t in data_times]).sum(axis=0)
				ch_msk = np.where(ch_cnt > ch_cnt.max()*0.33,1,0)
				for I,t in enumerate(rfi_counter):
					if rfi_counter[I].sum() > rfi_counter[I].size*0.33: rfi_counter[I] += np.ones(uv['nchan'])
					else: rfi_counter[I] += ch_msk		
				if i == j:
					bad_ints = a.rfi.flag_by_int(d)
					for I in np.where(bad_ints)[0]:
						rfi_counter[I] += np.ones(uv['nchan'])			
				
				rfi_counter = np.where(rfi_counter > 0,1,0)
			
			if opts.strategy == 'none':
				pass

			#Now, to count everything.
			if not infected.has_key(bl):
				infected[bl] = np.zeros(uv['nchan'])
				total[bl] = np.zeros(uv['nchan'])
			for ct in range(len(rfi_counter)):
				infected[bl] += rfi_counter[ct]
				total[bl] += np.ones(uv['nchan'])	

	del(uv)	

if opts.dwell:
	for bl in infected:
		plot_y[bl] = infected[bl] / total[bl]
		if opts.log: plot_y[bl] = np.log10(plot_y[bl])

#Housekeeping...

#Do some conversion of the frequency axis.
if opts.x_axis == 'freq': 
	plot_x = freqs
	xlabel = 'Frequency (MHz)'
if opts.x_axis == 'z': 
	plot_x = (1427.1/freqs)-1.
	xlabel = 'Redshift of EoR'
if opts.x_axis == 'chan':
	plot_x = np.arange(nchan)
	xlabel = 'Channel Number'

bls = plot_y.keys() #This keeps track of all the baselines.
if len(bls) == 0:
        print 'No data to plot'
        sys.exit(0)

def sort_func(a,b): #This guy will make sure that your baselines always show up in order.
	ai,aj = map(int,a.split(','))
	bi,bj = map(int,b.split(','))
	if bi > ai or (bi == ai and bj > aj): return -1
	return 1
bls.sort(cmp=sort_func)

m2 = int(math.sqrt(len(bls))) #These are to make the plots pretty -- to figure where to place subplots.
m1 = int(math.ceil(float(len(bls))) / m2)

if opts.max_hold and not opts.log: ylabel = 'Max Amplitude (counts)'
if opts.max_hold and opts.log: ylabel = 'Max Amplitude (log10(counts))'
if opts.dwell and not opts.log: ylabel = 'Fractional Dwell Time'
if opts.dwell and opts.log: ylabel = 'Fractional Dwell Time (log)'

#scan through everything if opts.best is chosen...
if opts.best:
	#recover the bl that corresponds the the median, hi, and lo of average dwell_time or max_hold.
	#First, i have to take the averages.
	averages = {}
	hi_bl,med_bl,lo_bl = '','',''
	for bl in plot_y: averages[bl] = np.average(plot_y[bl])
	#I have to do this, because if the number of baselines is even, the median is not contained within averages[]
	values = averages.values()
	values.sort()
	nmed = int(len(averages)/2)
	hi,med,lo = np.max(values),values[nmed],np.min(values)
	for count,bl in enumerate(averages):
		if averages[bl] == hi: hi_bl = bl
		if averages[bl] == med: med_bl = bl
		if averages[bl] == lo: lo_bl = bl
	#Now, plot the three on the same plot...
	p.plot(plot_x,plot_y[lo_bl],'g-',label='Best='+lo_bl)
	p.plot(plot_x,plot_y[med_bl],'b-',label='Median='+med_bl)
	p.plot(plot_x,plot_y[hi_bl],'r-',label='Worst='+hi_bl)
	p.ylabel(ylabel),p.xlabel(xlabel),p.legend()
	p.show()
	sys.exit(0)

if opts.sum:
	av = np.zeros(len(plot_x))
	num = float(len(plot_y))
	for bl in plot_y:
		av += plot_y[bl]
	av /= num
	p.plot(plot_x,av)
	p.xlabel(xlabel);p.ylabel('Average'+ylabel)
	p.show()
	sys.exit(0)

#Now to make the plots.
for counts,bl in enumerate(bls):
	p.subplot(m2,m1,counts+1)
	p.xlim((plot_x[0],plot_x[-1]))
	p.plot(plot_x,plot_y[bl])
	p.ylabel(ylabel);p.xlabel(xlabel);p.title(bl)
p.show()
