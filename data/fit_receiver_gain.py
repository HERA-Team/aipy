#! /usr/bin/env python
import numpy, pylab

POLY_ORDER = 8
F_START = 30
F_END = -20

f = open('Receiver_Gain.txt')
lines = [map(float, L.split()) for L in f.readlines() if L[0].isdigit()]
data = numpy.array(lines)
data[:,0] /= 1e9        # Put in GHz
x0 = data[:,0]
xf = x0[F_START:F_END]
y0 = 10**(data[:,1]/10) / 1e12
yf = y0[F_START:F_END]

pylab.plot(x0, y0, label='data')
for i in range(5, POLY_ORDER):
    p = numpy.polyfit(xf, yf, i+1)
    y1 = numpy.polyval(p, x0)
    pylab.plot(x0, y1.clip(0, y1.max()), label=str(i+1))
    if i+1 == 6: print p

pylab.legend()
pylab.show()

