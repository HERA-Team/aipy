#! /bin/bash
cat ../../src/ant.py | ./py2tex -n -d plain > ant.pt
cat ../../src/coord.py | ./py2tex -n -d plain > coord.pt
cat ../../src/deconv.py | ./py2tex -n -d plain > deconv.pt
cat ../../src/fit.py | ./py2tex -n -d plain > fit.pt
cat ../../src/healpix.py | ./py2tex -n -d plain > healpix.pt
cat ../../src/img.py | ./py2tex -n -d plain > img.pt
cat ../../src/loc.py | ./py2tex -n -d plain > loc.pt
cat ../../src/map.py | ./py2tex -n -d plain > map.pt
cat ../../src/miriad.py | ./py2tex -n -d plain > miriad.pt
cat ../../src/scripting.py | ./py2tex -n -d plain > scripting.pt
cat ../../src/sim.py | ./py2tex -n -d plain > sim.pt
cat ../../src/src.py | ./py2tex -n -d plain > src.pt
