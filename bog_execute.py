#!/bin/python
import sys
import h5py
import numpy as np
#from mayavi import mlab
import bog_fct as bfct
import sigclip

# run bog_execute.py <model> <run> <snap> 
# run bog_execute.py G00 0 100

# Units are GADGET units: kpc/h, km/s, M_sun/h

workpath="/Users/coralwheeler/Work/Dwarf/MetalMixing"

model=sys.argv[1]
run=sys.argv[2]
snap = sys.argv[3]
snapfull = 'Snapshot00%s'%(snap)
fname = "%s/%s_%s/irate/%s_%s_%s-irate.hdf5"%(workpath, model, run, model, run, snap)

# Check to see if file has stars
pf = h5py.File(fname,'r')
if len(pf[snapfull]['ParticleData'].keys()) < 3: # no stars case
    arestars = bool(False)
else:
    arestars = bool(True)

# Create dark matter particle snapshot object
dm = bfct.DarkParticles(fname, snapfull)
gas = bfct.GasParticles(fname, snapfull)
if arestars:
    stars = bfct.StarParticles(fname, snapfull)

# Recenter based on location of ALL particles if star particle exists
if arestars:
    pos = np.concatenate((dm.pos, stars.pos, gas.pos), axis=0)
    m = np.concatenate((dm.arrays.mass, stars.arrays.mass, gas.arrays.mass))
    v = np.concatenate((dm.vel, stars.vel, gas.vel), axis=0)
else: # Only use dark matter and gas if no stars
    pos = np.concatenate((dm.pos, gas.pos), axis=0)
    m = np.concatenate((dm.arrays.mass, gas.arrays.mass))
    v = np.concatenate((dm.vel, gas.vel), axis=0)

dm.recenter(pos,m,v)
gas.recenter(pos,m,v)
if arestars:
    stars.recenter(pos,m,v)

gasorder = bfct.get_order(gas.arrays.r)
rgasorder = gas.arrays.r[gasorder]
mgasorder = gas.arrays.mass[gasorder]