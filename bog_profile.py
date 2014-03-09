from __future__ import division
import h5py
import numpy as np
import math
from matplotlib import pyplot as plt
#from mayavi import mlab
import sigclip as clip
from scipy import stats
import scipy

class Particles:
    #initialize
    def __init__(self, fname, snap):
        self.fname = fname
        self.snap = snap
        self.file = h5py.File(self.fname,'r')
        self.n = 0 # will be set in child class
        self.arrays = np.ndarray(self.n) # will be set in child class
        self.hh = self.file['Cosmology'].attrs['HubbleParam'] / 100.

    def get_arrays_from_file(self):
        pass # these will be overwritten


class DarkParticles(Particles):
    def get_arrays_from_file(self):
        idhalo=self.file[self.snap]['ParticleData']['Dark_Halo']['ID'][:] # the ids of the halo particles
        mhalo=self.file[self.snap]['ParticleData']['Dark_Halo']['Mass'][:] / self.hh
        poshalo=self.file[self.snap]['ParticleData']['Dark_Halo']['Position'][:] / self.hh
        vhalo=self.file[self.snap]['ParticleData']['Dark_Halo']['Velocity'][:] / self.hh
        self.n=len(mhalo) # number of halo particles
        self.arrays = np.rec.array([idhalo, mhalo, poshalo[:,0], poshalo[:,1], poshalo[:,2], vhalo[:,0], \
            vhalo[:,1], vhalo[:,2]], dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8')])

class StarParticles(Particles):
    def get_arrays_from_file(self):
        age=self.file[self.snap]['ParticleData']['Star']['AGE'][:]
        idstar=self.file[self.snap]['ParticleData']['Star']['ID'][:] # the ids of the star particles
        mstar=self.file[self.snap]['ParticleData']['Star']['Mass'][:] / self.hh
        posstar=self.file[self.snap]['ParticleData']['Star']['Position'][:] / self.hh
        vstar=self.file[self.snap]['ParticleData']['Star']['Velocity'][:] / self.hh
        Zstar=self.file[self.snap]['ParticleData']['Star']['Z'][:]
        self.n=len(mstar) # number of star particles
        self.arrays = np.rec.array([idstar, mstar, posstar[:,0], posstar[:,1], posstar[:,2], \
            vstar[:,0], vstar[:,1], vstar[:,2], age, Zstar[:,0], Zstar[:,1], Zstar[:,2], \
            Zstar[:,3], Zstar[:,4], Zstar[:,5], Zstar[:,6], Zstar[:,7], Zstar[:,8], Zstar[:,9], \
            Zstar[:,10]], dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8'), ('age', '<f8'), \
            ('Z', '<f8'), ('He', '<f8'), ('C', '<f8'), ('N', '<f8'), ('O', '<f8'), ('Ne', '<f8'), \
            ('Mg', '<f8'), ('Si', '<f8'), ('S', '<f8'), ('Ca', '<f8'), ('Fe', '<f8')])

fname100 = '/Users/coralwheeler/Work/Dwarf/MetalMixing/G00/irate/G00_0_100-irate.hdf5'
snap100 = 'Snapshot00100'
dm100 = DarkParticles(fname100, snap100)
dm100.get_arrays_from_file()
print "There is/are %d dark matter particle(s)\n"%(dm100.n)
s100 = StarParticles(fname100, snap100)
s100.get_arrays_from_file()
print "There is/are %d star particle(s)\n"%(s100.n)