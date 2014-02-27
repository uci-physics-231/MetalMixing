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
        self.file = h5py.File(fname,'r')
        self.fname = fname
        self.snap = snap
        self.n = 0 # will be set in child class
        self.arrays = np.ndarray(self.n) # will be set in child class
        self.hh = file['Cosmology'].attrs['HubbleParam'] / 100.

    def get_arrays_from_file(self, file, snap):
        pass # these will be overwritten

class DarkParticles(Particles):
    def get_arrays_from_file(self, file, snap):
        idhalo=file[snap]['ParticleData']['Dark_Halo']['ID'][:] # the ids of the halo particles
        mhalo=file[snap]['ParticleData']['Dark_Halo']['Mass'][:] / self.hh
        poshalo=file[snap]['ParticleData']['Dark_Halo']['Position'][:] / self.hh
        vhalo=file[snap]['ParticleData']['Dark_Halo']['Velocity'][:] / self.hh
        self.n=len(mhalo) # number of halo particles
        self.arrays = rec.array([idhalo, mhalo, poshalo[:,0], poshalo[:,1], poshalo[:,2], vhalo[:,0], \
            vhalo[:,1], vhalo[:,2]], dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8')])

class StarParticles(Particles):
    def get_arrays_from_file(self, file, snap):
        age=file[snap]['ParticleData']['Star']['AGE'][:]
        idstar=file[snap]['ParticleData']['Star']['ID'][:] # the ids of the halo particles
        mstar=file[snap]['ParticleData']['Star']['Mass'][:] / hh
        posstar=file[snap]['ParticleData']['Star']['Position'][:] / hh
        vstar=file[snap]['ParticleData']['Star']['Velocity'][:] / hh
        Zstar=file[snap]['ParticleData']['Star']['Z'][:]
        self.n=len(mstar) # number of star particles
        self.arrays = rec.array([idstar, mstar, posstar[:,0], posstar[:,1], posstar[:,2], \
            vstar[:,0], vstar[:,1], vstar[:,2], age, Zstar], dtype = [('id', '<i8'), \
            ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('velx', '<f8'), \
            ('vely', '<f8'), ('velz', '<f8'), ('age', '<f8'), ('Z', '<f8')])
        
        

def get_arrays_halo(fname, initstring):
    """
    Gets simulation arrays of position, mass etc. from irate file. This is for Dark_Halo only
    
    ARGUMENTS:
    fname: irate FULL file name
    initstring: This is either the string 'ICs' to specify that we are looking at the initial
        conditions. It could also be a snap number, with 6 digits such as '000020'.
    
    RETURNS:
    idhalo,mhalo,poshalo,vhalo,nhalo
    
    """
     
    pf=h5py.File(fname,'r') # pf is now a pointer to all the data in the file but you \
    # still don't have anything in memory

    # Get hubble value assumed to run Gadget from the file (notice that it
    # is an attribute not a folder)
    hubble=pf['Cosmology'].attrs['HubbleParam']
    hh = hubble # hh = 0.7

    #### Dark Halo
    idhalo=pf[initstring]['ParticleData']['Dark_Halo']['ID'][:] # the ids of the halo particles
    mhalo=pf[initstring]['ParticleData']['Dark_Halo']['Mass'][:] / hh
    poshalo=pf[initstring]['ParticleData']['Dark_Halo']['Position'][:] / hh
    vhalo=pf[initstring]['ParticleData']['Dark_Halo']['Velocity'][:] / hh
    nhalo=len(mhalo) # number of halo particles

    return idhalo,mhalo,poshalo,vhalo,nhalo

def get_arrays_gas(fname, initstring):
    """
    Gets simulation arrays of position, mass etc. from irate file. This is for gas
    
    ARGUMENTS:
    fname: irate FULL file name
    initstring: This is either the string 'ICs' to specify that we are looking at the initial
        conditions. It could also be a snap number, with 6 digits such as '000020'.
    
    RETURNS:
    visc,rho,idgas,U,mgas,ne,nh,posgas,sfr,hsml,vgas,Zgas,ngas
    
    """
     
    pf=h5py.File(fname,'r') # pf is now a pointer to all the data in the file but you \
    # still don't have anything in memory

    # Get hubble value assumed to run Gadget from the file (notice that it
    # is an attribute not a folder)
    hubble=pf['Cosmology'].attrs['HubbleParam']
    hh = hubble # hh = 0.7

    #### Gas
    visc=pf[initstring]['ParticleData']['Gas']['ArtificialViscosity'][:] 
    rho=pf[initstring]['ParticleData']['Gas']['Density'][:]
    idgas=pf[initstring]['ParticleData']['Gas']['ID'][:]
    U=pf[initstring]['ParticleData']['Gas']['InternalEnergy'][:]
    mgas=pf[initstring]['ParticleData']['Gas']['Mass'][:] / hh
    ne=pf[initstring]['ParticleData']['Gas']['NE'][:]
    nh=pf[initstring]['ParticleData']['Gas']['NH'][:]
    posgas=pf[initstring]['ParticleData']['Gas']['Position'][:] / hh
    sfr=pf[initstring]['ParticleData']['Gas']['SFR'][:]
    hsml=pf[initstring]['ParticleData']['Gas']['SmoothingLength'][:]
    vgas=pf[initstring]['ParticleData']['Gas']['Velocity'][:] / hh
    Z=pf[initstring]['ParticleData']['Gas']['Z'][:]
    ngas=len(mgas)

    return visc,rho,idgas,U,mgas,ne,nh,posgas,sfr,hsml,vgas,Zgas,ngas
    
def get_arrays_star(fname, initstring):
    """
    Gets simulation arrays of position, mass etc. from irate file. This is for stars only
    
    ARGUMENTS:
    fname: irate FULL file name
    initstring: This is either the string 'ICs' to specify that we are looking at the initial
        conditions. It could also be a snap number, with 6 digits such as '000020'.
    
    RETURNS:
    age,idstar,mstar,posstar,vstar,Zstar,nstar
    
    
    """
     
    pf=h5py.File(fname,'r') # pf is now a pointer to all the data in the file but you \
    # still don't have anything in memory

    # Get hubble value assumed to run Gadget from the file (notice that it
    # is an attribute not a folder)
    hubble=pf['Cosmology'].attrs['HubbleParam']
    hh = hubble # hh = 0.7

    #### Dark Halo
    age=pf[initstring]['ParticleData']['Star']['AGE'][:]
    idstar=pf[initstring]['ParticleData']['Star']['ID'][:] # the ids of the halo particles
    mstar=pf[initstring]['ParticleData']['Star']['Mass'][:] / hh
    posstar=pf[initstring]['ParticleData']['Star']['Position'][:] / hh
    vstar=pf[initstring]['ParticleData']['Star']['Velocity'][:] / hh
    Zstar=pf[initstring]['ParticleData']['Star']['Z'][:]
    nstar=len(mstar) # number of halo particles

    return age,idstar,mstar,posstar,vstar,Zstar,nstar

def get_r(pos):
    """
    Converts Nx3 position vector to distance from center
    
    """
    return np.sqrt(pos[:,0]**2+pos[:,1]**2+pos[:,2]**2)

def plot_3D(pos, stride=10, clr=(0,1,0)):
    """
    Using myavi plots the array "pos" in 3D using a stride of "stride" in color 'k'
    
    """
    pts = mlab.points3d(pos[:,0][::stride],pos[:,1][::stride],pos[:,2]\
            [::stride],mode='point',scale_factor=0,opacity=.25,color=clr)



def get_center(poshalo, posgas, posstar, vhalo, vgas, vstar, mhalo, mgas, mstar, which=3):
    """
    Recenter based on location of:
        which=0: halo particles
        which=1: gas particles
        which=2: star particles
        which=3: all particles - default
    
    """
    if which == 3: # DEFAULT
        pos = np.concatenate((poshalo, posgas, posstar), axis=0)
        m = np.concatenate((mhalo, mgas, mstar))
        v = np.concatenate((vhalo, vgas, vstar))
        
    elif which == 0:
        pos = poshalo
        m = mhalo
        v = vhalo
    
    elif which == 1:
        pos = posgas
        m = mgas
        v = vgas
    
    else:
        pos = posstar
        m = mstar
        v = vstar


    return clip.sigclipm(pos, v, m)

def recenter(pos, vel, cen):
    """
    Recenters position and velocity vectors based on new center
    
    """
    poscen = pos - cen[0:3]
    velcen = vel - cen[3:6]
    return poscen, velcen

