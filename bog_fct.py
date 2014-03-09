from __future__ import division
import h5py
import numpy as np
from numpy import rec
import math
from matplotlib import pyplot as plt
#from mayavi import mlab
import sigclip 
from scipy import stats
import scipy

# Nice plots - thick axes
plt.rc('axes',linewidth=1.5)

class Particles:
    #initialize
    def __init__(self, fname, snap):
        self.fname = fname
        self.snap = snap
        self.file = h5py.File(self.fname,'r')
        self.n = 0 # will be set in child class
        self.arrays = np.ndarray(self.n) # will be set in child class
        self.pos = [] # will be set in child class
        self.vel = [] # will be set in child class
        self.hh = self.file['Cosmology'].attrs['HubbleParam'] / 100.
        self.get_arrays_from_file()

    def get_arrays_from_file(self):
        pass # these will be overwritten

    def recenter(self):
        pass # these will be overwritten

class DarkParticles(Particles):
    # As a part of the initialization, get all arrays
    def get_arrays_from_file(self):
        idhalo=self.file[self.snap]['ParticleData']['Dark_Halo']['ID'][:] # the ids of the halo particles
        mhalo=self.file[self.snap]['ParticleData']['Dark_Halo']['Mass'][:] / self.hh
        self.pos=self.file[self.snap]['ParticleData']['Dark_Halo']['Position'][:] / self.hh
        self.vel=self.file[self.snap]['ParticleData']['Dark_Halo']['Velocity'][:] / self.hh
        rhalo = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.n=len(mhalo) # number of halo particles
        self.arrays = rec.array([idhalo, mhalo, self.pos[:,0], self.pos[:,1], self.pos[:,2], rhalo, self.vel[:,0], \
            self.vel[:,1], self.vel[:,2]], dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'), ('r', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8')])
        
    # Recenter can be passed r,m,v of all particles, or just a certain type
    # this will be decided in the execute file. Arrays must be reset with new values.
    def recenter(self,r,m,v):
        centot = sigclip.sigclipm(r, v, m)
        self.pos = self.pos - centot[0:3]
        self.vel = self.vel - centot[3:6]
        rhalo = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.arrays.x = self.pos[:,0]
        self.arrays.y = self.pos[:,1]
        self.arrays.z = self.pos[:,2]
        self.arrays.velx = self.vel[:,0]
        self.arrays.vely = self.vel[:,0]
        self.arrays.velz = self.vel[:,0]
        self.arrays.r = rhalo
        
class StarParticles(Particles):
    def get_arrays_from_file(self):
        age=self.file[self.snap]['ParticleData']['Star']['AGE'][:]
        idstar=self.file[self.snap]['ParticleData']['Star']['ID'][:] # the ids of the star particles
        mstar=self.file[self.snap]['ParticleData']['Star']['Mass'][:] / self.hh
        self.pos=self.file[self.snap]['ParticleData']['Star']['Position'][:] / self.hh
        self.vel=self.file[self.snap]['ParticleData']['Star']['Velocity'][:] / self.hh
        rstar = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        Zstar=self.file[self.snap]['ParticleData']['Star']['Z'][:]
        self.n=len(mstar) # number of star particles
        self.arrays = rec.array([idstar, mstar, self.pos[:,0], self.pos[:,1], self.pos[:,2], rstar, \
            self.vel[:,0], self.vel[:,1], self.vel[:,2], age, Zstar[:,0], Zstar[:,1], Zstar[:,2], \
            Zstar[:,3], Zstar[:,4], Zstar[:,5], Zstar[:,6], Zstar[:,7], Zstar[:,8], Zstar[:,9], \
            Zstar[:,10]], dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'),  ('r', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8'), ('age', '<f8'), \
            ('Z', '<f8'), ('He', '<f8'), ('C', '<f8'), ('N', '<f8'), ('O', '<f8'), ('Ne', '<f8'), \
            ('Mg', '<f8'), ('Si', '<f8'), ('S', '<f8'), ('Ca', '<f8'), ('Fe', '<f8')])
    
    def recenter(self,r,m,v):
        centot = sigclip.sigclipm(r, v, m)
        self.pos = self.pos - centot[0:3]
        self.vel = self.vel - centot[3:6]
        rstar = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.arrays.x = self.pos[:,0]
        self.arrays.y = self.pos[:,1]
        self.arrays.z = self.pos[:,2]
        self.arrays.velx = self.vel[:,0]
        self.arrays.vely = self.vel[:,0]
        self.arrays.velz = self.vel[:,0]
        self.arrays.r = rstar

class GasParticles(Particles):
    def get_arrays_from_file(self):
        visc=self.file[self.snap]['ParticleData']['Gas']['ArtificialViscosity'][:] 
        rho=self.file[self.snap]['ParticleData']['Gas']['Density'][:]
        idgas=self.file[self.snap]['ParticleData']['Gas']['ID'][:]
        U=self.file[self.snap]['ParticleData']['Gas']['InternalEnergy'][:]
        mgas=self.file[self.snap]['ParticleData']['Gas']['Mass'][:] / self.hh
        ne=self.file[self.snap]['ParticleData']['Gas']['NE'][:]
        nh=self.file[self.snap]['ParticleData']['Gas']['NH'][:]
        self.pos=self.file[self.snap]['ParticleData']['Gas']['Position'][:] / self.hh
        sfr=self.file[self.snap]['ParticleData']['Gas']['SFR'][:]
        hsml=self.file[self.snap]['ParticleData']['Gas']['SmoothingLength'][:]
        self.vel=self.file[self.snap]['ParticleData']['Gas']['Velocity'][:] / self.hh
        rgas = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        Zgas=self.file[self.snap]['ParticleData']['Gas']['Z'][:]
        self.n=len(mgas) # number of gas particles
        self.arrays = rec.array([idgas, mgas, self.pos[:,0], self.pos[:,1], self.pos[:,2], rgas, \
            self.vel[:,0], self.vel[:,1], self.vel[:,2], visc, rho, U, ne, nh, sfr, hsml, Zgas[:,0], \
            Zgas[:,1], Zgas[:,2], Zgas[:,3], Zgas[:,4], Zgas[:,5], Zgas[:,6], Zgas[:,7], \
            Zgas[:,8], Zgas[:,9], Zgas[:,10]], dtype = [('id', '<i8'), ('mass', '<f8'), \
            ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('r', '<f8'),  ('velx', '<f8'), ('vely', '<f8'), \
            ('velz', '<f8'), ('visc', '<f8'), ('rho', '<f8'), ('U', '<f8'), ('ne', '<f8'), \
            ('nh', '<f8'), ('sfr', '<f8'), ('hsml', '<f8'), ('Z', '<f8'), ('He', '<f8'), \
            ('C', '<f8'), ('N', '<f8'), ('O', '<f8'), ('Ne', '<f8'), ('Mg', '<f8'), \
            ('Si', '<f8'), ('S', '<f8'), ('Ca', '<f8'), ('Fe', '<f8')])
    
    def recenter(self,r,m,v):
        centot = sigclip.sigclipm(r, v, m)
        self.pos = self.pos - centot[0:3]
        self.vel = self.vel - centot[3:6]
        rgas = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.arrays.x = self.pos[:,0]
        self.arrays.y = self.pos[:,1]
        self.arrays.z = self.pos[:,2]
        self.arrays.velx = self.vel[:,0]
        self.arrays.vely = self.vel[:,0]
        self.arrays.velz = self.vel[:,0]
        self.arrays.r = rgas






# New Scrips 

  
def get_order(r):
    """
    Returns ordered indices of r
    
    """    
    return np.argsort(r) # array indices in order of increasing r

def gen_bins(rad, nbins, logbins = True):
    """
    Takes an array and a number of bins and returns the bin edges
    
    """
    # Use logarithmically spaced bins?
    if logbins == True:
        radius = np.log10(rad)
        radmax = np.amax(radius)
        radmin = np.amin(radius)
        rbinslog = np.linspace(radmin, radmax, nbins+1)
        rbins = 10 ** rbinslog
    else:
        rmax = np.amax(rad) # largest radius
        rmin = np.amin(rad) # smallest radius
        rbins = np.linspace(rmin, rmax, nbins+1) # bins in r
        
    return rbins

def plot_rho(rad, mpart, nbins, clr = 'k', logbins = True):
    """
    This function takes an ordered array of radii, an ordered array of particle mass and 
    a pre-calculated array of bins and plots the mass density (3D) as a function of radius
    * Each particle must have the same mass
    This currently plots density at left endpoint of each bin. Switch rplot to midbin to get 
    midpoint instead
    
    """
    PI = 3.14159265359
    
    rbins = gen_bins(rad, nbins)
    
    vperbin = np.ndarray(nbins) # will hold volume per bin
    mperbin = np.ndarray(nbins) # will hold mass per bin
    nperbin = np.ndarray(nbins) # will hold number of particles per bin
    midbin = np.ndarray(nbins) # will store middle of r bin
    masspart = mpart[0] # each particle has the same mass
    
    # Calculate volume per bin
    for c in range(nbins):
        vperbin[c] = (rbins[c+1] ** 3 - rbins[c] ** 3) * (4./3) * PI # volume per bin (kpc^3)
        radinbin = rad[((rad >= rbins[c]) & (rad < rbins[c+1]))]
        #massinbin = mpart[((rad >= rbins[c]) & (rad < rbins[c+1]))] # not really nec until diff masses
        nperbin[c] = np.size(radinbin) # number of particles in each bin
        mperbin[c] = nperbin[c] * masspart * 10 ** 10 # Want in units of Msun (not 10^10 Msun)
        midbin[c] = (rbins[c+1] + rbins[c]) / 2

    # Volume density = npart * mpart / vbin
    rho = mperbin / vperbin    
    
    rplot = rbins[0:-1]
    plt.figure(1, figsize=(7, 6))
    plt.plot(np.log10(rplot), np.log10(rho), color = clr, linewidth = 2)
    plt.xlim(-2,1)
    plt.ylim(2,9)
    plt.xlabel(r'$\rm log \,R (kpc)$', fontsize = 22)
    plt.ylabel(r'$\rm log \,\rho(R) (M_{\odot}/kpc^{3})$', fontsize = 22)
    plt.xticks(weight = 'bold')
    plt.yticks(weight = 'bold')
    
    return vperbin, mperbin, nperbin







# 
# 
# 
# 
# def get_arrays_halo(fname, initstring):
#     """
#     Gets simulation arrays of position, mass etc. from irate file. This is for Dark_Halo only
#     
#     ARGUMENTS:
#     fname: irate FULL file name
#     initstring: This is either the string 'ICs' to specify that we are looking at the initial
#         conditions. It could also be a snap number, with 6 digits such as '000020'.
#     
#     RETURNS:
#     idhalo,mhalo,poshalo,vhalo,nhalo
#     
#     """
#      
#     pf=h5py.File(fname,'r') # pf is now a pointer to all the data in the file but you \
#     # still don't have anything in memory
# 
#     # Get hubble value assumed to run Gadget from the file (notice that it
#     # is an attribute not a folder)
#     hubble=pf['Cosmology'].attrs['HubbleParam'] / 100.
#     hh = hubble # hh = 0.7
# 
#     #### Dark Halo
#     idhalo=pf[initstring]['ParticleData']['Dark_Halo']['ID'][:] # the ids of the halo particles
#     mhalo=pf[initstring]['ParticleData']['Dark_Halo']['Mass'][:] / hh
#     poshalo=pf[initstring]['ParticleData']['Dark_Halo']['Position'][:] / hh
#     vhalo=pf[initstring]['ParticleData']['Dark_Halo']['Velocity'][:] / hh
#     nhalo=len(mhalo) # number of halo particles
# 
#     return idhalo,mhalo,poshalo,vhalo,nhalo
# 
# def get_arrays_gas(fname, initstring):
#     """
#     Gets simulation arrays of position, mass etc. from irate file. This is for gas
#     
#     ARGUMENTS:
#     fname: irate FULL file name
#     initstring: This is either the string 'ICs' to specify that we are looking at the initial
#         conditions. It could also be a snap number, with 6 digits such as '000020'.
#     
#     RETURNS:
#     visc,rho,idgas,U,mgas,ne,nh,posgas,sfr,hsml,vgas,Zgas,ngas
#     
#     """
#      
#     pf=h5py.File(fname,'r') # pf is now a pointer to all the data in the file but you \
#     # still don't have anything in memory
# 
#     # Get hubble value assumed to run Gadget from the file (notice that it
#     # is an attribute not a folder)
#     hubble=pf['Cosmology'].attrs['HubbleParam'] / 100.
#     hh = hubble # hh = 0.7
# 
#     #### Gas
#     visc=pf[initstring]['ParticleData']['Gas']['ArtificialViscosity'][:] 
#     rho=pf[initstring]['ParticleData']['Gas']['Density'][:]
#     idgas=pf[initstring]['ParticleData']['Gas']['ID'][:]
#     U=pf[initstring]['ParticleData']['Gas']['InternalEnergy'][:]
#     mgas=pf[initstring]['ParticleData']['Gas']['Mass'][:] / hh
#     ne=pf[initstring]['ParticleData']['Gas']['NE'][:]
#     nh=pf[initstring]['ParticleData']['Gas']['NH'][:]
#     posgas=pf[initstring]['ParticleData']['Gas']['Position'][:] / hh
#     sfr=pf[initstring]['ParticleData']['Gas']['SFR'][:]
#     hsml=pf[initstring]['ParticleData']['Gas']['SmoothingLength'][:]
#     vgas=pf[initstring]['ParticleData']['Gas']['Velocity'][:] / hh
#     Z=pf[initstring]['ParticleData']['Gas']['Z'][:]
#     ngas=len(mgas)
# 
#     return visc,rho,idgas,U,mgas,ne,nh,posgas,sfr,hsml,vgas,Zgas,ngas
#     
# def get_arrays_star(fname, initstring):
#     """
#     Gets simulation arrays of position, mass etc. from irate file. This is for stars only
#     
#     ARGUMENTS:
#     fname: irate FULL file name
#     initstring: This is either the string 'ICs' to specify that we are looking at the initial
#         conditions. It could also be a snap number, with 6 digits such as '000020'.
#     
#     RETURNS:
#     age,idstar,mstar,posstar,vstar,Zstar,nstar
#     
#     
#     """
#      
#     # pf is now a pointer to all the data in the file but you \
#     # still don't have anything in memory
# 
#     # Get hubble value assumed to run Gadget from the file (notice that it
#     # is an attribute not a folder)
#     hubble=pf['Cosmology'].attrs['HubbleParam'] / 100.
#     hh = hubble # hh = 0.7
# 
#     #### Star
#     age=pf[initstring]['ParticleData']['Star']['AGE'][:]
#     idstar=pf[initstring]['ParticleData']['Star']['ID'][:] # the ids of the halo particles
#     mstar=pf[initstring]['ParticleData']['Star']['Mass'][:] / hh
#     posstar=pf[initstring]['ParticleData']['Star']['Position'][:] / hh
#     vstar=pf[initstring]['ParticleData']['Star']['Velocity'][:] / hh
#     Zstar=pf[initstring]['ParticleData']['Star']['Z'][:]
#     nstar=len(mstar) # number of halo particles
# 
#     return age,idstar,mstar,posstar,vstar,Zstar,nstar
# 
# def get_r(pos):
#     """
#     Converts Nx3 position vector to distance from center
#     
#     """
#     return np.sqrt(pos[:,0]**2+pos[:,1]**2+pos[:,2]**2)
# 
# def plot_3D(pos, stride=10, clr=(0,1,0)):
#     """
#     Using myavi plots the array "pos" in 3D using a stride of "stride" in color 'k'
#     
#     """
#     pts = mlab.points3d(pos[:,0][::stride],pos[:,1][::stride],pos[:,2]\
#             [::stride],mode='point',scale_factor=0,opacity=.25,color=clr)
# 
# 
# 
# def get_center(poshalo, posgas, posstar, vhalo, vgas, vstar, mhalo, mgas, mstar, which=3):
#     """
#     Recenter based on location of:
#         which=0: halo particles
#         which=1: gas particles
#         which=2: star particles
#         which=3: all particles - default
#     
#     """
#     if which == 3: # DEFAULT
#         pos = np.concatenate((poshalo, posgas, posstar), axis=0)
#         m = np.concatenate((mhalo, mgas, mstar))
#         v = np.concatenate((vhalo, vgas, vstar))
#         
#     elif which == 0:
#         pos = poshalo
#         m = mhalo
#         v = vhalo
#     
#     elif which == 1:
#         pos = posgas
#         m = mgas
#         v = vgas
#     
#     else:
#         pos = posstar
#         m = mstar
#         v = vstar
# 
# 
#     return clip.sigclipm(pos, v, m)
# 
# def recenter(pos, vel, cen):
#     """
#     Recenters position and velocity vectors based on new center
#     
#     """
#     poscen = pos - cen[0:3]
#     velcen = vel - cen[3:6]
#     return poscen, velcen
# 
