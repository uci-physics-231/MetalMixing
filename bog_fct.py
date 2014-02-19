from __future__ import division
import h5py
import numpy as np
import math
from matplotlib import pyplot as plt
from mayavi import mlab
from scipy import stats
import scipy


def bog_arrays(fname, initstring):
    """
    Gets simulation arrays of position, mass etc. from irate file. This is for ball of gas = BOG
    
    ARGUMENTS:
    fname: irate FULL file name
    initstring: This is either the string 'ICs' to specify that we are looking at the initial
        conditions. It could also be a snap number, with 6 digits such as '000020'.
    
    RETURNS:
    mhalo, 
    rhalo, 
    vhalo, 
    idhalo, 
    
    
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

    return idhalo,mhalo,poshalo,vhalo,nhalo,visc,rho,idgas,U,mgas,ne,nh,posgas,sfr,hsml,vgas,Z,ngas


def plot_3D(pos, stride=10, clr=(0,1,0)):
    """
    Using myavi plots the array "pos" in 3D using a stride of "stride" in color 'k'
    
    """
    pts = mlab.points3d(pos[:,0][::stride],pos[:,1][::stride],pos[:,2]\
            [::stride],mode='point',scale_factor=0,opacity=.25,color=clr)



# #Recenter based on location of ALL particles
# r = np.concatenate((rhalo, rbulge, rdisk, rgasdisk, rgashalo), axis=0)
# m = np.concatenate((mh, mb, md, mgd, mgh))
# v = np.concatenate((vhalo, vbulge, vdisk, vgasdisk, vgashalo))
# 
# #Recenter based on location of bulge only
# r = rbulge
# m = mb
# v = vbulge
# 
# centot = sigclip.sigclipm(r, v, m)
# 
# rhnew = rhalo - centot[0:3] 
# rbnew = rbulge - centot[0:3] 
# rdnew = rdisk - centot[0:3] 
# rgdnew = rgasdisk - centot[0:3]
# rghnew = rgashalo - centot[0:3]
# vhnew = vhalo - centot[3:6] 
# vbnew = vbulge - centot[3:6] 
# vdnew = vdisk - centot[3:6] 
# vgdnew = vgasdisk - centot[3:6]
# vghnew = vgashalo - centot[3:6]


    
# #recenter and recalculate position and velocity
#     poshalo,vhalo = recenter
#     rhalo=np.sqrt(poshalo[:,0]**2+poshalo[:,1]**2+poshalo[:,2]**2)

