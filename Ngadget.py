from __future__ import division
import numpy as np

def gadget_hdf5_to_irate(inname,outname,snapnum,t0_name="Gas",t1_name="Dark_Halo",
                         t2_name="Dark_Disk",t3_name="Dark_Bulge",t4_name="Star",
                         t5_name="Dark_Boundary",s8=None,ns=None,omegaB=None,
                         lname="comoving Mpc/h",lunits=np.array([3.08568025e24,-1,1]),
                         vname="(km/s)*sqrt(a)",vunits=np.array([1e5,0,0.5]),
                         mname="1e10 M_sun/h",munits=np.array([1.98892e43,-1,0])):
    """
    Transforms Gadget type 3 (HDF5) snapshots into a format that meets the
    IRATE specifications.  This will automatically check for Makefile
    enabled blocks.
        
    :param inname:  
        Gadget type 3 snapshot to convert to IRATE
    :param outname:  
        IRATE file to output
    :param int snapnum:  
        The number to add to 'Snapshot' that becomes the name of
        the group that particle data is added to.
    :param #name:  
        Name of group that particles of type # are given
    :param float s8: 
        sigma_8, for the purposes of adding it to the Cosmology group
    :param float ns: 
        n_s, for the purposes of adding it to the Cosmology group
    :param float omegaB:
        omegaB, for the purposes of adding it to the Cosmology group
    :param ___name: 
        A human readable string that defines the units for either
        length (l), velocity (v), or mass (m)
    :param ___units: 
        A three-element array.  The first entry defines the
        conversion factor to convert either length (l), velocity (v),
        or mass (m) to CGS units.  The second element is the exponent on the
        reduced Hubble Parameter that appears in that unit, and the third
        is the exponent on the scale factor that appears in that unit. 
        
    # refers to the same names as Gadget2:
    
    * 0 = gas particles
    * 1 = halo particles
    * 2 = disk particles
    * 3 = bulge particles
    * 4 = star particles
    * 5 = bndry particles
    """
    import h5py,os,sys
    from numpy import empty
    from Ncore import add_cosmology,check_version
    print "Opening GADGET format HDF5 file"+inname
    
    # Added by CRW
    testnum = int(inname[-8:-5])
    print testnum
    print snapnum
    if testnum != snapnum:
        raise NameError('Error!! Filename does not match snapnum')
        print "This occured at snapnum %d and filename %d \n"%(snapnum, testnum)
    
    fin = h5py.File(inname,'r')
    inhead = fin['Header']
    #Grab the cosmology to check that it's the same
    omegaM = inhead.attrs['Omega0']
    omegaL = inhead.attrs['OmegaLambda']
    hubble = inhead.attrs['HubbleParam']

    
    if os.path.isfile(outname):
        new = False
        irate = h5py.File(outname,'a')
        print "Opening {0} to add data to it.".format(outname)
        
#         try:
#             check_version(irate,update=False)
#         except ValueError,e:
#             print "\nWARNING:  ",e,"\n" 
        
        if "Cosmology" not in irate.keys():
            add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)
        else:
            try:
                add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=False)
            except KeyError:    #Then ns, s8, or omegaB was passed, but wasn't in the file before (in which case I want to add it)
                add_cosmology(outname,s8=s8,ns=ns,omegaB=omegaB,update=True)
                print("Key Error")
    else:
        new = True
        print "Creating new IRATE file {0}".format(outname)
        irate = h5py.File(outname,'w')
        #check_version(irate,update=True)
        add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)


    #Create the simulation properties group in the output file
    try:
        simprops = irate.create_group("SimulationProperties")
        addprops = True
        print("no simprops error")
    except ValueError:
        #Then it already exists, so let's check that what I was going to put there matches
        simprops = irate['SimulationProperties']
        addprops = False
        print("simprops error")
        try:
            msg = "\n"
            exitflag = False
            if inhead.attrs['BoxSize'] != simprops.attrs['Boxsize']:
                msg = msg + "Box size from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if inhead.attrs['Flag_Sfr'] != simprops.attrs['FlagSFR']:
                msg = msg + "Flag SFR from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if inhead.attrs['Flag_Feedback'] != simprops.attrs['FlagFeedback']:
                msg = msg + "Flag Feedback from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if inhead.attrs['Flag_Cooling'] != simprops.attrs['FlagCooling']:
                msg = msg + "Flag Cooling from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if inhead.attrs['Flag_StellarAge'] != simprops.attrs['FlagAge']:
                msg = msg + "Flag Age from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if inhead.attrs['Flag_Metals'] != simprops.attrs['FlagMetals']:
                msg = msg + "Flag Metals from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if (inhead.attrs['Flag_Entropy_ICs'] != simprops.attrs['FlagEntropyICs']).all():
                msg = msg + "Flag Entropy ICs from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if (inhead.attrs['Flag_IC_Info'] != simprops.attrs['Flag_IC_Info']).all():
                msg = msg + "Flag Entropy ICs from {0} doesn't match that in {1}.\n".format(inname,outname)
            if exitflag:
                raise ValueError(msg)
        except KeyError:
            raise KeyError("There is a group for simulation properties, but it's not filled with all the relevant values.  This is considered odd.")
            
    if addprops:
        simprops.attrs['NumFilesPerSnapshot'] = inhead.attrs['NumFilesPerSnapshot']
        simprops.attrs['Boxsize'] = inhead.attrs['BoxSize']
        simprops.attrs['FlagSFR'] = inhead.attrs['Flag_Sfr']
        simprops.attrs['FlagFeedback'] = inhead.attrs['Flag_Feedback']
        simprops.attrs['FlagCooling'] = inhead.attrs['Flag_Cooling']
        simprops.attrs['FlagAge'] = inhead.attrs['Flag_StellarAge']
        simprops.attrs['FlagMetals'] = inhead.attrs['Flag_Metals']
        try:
            simprops.attrs['FlagEntropyICs'] = inhead.attrs['Flag_Entropy_ICs']
        except KeyError:
            print("FlagEntropyICs header flag does not exist in the file")
            simprops.attrs['FlagEntropyICs'] = 0


    masses = inhead.attrs['MassTable']
    
    print"Reading Mass Table: %g %g %g %g %g %g"%(masses[0],masses[1],masses[2],masses[3],masses[4],masses[5])
    #Create the tree structure of the output file
    nfile = inhead.attrs['NumPart_ThisFile']
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]
        
    ztime=inhead.attrs['Redshift']
    atime=inhead.attrs['Time']
    try:
        snap = irate.create_group("Snapshot{0:05}".format(snapnum))
        print "Created new group for Snapshot{0:05}".format(snapnum)
    except ValueError:
        snap = irate['Snapshot{0:05}'.format(snapnum)]
        print "Adding data to existing group Snapshot{0:05}".format(snapnum)
        
    ## ADDED
    print snap
    print snapnum

    if "Redshift" in snap.attrs.keys():
        if snap.attrs['Redshift'] != ztime:
            msg = "The existing redshift for group /Snapshot{0:05} doesn't match that in the Gadget file.  Please specify the snapshot group manually.".format(snapnum)
            raise ValueError(msg)
    else:
        snap.attrs['Redshift'] = ztime
    
    if "ScaleFactor" in snap.attrs.keys():
        if snap.attrs['ScaleFactor'] != atime:
            msg = "The existing scale factor for group /Snapshot{0:05} doesn't match that in the Gadget file.  Please specify the snapshot group manually".format(snapnum)
            raise ValueError(msg)
    else:
        snap.attrs['ScaleFactor'] = atime

    ## Add Total num particles. 
    ##Must be in snapshot because ngas and nstar varies
    snap.attrs['NumPart_Total'] = inhead.attrs['NumPart_Total']
    snap.attrs['NumPart_Total_HighWord'] = inhead.attrs['NumPart_Total_HighWord']


    capbet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    try:
        pdata = snap.create_group('ParticleData')
        #Add units to pdata
        pdata.attrs['Positionunitname'] = lname
        pdata.attrs['Positionunitcgs'] = lunits
        pdata.attrs['Velocityunitname'] = vname
        pdata.attrs['Velocityunitcgs'] = vunits
        pdata.attrs['Massunitname'] = mname
        pdata.attrs['Massunitcgs'] = munits
    except ValueError:
        pdata = snap['ParticleData']
        #If 'ParticleData' already exists, check that the units that I'm adding match up (if it has units)
        try:
            #Check that the units that are there match what I have
            msg = "\n"
            exitflag = False
            if pdata.attrs['Positionunitname'] != lname or (pdata.attrs['Positionunitcgs'] != lunits).all():
                msg = msg + "Length units in existing file are different than those provided.\n"
                exitflag = True
            if pdata.attrs['Velocityunitname'] != vname or (pdata.attrs['Velocityunitcgs'] != vunits).all():
                msg = msg + "Velocity units in existing file are different than those provided.\n"
                exitflag = True
            if pdata.attrs['Massunitname'] != mname or (pdata.attrs['Massunitcgs'] != munits).all():
                msg = msg + "Mass units in existing file are different than those provided.\n"
                exitflag = True
            if exitflag:
                raise ValueError(msg)
        except KeyError:
            #add the units, since they're not there 
            pdata.attrs['Positionunitname'] = lname
            pdata.attrs['Positionunitcgs'] = lunits
            pdata.attrs['Velocityunitname'] = vname
            pdata.attrs['Velocityunitcgs'] = vunits
            pdata.attrs['Massunitname'] = mname
            pdata.attrs['Massunitcgs'] = munits
    
    exitflag = False
    if ngas > 0:
        if t0_name[0] not in capbet:  t0_name = t0_name.capitalize()
        try:
            gas = pdata.create_group(t0_name)
            print "Saving gas data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,t0_name)
        except ValueError:
            msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save gas data.\n".format(snapnum,t0_name)
            exitflag = True
        
    if nhalo > 0:
        if t1_name[0] not in capbet:  t1_name = t1_name.capitalize()
        try:
            halo = pdata.create_group(t1_name)
            print "Saving halo data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,t1_name)
        except ValueError:
            msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save halo data.\n".format(snapnum,t1_name)
            exitflag = True
        
    if ndisk > 0:
        if t2_name[0] not in capbet:  t2_name = t2_name.capitalize()
        try:
            disk = pdata.create_group(t2_name)
            print "Saving disk data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,t2_name)
        except ValueError:
            msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save disk data.\n".format(snapnum,t2_name)
            exitflag = True
        
    if nbulge > 0:
        if t3_name[0] not in capbet:  t3_name = t3_name.capitalize()
        try:
            bulge = pdata.create_group(t3_name)
            print "Saving bulge data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,t3_name)
        except ValueError:
            msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save bulge data.\n".format(snapnum,t3_name)
            exitflag = True
        
    if nstar > 0:
        if t4_name[0] not in capbet:  t4_name = t4_name.capitalize()
        try:
            star = pdata.create_group(t4_name)
            print "Saving star data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,t4_name)
        except ValueError:
            msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save star data.\n".format(snapnum,t4_name)
            exitflag = True
        
    if nbndry > 0:
        if t5_name[0] not in capbet:  t5_name = t5_name.capitalize()
        try:
            bndry = pdata.create_group(t5_name)
            print "Saving boundary data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,t5_name)
        except ValueError:
            msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save boundary data.\n".format(snapnum,t5_name)
            exitflag = True
        
    if exitflag:
        raise ValueError(msg)   
            
    if ngas > 0:
        path = fin['PartType0']
        for block in path.keys():
            if block == 'Coordinates':
                gas.create_dataset('Position',data=path[block])
            elif block == 'Velocities':
                gas.create_dataset('Velocity',data=path[block])
            elif block == 'ParticleIDs':
                gas.create_dataset('ID',data=path[block])
            elif block == 'Masses':
                gas.create_dataset('Mass',data=path[block])
            elif block == 'ArtificialViscosity':
                gas.create_dataset('ArtificialViscosity',data=path[block])
            elif block == 'DelayTime':
                gas.create_dataset('DelayTime',data=path[block])
            elif block == 'Density':
                gas.create_dataset('Density',data=path[block])
            elif block == 'ElectronAbundance':
                gas.create_dataset('NE',data=path[block])
            elif block == 'InternalEnergy':
                gas.create_dataset('InternalEnergy',data=path[block])
            elif block == 'Metallicity':
                gas.create_dataset('Z',data=path[block])
            elif block == 'NeutralHydrogenAbundance':
                gas.create_dataset('NH',data=path[block])
            elif block == 'SmoothingLength':
                gas.create_dataset('SmoothingLength',data=path[block])
            elif block == 'StarFormationRate':
                gas.create_dataset('SFR',data=path[block])
            else:            
                gas.create_dataset(block,data=path[block])
        #If there's not currently an array for the mass, I need to make one
        if 'Masses' not in path.keys():
            marray = empty(nfile[0])
            marray.fill(masses[0])
            gas.create_dataset('Mass',data=marray)
            del marray
    
    if nhalo > 0:
        path = fin['PartType1']
        for block in path.keys():
            if block == 'Coordinates':
                halo.create_dataset('Position',data=path[block])
            elif block == 'Velocities':
                halo.create_dataset('Velocity',data=path[block])
            elif block == 'ParticleIDs':
                halo.create_dataset('ID',data=path[block])
            elif block == 'Masses':
                halo.create_dataset('Mass',data=path[block])
            else:            
                halo.create_dataset(block,data=path[block])
        #If there's not currently an array for the mass, I need to make one
        if 'Masses' not in path.keys():
            marray = empty(nfile[1])
            marray.fill(masses[1])
            halo.create_dataset('Mass',data=marray)
            del marray
            
    if ndisk > 0:
        path = fin['PartType2']
        for block in path.keys():
            if block == 'Coordinates':
                disk.create_dataset('Position',data=path[block])
            elif block == 'Velocities':
                disk.create_dataset('Velocity',data=path[block])
            elif block == 'ParticleIDs':
                disk.create_dataset('ID',data=path[block])
            elif block == 'Masses':
                disk.create_dataset('Mass',data=path[block])
            else:            
                disk.create_dataset(block,data=path[block])
        #If there's not currently an array for the mass, I need to make one
        if 'Masses' not in path.keys():
            marray = empty(nfile[2])
            marray.fill(masses[2])
            disk.create_dataset('Mass',data=marray)
            del marray
            
    if nbulge > 0:
        path = fin['PartType3']
        for block in path.keys():
            if block == 'Coordinates':
                bulge.create_dataset('Position',data=path[block])
            elif block == 'Velocities':
                bulge.create_dataset('Velocity',data=path[block])
            elif block == 'ParticleIDs':
                bulge.create_dataset('ID',data=path[block])
            elif block == 'Masses':
                bulge.create_dataset('Mass',data=path[block])
            else:            
                bulge.create_dataset(block,data=path[block])
        #If there's not currently an array for the mass, I need to make one
        if 'Masses' not in path.keys():
            marray = empty(nfile[3])
            marray.fill(masses[3])
            bulge.create_dataset('Mass',data=marray)
            del marray
            
    if nstar > 0:
        path = fin['PartType4']
        for block in path.keys():
            if block == 'Coordinates':
                star.create_dataset('Position',data=path[block])
            elif block == 'Velocities':
                star.create_dataset('Velocity',data=path[block])
            elif block == 'ParticleIDs':
                star.create_dataset('ID',data=path[block])
            elif block == 'Masses':
                star.create_dataset('Mass',data=path[block])
            elif block == 'Metallicity':
                star.create_dataset('Z',data=path[block])
            elif block == 'StellarFormationTime':
                star.create_dataset('AGE',data=path[block])
            else:            
                star.create_dataset(block,data=path[block])
        #If there's not currently an array for the mass, I need to make one
        if 'Masses' not in path.keys():
            marray = empty(nfile[4])
            marray.fill(masses[4])
            star.create_dataset('Mass',data=marray)
            del marray
            
    if nbndry > 0:
        path = fin['PartType5']
        for block in path.keys():
            if block == 'Coordinates':
                bndry.create_dataset('Position',data=path[block])
            elif block == 'Velocities':
                bndry.create_dataset('Velocity',data=path[block])
            elif block == 'ParticleIDs':
                bndry.create_dataset('ID',data=path[block])
            elif block == 'Masses':
                bndry.create_dataset('Mass',data=path[block])
            else:            
                bndry.create_dataset(block,data=path[block])
        #If there's not currently an array for the mass, I need to make one
        if 'Masses' not in path.keys():
            marray = empty(nfile[5])
            marray.fill(masses[5])
            bndry.create_dataset('Mass',data=marray)
            del marray
            
    fin.close()
    irate.close()


