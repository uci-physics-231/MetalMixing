from __future__ import division
import numpy as np
from .validate import Validator,add_validator_type

def gbin_to_irate(inname,outname,snapnum,potential=False,accel=False,entropy=False,timestep=False,
    t0_name="Gas",t1_name="Dark_Halo",t2_name="Dark_Disk",t3_name="Dark_Bulge",
    t4_name="Star",t5_name="Dark_Boundary",
    ics=False,s8=None,ns=None,omegaB=None,
    lname="comoving Mpc/h",lunits=[3.08568025e24,-1,1],
    vname="(km/s)*sqrt(a)",vunits=[1e5,0,0.5],
    mname="1e10 M_sun/h",munits=[1.98892e43,-1,0]):
    """
    Reads a GADGET format file block by block (e.g. coordinate block for gas
    particles), writes the block to an IRATE formate HDF5 file, and then 
    deletes that block from memory. 

    :param inname: 
        The name of the gadget binary file to be read
    :param  outname: 
        The name of the IRATE file to be written
    :param int snapnum:  
        The number to add to 'Snapshot' that becomes the name of
        the group that particle data is added to. 
    :param bool potential: 
        True to read the Makefile enabled potential block
    :param bool accel: 
        True to read the Makefile enabled acceleration block
    :param bool entropy: 
        True to read the Makefile enabled dA/dt block
    :param bool timestep: 
        True to read the Makefile enabled timestep block 
    :param t#_name: 
        Determines the name of the group that contains the data from #
    :param bool ics:  
        True if the file being converted is an initial conditions
        file, in which case the gas density and smoothing length blocks won't
        be looked for.
    :param float s8: 
        sigma_8, for the purposes of adding it to the Cosmology group
    :param float ns: 
        n_s, for the purposes of adding it to the Cosmology group
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
    import h5py,struct,os
    from numpy import fromstring,empty
    from .core import add_cosmology,check_version
    try:
        snapnum = int(snapnum)
    except ValueError:
        raise ValueError("The number used for the snapshot is not an integer--please check and provide it by hand.")
    print "\nOpening "+inname    
    f = open(inname,'rb')
    
    #Grab the cosmology to check that it's the same
    f.seek(140)
    omegaM = struct.unpack('<d',f.read(8))[0]
    omegaL = struct.unpack('<d',f.read(8))[0]
    hubble = struct.unpack('<d',f.read(8))[0]
    f.seek(0)
    
    if os.path.isfile(outname):
        new = False
        print "Opening {0} to add data to it.".format(outname)
        irate = h5py.File(outname,'a')
        
        try:
            check_version(irate,update=False)
        except ValueError,e:
            print "\nWARNING:  ",e,"\n" 
        
        if 'Cosmology' not in irate.keys():
            add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)
        else:
            try:
                add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=False)
            except KeyError:    #Then ns, s8, or omegaB was passed, but wasn't in the file before (in which case I want to add it)
                add_cosmology(outname,s8=s8,ns=ns,omegaB=omegaB,update=True)
    else:
        new = True
        print "Creating new IRATE file {0}".format(outname)
        irate = h5py.File(outname,'w')
        check_version(irate,update=True)
        add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)


    #Create the file structure
    try:
        snap = irate.create_group("Snapshot{0:05}".format(snapnum))
        print "Created new group for Snapshot{0:05}".format(snapnum)
    except ValueError:
        snap = irate["Snapshot{0:05}".format(snapnum)]
        print "Adding data to existing group Snapshot{0:05}".format(snapnum)
        
    #Now I want to know what particles I have, so I can only create groups for those that have particles
    #First read the header so I know how many particles of each type I have
    header_size = struct.unpack('<I',f.read(4))[0]

    #number of particles of each type in this file
    nfile = struct.unpack('<6I',f.read(24)) #Number of particles in this file

    masstable = struct.unpack('<6d',f.read(48))  #masses of the particle groups
        
    a = struct.unpack('<d',f.read(8))[0]        #expansion factor
    z = struct.unpack('<d',f.read(8))[0]        #redshift
    
    if "Redshift" in snap.attrs.keys():
        if snap.attrs['Redshift'] != z:
            msg = "The existing redshift for group /Snapshot{0:05} doesn't match that in the Gadget file.  Please specify the snapshot group manually.".format(snapnum)
            raise ValueError(msg)
    else:
        snap.attrs['Redshift'] = z
    
    if "ScaleFactor" in snap.attrs.keys():
        if snap.attrs['ScaleFactor'] != a:
            msg = "The existing scale factor for group /Snapshot{0:05} doesn't match that in the Gadget file.  Please specify the snapshot group manually".format(snapnum)
            raise ValueError(msg)
    else:
        snap.attrs['ScaleFactor'] = a

    flag_sfr = struct.unpack('<i',f.read(4))[0] #star formation included?
    flag_feed = struct.unpack('<i',f.read(4))[0] #feedback included?

    ntot = struct.unpack('<6i',f.read(24))      #total number of particles in the simulation (= nfile if numfiles == 1)
        
    flag_cool = struct.unpack('<i',f.read(4))[0]  #cooling included?
    numfiles = struct.unpack('<i',f.read(4))[0]   #number of files in each snapshot
    boxsize = struct.unpack('<d',f.read(8))[0] #Size of the box, if periodic
    omega0 = struct.unpack('<d',f.read(8))[0]  #matter density at z = 0
    omegaL = struct.unpack('<d',f.read(8))[0]  #vacuum energy density at z = 0
    h = struct.unpack('<d',f.read(8))[0] #hubble parameter in units of 100 km/s/Mpc
    flag_age = struct.unpack('<i',f.read(4))[0]  #stellar age included?
    flag_metals = struct.unpack('<i',f.read(4))[0]  #use metals?
    nhighword = struct.unpack('<6i',f.read(24))   #contains the most significant word of 64-bit particle numbers (if npart > 2^32)

    flag_entropy = struct.unpack('<i',f.read(4))[0] #entropy instead of thermal energy in initial conditions?

    f.seek(264,0)   #Moves to the end of the header (and block that tells you size of header)

    #Now I want to create the simulation group and the cosmology group (if they don't exist)
    try:
        simprops = irate.create_group("SimulationProperties")
        addprops = True
    except ValueError:
        #Then it already exists, so let's check that what I was going to put there matches
        simprops = irate['SimulationProperties']
        addprops = False
        try:
            msg = "\n"
            exitflag = False
            if boxsize != simprops.attrs['Boxsize']:
                msg = msg + "Box size from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_sfr != simprops.attrs['FlagSFR']:
                msg = msg + "Flag SFR from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_feed != simprops.attrs['FlagFeedback']:
                msg = msg + "Flag Feedback from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_cool != simprops.attrs['FlagCooling']:
                msg = msg + "Flag Cooling from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_age != simprops.attrs['FlagAge']:
                msg = msg + "Flag Age from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_metals != simprops.attrs['FlagMetals']:
                msg = msg + "Flag Metals from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_entropy != simprops.attrs['FlagEntropyICs']:
                msg = msg + "Flag Entropy ICs from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if exitflag:
                raise ValueError(msg)
        except KeyError:
            raise KeyError("There is a group for simulation properties, but it's not filled with all the relevant values.  This is considered odd.")
    
    if addprops:
        simprops.attrs['Boxsize'] = boxsize
        simprops.attrs['FlagSFR'] = flag_sfr
        simprops.attrs['FlagFeedback'] = flag_feed
        simprops.attrs['FlagCooling'] = flag_cool
        simprops.attrs['FlagAge'] = flag_age
        simprops.attrs['FlagMetals'] = flag_metals
        simprops.attrs['FlagEntropyICs'] = flag_entropy
        
    #Now I'm done with the IRATE file, and instead I need to create groups inside the Snapshot group
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]

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
            msg = "\n"
            #Check that the units that are there match what I have
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
    msg = "\n"
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

    #Ok, now to read in the blocks and immediately save them in the file
    #Read in the coordinates, starting with the size of the block    
    print "Reading coordinates"
    coord_size = struct.unpack('<I',f.read(4))[0]

    if ngas > 0:    #Only try to do something for a group of particles if they exist
        #Read in the binary data and convert to the appropriate type and save in the appropriate place
        gas.create_dataset("Position",data=fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
    if nhalo > 0:
        halo.create_dataset("Position",data=fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
    if ndisk > 0:
        disk.create_dataset("Position",data=fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
    print "Read coordinates for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    if nbulge > 0:
        bulge.create_dataset("Position",data=fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
    if nstar > 0:
        star.create_dataset("Position",data=fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
    if nbndry > 0:
        bndry.create_dataset("Position",data=fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))
        
    #And read the size of the coordinate block again.
    if struct.unpack('<I',f.read(4))[0] != coord_size:
        raise IOError("The block size at the end of the coordinate block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")

    #next up is velocities.  pretty identical to the coordinates.
    vel_size = struct.unpack('<I',f.read(4))[0]
    print "Reading velocities"
    
    if ngas > 0:    #Only try to do something for a group of particles if they exist
        gas.create_dataset("Velocity",data=fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
    if nhalo > 0:
        halo.create_dataset("Velocity",data=fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
    if ndisk > 0:
        disk.create_dataset("Velocity",data=fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
    print "Read velocities for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    if nbulge > 0:
        bulge.create_dataset("Velocity",data=fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
    if nstar > 0:
        star.create_dataset("Velocity",data=fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
    if nbndry > 0:
        bndry.create_dataset("Velocity",data=fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))

    if struct.unpack('<I',f.read(4))[0] != vel_size:   #And read the size of the block again.
        raise IOError("The block size at the end of the velocity block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")

    #Next up are the particle IDs, which are unsigned integers.
    id_size = struct.unpack('<I',f.read(4))[0]
    print "Reading particle IDs"

    if ngas > 0:
        gas.create_dataset("ID",data=fromstring(f.read(4*ngas),dtype='I'))
    if nhalo > 0:
        halo.create_dataset("ID",data=fromstring(f.read(4*nhalo),dtype='I'))
    if ndisk > 0:
        disk.create_dataset("ID",data=fromstring(f.read(4*ndisk),dtype='I'))
    print "Read IDs for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    if nbulge > 0:
        bulge.create_dataset("ID",data=fromstring(f.read(4*nbulge),dtype='I'))
    if nstar > 0:
        star.create_dataset("ID",data=fromstring(f.read(4*nstar),dtype='I'))
    if nbndry > 0:
        bndry.create_dataset("ID",data=fromstring(f.read(4*nbndry),dtype='I'))

    if struct.unpack('<I',f.read(4))[0] != id_size:   #And read the size of the block again.
        raise IOError("The block size at the end of the IDs block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")

    #Now I have to do the mass block.  Do this by checking if there are particles in a block that have a zero in the mass table.    
    if (ngas > 0 and masstable[0] == 0) or (nhalo > 0 and masstable[1] == 0) or (ndisk > 0 and masstable[2] == 0) or (nbulge > 0 and masstable[3] == 0) or (nstar > 0 and masstable[4] == 0) or (nbndry > 0 and masstable[5] == 0):
    #In other words, only read the size of the mass block if there is a mass block (for any of the groups)
        mass_size = struct.unpack('<I',f.read(4))[0]
    
        if ngas > 0 and masstable[0] == 0:    #There are particles in the group, but their masses aren't in the header (so they must be in the file)
            print "Reading variable masses for gas group"
            gas.create_dataset("Mass",data=fromstring(f.read(4*ngas),dtype='f'))
        elif ngas > 0 and masstable[0] > 0:    #There are particles in the group, and their masses are in the header (so I have to fill in a variable block)
            marray = empty(ngas,dtype='f')
            marray[:] = masstable[0]
            gas.create_dataset("Mass",data=marray)
        
        if nhalo > 0 and masstable[1] == 0:
            print "Reading variable masses for halo group"
            halo.create_dataset("Mass",data=fromstring(f.read(4*nhalo),dtype='f'))
        elif nhalo > 0 and masstable[1] > 0:
            marray = empty(nhalo,dtype='f')
            marray[:] = masstable[1]
            halo.create_dataset("Mass",data=marray)
            
        if ndisk > 0 and masstable[2] == 0:
            print "Reading variable masses for disk group"
            disk.create_dataset("Mass",data=fromstring(f.read(4*ndisk),dtype='f'))
        elif ndisk > 0 and masstable[2] > 0:
            marray = empty(ndisk,dtype='f')
            marray[:] = masstable[2]
            disk.create_dataset("Mass",data=marray)
            
        if nbulge > 0 and masstable[3] == 0:
            print "Reading variable masses for bulge group"
            bulge.create_dataset("Mass",data=fromstring(f.read(4*nbulge),dtype='f'))
        elif nbulge > 0 and masstable[3] > 0:
            marray = empty(nbulge,dtype='f')
            marray[:] = masstable[3]
            bulge.create_dataset("Mass",data=marray)
            
        if nstar > 0 and masstable[4] == 0:
            print "Reading variable masses for star group"
            star.create_dataset("Mass",data=fromstring(f.read(4*nstar),dtype='f'))
        elif nstar > 0 and masstable[4] > 0:
            marray = empty(nstar,dtype='f')
            marray[:] = masstable[4]
            star.create_dataset("Mass",data=marray)
            
        if nbndry > 0 and masstable[5] == 0:
            print "Reading variable masses for boundary group"
            bndry.create_dataset("Mass",data=fromstring(f.read(4*nbndry),dtype='f'))
        elif nbndry > 0 and masstable[5] > 0:
            marray = empty(nbndry,dtype='f')
            marray[:] = masstable[5]
            bndry.create_dataset("Mass",data=marray)
        
        try:    
            del marray      #If there are no masses in the mass table, marray will have never been defined
        except NameError:
            pass
            
        if struct.unpack('<I',f.read(4))[0] != mass_size:   #And read the size of the block again.
            raise IOError("The block size at the end of the mass block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")
    
    else:       #Then all the particles in the file have their mass defined in the header, so I need to fill arrays
        if ngas > 0:
            marray = empty(ngas)
            marray.fill(masstable[0])
            gas.create_dataset("Mass",data=marray)
        if nhalo > 0:
            marray = empty(nhalo)
            marray.fill(masstable[1])
            halo.create_dataset("Mass",data=marray)
        if ndisk > 0:
            marray = empty(ndisk)
            marray.fill(masstable[2])
            disk.create_dataset("Mass",data=marray)
        if nbulge > 0:
            marray = empty(nbulge)
            marray.fill(masstable[3])
            bulge.create_dataset("Mass",data=marray)
        if nstar > 0:
            marray = empty(nstar)
            marray.fill(masstable[4])
            star.create_dataset("Mass",data=marray)
        if nbndry > 0:
            marray = empty(nbndry)
            marray.fill(masstable[5])
            bndry.create_dataset("Mass",data=marray)
            
        del marray      #marray had to have been defined if there are any particles in the file
            
    #Next up is gas specific stuff:
    if ngas > 0:
        print "Reading gas specific data."
        
        #Put all this inside the if statement because I don't want it reading block sizes for gas blocks if there is no gas data (cause then there will be no block size and I will accidentally read into other data or past the end of the file.)
        
        #Internal energy:
        u_size = struct.unpack('<I',f.read(4))[0]
        gas.create_dataset("InternalEnergy",data=fromstring(f.read(4*ngas),dtype='f'))
        if struct.unpack('<I',f.read(4))[0] != u_size:
            raise IOError("The block size at the end of the internal energy block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")
        print "Read gas internal energy."
        
        if not ics:     #These two blocks aren't required in IC files.
            #Density:
            rho_size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset("Density",data=fromstring(f.read(4*ngas),dtype='f'))  
            if struct.unpack('<I',f.read(4))[0] != rho_size:
                raise IOError("The block size at the end of the density block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")
            
            #Smoothing length:
            hsml_size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset("SmoothingLength",data=fromstring(f.read(4*ngas),dtype='f'))       
            if struct.unpack('<I',f.read(4))[0] != hsml_size:
                raise IOError("The block size at the end of the HSML block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.") 
            
            print "Read gas density and smoothing lengths."
            
        else:
            print "Skipping blocks for gas density and smoothing lengths because this is an initial conditions file."
        
    #Now for the things that have to be specifically enabled in the makefile
    if potential:
        print "Reading gravitational potentials."
        phi_size = struct.unpack('<I',f.read(4))[0]
        
        if ngas > 0:
            gas.create_dataset("Potential",data=fromstring(f.read(4*ngas),dtype='f'))     
        if nhalo > 0:
            halo.create_dataset("Potential",data=fromstring(f.read(4*nhalo),dtype='f'))
        if ndisk > 0:
            disk.create_dataset("Potential",data=fromstring(f.read(4*ndisk),dtype='f'))
        print "Read gravitational potentials for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
        if nbulge > 0:
            bulge.create_dataset("Potential",data=fromstring(f.read(4*nbulge),dtype='f'))
        if nstar > 0:
            star.create_dataset("Potential",data=fromstring(f.read(4*nstar),dtype='f'))
        if nbndry > 0:
            bndry.create_dataset("Potential",data=fromstring(f.read(4*nbndry),dtype='f'))
        
        if struct.unpack('<I',f.read(4))[0] != phi_size:
            raise IOError("The block size at the end of the gravitational potential block doesn't match that at the beginning.  This may be because the file contains extra unexpected blocks.  It is suggested that you rerun, ignoring the extra optional blocks, in the hopes that the standard blocks haven't been modified.")

    
    #Next, acceleration, which is the same as velocity:
    if accel:
        print "Reading accelerations."
        accel_size = struct.unpack('<I',f.read(4))[0]
        
        if ngas > 0:
            gas.create_dataset("Acceleration",data=fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
        if nhalo > 0:
            halo.create_dataset("Acceleration",data=fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
        if ndisk > 0:
            disk.create_dataset("Acceleration",data=fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
        print "Read accelerations for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
        if nbulge > 0:
            bulge.create_dataset("Acceleration",data=fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
        if nstar > 0:
            star.create_dataset("Acceleration",data=fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
        if nbndry > 0:
            bndry.create_dataset("Acceleration",data=fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))
        
        if struct.unpack('<I',f.read(4))[0] != accel_size:
            raise IOError("The block size at the end of the acceleration block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")
    
    if entropy and ngas > 0:    #i.e. neglect their call for entropy data if there are no gas particles
        print "Reading rate of change of entropy for gas data"
        dsdt_size = struct.unpack('<I',f.read(4))[0]
        gas.create_dataset("RateofChangeofEntropy",data=fromstring(f.read(4*ngas),dtype='f'))
        if struct.unpack('<I',f.read(4))[0] != dsdt_size:
            raise IOError("The block size at the end of the entropy block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")    
    else:
        if entropy:  print "Warning:  You specified that entropy is included, but there are no gas particles."
        
    if timestep:
        print "Reading timesteps"
        timestep_size = struct.unpack('<I',f.read(4))[0]
        
        if ngas > 0:
            gas.create_dataset("TimeStep",data=fromstring(f.read(4*ngas),dtype='f'))
        if nhalo > 0:
            halo.create_dataset("TimeStep",data=fromstring(f.read(4*nhalo),dtype='f'))
        if ndisk > 0:
            disk.create_dataset("TimeStep",data=fromstring(f.read(4*ndisk),dtype='f'))     
        print "Read timesteps for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))    
        if nbulge > 0:
            bulge.create_dataset("TimeStep",data=fromstring(f.read(4*nbulge),dtype='f')) 
        if nstar > 0:
            star.create_dataset("TimeStep",data=fromstring(f.read(4*nstar),dtype='f'))
        if nbndry > 0:
            bndry.create_dataset("TimeStep",data=fromstring(f.read(4*nbndry),dtype='f'))
        
                
        if struct.unpack('<I',f.read(4))[0] != timestep_size:
            raise IOError("The block size at the end of the timestep block doesn't match that at the beginning.  There is something wrong with the file, or my reading of it.")
    
    current_pos = f.tell()
    f.seek(0,2) #Jump to the end of the file
    if current_pos == f.tell():
        print "Read "+inname
    else:
        print "Completed reading "+inname+" but there remain {0} bytes at the end of the file unread.".format(f.tell()-current_pos)
    f.close()
    if new:  print "Created IRATE format file "+outname
    else:  print "Added data to IRATE format file "+outname
    irate.close()        

def gbin2_to_irate(inname,outname,snapnum,t0_name="Gas",t1_name="Dark_Halo",t2_name="Dark_Disk",
    t3_name="Dark_Bulge",t4_name="Star",t5_name="Dark_Boundary",
    s8=None,ns=None,omegaB=None,
    fprec=4,zfields=3,
    lname="comoving Mpc/h",lunits=[3.08568025e24,-1,1],
    vname="(km/s)*sqrt(a)",vunits=[1e5,0,0.5],
    mname="1e10 M_sun/h",munits=[1.98892e43,-1,0]):
    """
    Reads a GADGET type 2 snapshot file block by block (e.g. coordinate 
    block for gas particles), writes the block to an IRATE formate HDF5 file,
    and then deletes that block from memory. If a given block identifier isn't
    recognized, one of two things will happen:  either the script will figure
    out which particles that data belongs to and add it automatically with
    the dataset named according to the label used, or if it's not obvious
    what particles the data belongs to, the user will be given the option to
    skip it or to exit entirely.

    :param inname: 
        The name of the gadget binary file to be read
    :param  outname: 
        The name of the IRATE file to be written, 
    :param int snapnum:  
        The number to add to 'Snapshot' that becomes the name of
        the group that particle data is added to.
    :param t#_name: 
        Determines the name of the group that contains the data 
        from #
    :param float s8: 
        sigma_8, for the purposes of adding it to the Cosmology group
    :param float ns: 
        n_s, for the purposes of adding it to the Cosmology group
    :param float omegaB:
        omegaB, for the purposes of adding it to the Cosmology group
    :param fprec:
        Precision of the file. Single (4) or double (8). Assumed
        single precision. Warning! Right now changing this parameter
        to double prec (8) will probably not work but I hope it is a
        good start point.
    :param zfields:
        Fields in the Z block for each particle.
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
    import h5py,struct,os,sys
    from numpy import fromstring,empty
    from .core import add_cosmology,check_version
    try:
        snapnum = int(snapnum)
    except ValueError:
        raise ValueError("The number used for the snapshot is not an integer--please check and provide it by hand.")
    print "\nOpening "+inname    
    f = open(inname,'rb')
    #First let's find out the total length of the file
    f.seek(0,2)
    flen = f.tell()
    f.seek(0,0)
    
    #Grab the cosmology to check that it's the same
    f.seek(156)
    omegaM = struct.unpack('<d',f.read(8))[0]
    omegaL = struct.unpack('<d',f.read(8))[0]
    hubble = struct.unpack('<d',f.read(8))[0]
    f.seek(0)
    
    if os.path.isfile(outname):
        new = False
        print "Opening {0} to add data to it.".format(outname)
        irate = h5py.File(outname,'a')
        
        try:
            check_version(irate,update=False)
        except ValueError,e:
            print "\nWARNING:  ",e,"\n" 
            
        if 'Cosmology' not in irate.keys():
            add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)
        else:
            try:
                add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=False)
            except KeyError:    #Then ns, s8, or omegaB was passed, but wasn't in the file before (in which case I want to add it)
                add_cosmology(outname,s8=s8,ns=ns,omegaB=omegaB,update=True)
    else:
        new = True
        print "Creating new IRATE file {0}".format(outname)
        irate = h5py.File(outname,'w')
        check_version(irate,update=True)
        add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)
        
    #Create the file structure
    try:
        snap = irate.create_group("Snapshot{0:05}".format(snapnum))
        print "Created new group for Snapshot{0:05}".format(snapnum)
    except ValueError:
        snap = irate["Snapshot{0:05}".format(snapnum)]
        print "Adding data to existing group Snapshot{0:05}".format(snapnum)
        
    #First read the header so I know how many particles of each type I have
    #Since this is a type 2 binary, first I have to read the size of a block that tells me that the header block is next, etc.
    bsize = struct.unpack('<I',f.read(4))[0]   #This should be an 8
    if f.read(4) != 'HEAD':      #Most of the time I'm going to be using the labels to tell me what something is, but here I know that I have to read the header
        print "You're reading the wrong file type or HEAD isn't the first block.  Exiting..."
        sys.exit(1337)
    hblocksize = struct.unpack('<I',f.read(4))[0]      #So this is the size of the header block, including the numbers before and after that tell you how big it is
    assert struct.unpack('<I',f.read(4))[0] == bsize   #So this is the close of the block that includes nothing but the information about what block is next
    
    #Now I can actually read the header:
    header_size = struct.unpack('<I',f.read(4))[0]

    #number of particles of each type in this file
    nfile = struct.unpack('<6I',f.read(24)) #Number of particles in this file

    masstable = struct.unpack('<6d',f.read(48))  #masses of the particle groups
        
    a = struct.unpack('<d',f.read(8))[0]        #expansion factor
    z = struct.unpack('<d',f.read(8))[0]        #redshift
    
    if "Redshift" in snap.attrs.keys():
        if snap.attrs['Redshift'] != z:
            msg = "The existing redshift for group /Snapshot{0:05} doesn't match that in the Gadget file.  Please specify the snapshot group manually.".format(snapnum)
            raise ValueError(msg)
    else:
        snap.attrs['Redshift'] = z
    
    if "ScaleFactor" in snap.attrs.keys():
        if snap.attrs['ScaleFactor'] != a:
            msg = "The existing scale factor for group /Snapshot{0:05} doesn't match that in the Gadget file.  Please specify the snapshot group manually".format(snapnum)
            raise ValueError(msg)
    else:
        snap.attrs['ScaleFactor'] = a

    flag_sfr = struct.unpack('<i',f.read(4))[0] #star formation included?
    flag_feed = struct.unpack('<i',f.read(4))[0] #feedback included?

    ntot = struct.unpack('<6i',f.read(24))      #total number of particles in the simulation (= nfile if numfiles == 1)
        
    flag_cool = struct.unpack('<i',f.read(4))[0]  #cooling included?
    numfiles = struct.unpack('<i',f.read(4))[0]   #number of files in each snapshot
    boxsize = struct.unpack('<d',f.read(8))[0] #Size of the box, if periodic
    omega0 = struct.unpack('<d',f.read(8))[0]  #matter density at z = 0
    omegaL = struct.unpack('<d',f.read(8))[0]  #vacuum energy density at z = 0
    h = struct.unpack('<d',f.read(8))[0] #hubble parameter in units of 100 km/s/Mpc
    flag_age = struct.unpack('<i',f.read(4))[0]  #stellar age included?
    flag_metals = struct.unpack('<i',f.read(4))[0]  #use metals?
    nhighword = struct.unpack('<6i',f.read(24))   #contains the most significant word of 64-bit particle numbers (if npart > 2^32)

    flag_entropy = struct.unpack('<i',f.read(4))[0] #entropy instead of thermal energy in initial conditions?

    f.seek(280,0)   #Moves to the end of the header (and block that tells you size of header)

    #Now I want to create the simulation group and the cosmology group (if they don't exist)
    try:
        simprops = irate.create_group("SimulationProperties")
        addprops = True
    except ValueError:
        #Then it already exists, so let's check that what I was going to put there matches
        simprops = irate['SimulationProperties']
        addprops = False
        try:
            msg = "\n"
            exitflag = False
            if boxsize != simprops.attrs['Boxsize']:
                msg = msg + "Box size from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_sfr != simprops.attrs['FlagSFR']:
                msg = msg + "Flag SFR from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_feed != simprops.attrs['FlagFeedback']:
                msg = msg + "Flag Feedback from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_cool != simprops.attrs['FlagCooling']:
                msg = msg + "Flag Cooling from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_age != simprops.attrs['FlagAge']:
                msg = msg + "Flag Age from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_metals != simprops.attrs['FlagMetals']:
                msg = msg + "Flag Metals from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if flag_entropy != simprops.attrs['FlagEntropyICs']:
                msg = msg + "Flag Entropy ICs from {0} doesn't match that in {1}.\n".format(inname,outname)
                exitflag = True
            if exitflag:
                raise ValueError(msg)
        except KeyError:
            raise KeyError("There is a group for simulation properties, but it's not filled with all the relevant values.  This is considered odd.")
    
    if addprops:
        simprops.attrs['Boxsize'] = boxsize
        simprops.attrs['FlagSFR'] = flag_sfr
        simprops.attrs['FlagFeedback'] = flag_feed
        simprops.attrs['FlagCooling'] = flag_cool
        simprops.attrs['FlagAge'] = flag_age
        simprops.attrs['FlagMetals'] = flag_metals
        simprops.attrs['FlagEntropyICs'] = flag_entropy
                
    #Now I'm done with the IRATE file, and instead I need to create groups inside the Snapshot group
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]

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
            msg = "\n"
            #Check that the units that are there match what I have
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
    msg = "\n"
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
            
    posflag = False
    velflag = False
    idflag = False
    massflag = False
    
    while f.tell() < flen:
        #If I have this loop right, it'll check whether it's at the end of the file
        #If it's not, it'll read in an 8 byte block that includes a label and the size of the next block
        #Then, it'll go through and compare the labels to everything that's known, find where it goes, then read that block and save it to the appropriate tree
        bsize = struct.unpack('<I',f.read(4))
        label = f.read(4)
        btot = struct.unpack('<I',f.read(4))[0]
        assert struct.unpack('<I',f.read(4)) == bsize
        
        if label == 'POS ':
            posflag = True
            print "Reading coordinates for "+str(ngas+nhalo+ndisk+nbulge+nstar+nbndry)+" particles with "+str(fprec)+" bytes per float"
            coord_size = struct.unpack('<I',f.read(4))[0]

            if ngas > 0:    #Only try to do something for a group of particles if they exist
                #Read in the binary data and convert to the appropriate type and save in the appropriate place
                gas.create_dataset("Position",data=fromstring(f.read(3*fprec*ngas),dtype='f').reshape((-1,3)))
            if nhalo > 0:
                halo.create_dataset("Position",data=fromstring(f.read(3*fprec*nhalo),dtype='f').reshape((-1,3)))
            if ndisk > 0:
                disk.create_dataset("Position",data=fromstring(f.read(3*fprec*ndisk),dtype='f').reshape((-1,3)))
            #print "Read coordinates for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
            if nbulge > 0:
                bulge.create_dataset("Position",data=fromstring(f.read(3*fprec*nbulge),dtype='f').reshape((-1,3)))
            if nstar > 0:
                star.create_dataset("Position",data=fromstring(f.read(3*fprec*nstar),dtype='f').reshape((-1,3)))
            if nbndry > 0:
                bndry.create_dataset("Position",data=fromstring(f.read(3*fprec*nbndry),dtype='f').reshape((-1,3)))
        
            #And read the size of the coordinate block again.
            if struct.unpack('<I',f.read(4))[0] != coord_size:
                raise IOError("The block size at the end of the coordinate block doesn't match that at the beginning: "+str(coord_size)+". There is something wrong with the file, or my reading of it.")
            continue
                
        elif label == 'VEL ':
            velflag = True
            vel_size = struct.unpack('<I',f.read(4))[0]
            print "Reading velocities for "+str(ngas+nhalo+ndisk+nbulge+nstar+nbndry)+" particles with  "+str(fprec)+" bytes per float"
            
            if ngas > 0:    #Only try to do something for a group of particles if they exist
                gas.create_dataset("Velocity",data=fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
            if nhalo > 0:
                halo.create_dataset("Velocity",data=fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
            if ndisk > 0:
                disk.create_dataset("Velocity",data=fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
            #print "Read velocities for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
            if nbulge > 0:
                bulge.create_dataset("Velocity",data=fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
            if nstar > 0:
                star.create_dataset("Velocity",data=fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
            if nbndry > 0:
                bndry.create_dataset("Velocity",data=fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))

            if struct.unpack('<I',f.read(4))[0] != vel_size:   #And read the size of the block again.
                raise IOError("The block size at the end of the velocity block doesn't match that at the beginning: "+str(vel_size)+". There is something wrong with the file, or my reading of it.")
            continue
                
        elif label == 'ID  ':
            idflag = True
            id_size = struct.unpack('<I',f.read(4))[0]
            print "Reading IDs for "+str(ngas+nhalo+ndisk+nbulge+nstar+nbndry)+" particles with "+str(fprec)+" bytes per integer"

            if ngas > 0:
                gas.create_dataset("ID",data=fromstring(f.read(4*ngas),dtype='I'))
            if nhalo > 0:
                halo.create_dataset("ID",data=fromstring(f.read(4*nhalo),dtype='I'))
            if ndisk > 0:
                disk.create_dataset("ID",data=fromstring(f.read(4*ndisk),dtype='I'))
            #print "Read IDs for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
            if nbulge > 0:
                bulge.create_dataset("ID",data=fromstring(f.read(4*nbulge),dtype='I'))
            if nstar > 0:
                star.create_dataset("ID",data=fromstring(f.read(4*nstar),dtype='I'))
            if nbndry > 0:
                bndry.create_dataset("ID",data=fromstring(f.read(4*nbndry),dtype='I'))

            if struct.unpack('<I',f.read(4))[0] != id_size:   #And read the size of the block again.
                raise IOError("The block size at the end of the IDs block doesn't match that at the beginning: "+str(id_size)+".  There is something wrong with the file, or my reading of it.")
            continue
                
        elif label == 'MASS':
            massflag = True
            mass_size = struct.unpack('<I',f.read(4))[0]
            
            if ngas > 0 and masstable[0] == 0:    #There are particles in the group, but their masses aren't in the header (so they must be in the file)
                print "Reading variable masses for gas group: "+str(ngas)+" particles with "+str(fprec)+" bytes per float"
                gas.create_dataset("Mass",data=fromstring(f.read(fprec*ngas),dtype='f'))
            #Otherwise, there are either no particles in the group (in which case I don't have to worry about it) or the mass is in the header
            elif ngas > 0 and masstable[0] != 0:    #Then the mass is in the header
                marray = empty(ngas)
                marray.fill(masstable[0])
                gas.create_dataset("Mass",data=marray)
                
            if nhalo > 0 and masstable[1] == 0:
                print "Reading variable masses for halo group: "+str(nhalo)+" particles with "+str(fprec)+" bytes per float"
                halo.create_dataset("Mass",data=fromstring(f.read(fprec*nhalo),dtype='f'))
            elif nhalo > 0 and masstable[1] > 0:
                    marray = empty(nhalo)
                    marray.fill(masstable[1])
                    halo.create_dataset("Mass",data=marray)
                    
            if ndisk > 0 and masstable[2] == 0:
                print "Reading variable masses for disk group: "+str(ndisk)+" particles with "+str(fprec)+" bytes per float"
                disk.create_dataset("Mass",data=fromstring(f.read(fprec*ndisk),dtype='f'))
            elif ndisk > 0 and masstable[2] > 0:
                    marray = empty(ndisk)
                    marray.fill(masstable[2])
                    disk.create_dataset("Mass",data=marray)
                
            if nbulge > 0 and masstable[3] == 0:
                print "Reading variable masses for bulge group: "+str(nbulge)+" particles with "+str(fprec)+" bytes  per float"
                bulge.create_dataset("Mass",data=fromstring(f.read(fprec*nbulge),dtype='f'))
            elif nbulge > 0 and masstable[3] > 0:
                    marray = empty(nbulge)
                    marray.fill(masstable[3])
                    bulge.create_dataset("Mass",data=marray)
                
            if nstar > 0 and masstable[4] == 0:
                print "Reading variable masses for star group: "+str(nstar)+" particles with "+str(fprec)+" bytes per float"
                star.create_dataset("Mass",data=fromstring(f.read(fprec*nstar),dtype='f'))
            elif nstar > 0 and masstable[4] > 0:
                    marray = empty(nstar)
                    marray.fill(masstable[4])
                    star.create_dataset("Mass",data=marray)                
                
            if nbndry > 0 and masstable[5] == 0:
                print "Reading variable masses for boundary group: "+str(nbndry)+" particles with "+str(fprec)+" bytes per float"
                bndry.create_dataset("Mass",data=fromstring(f.read(fprec*nbndry),dtype='f'))
            elif nbndry > 0 and masstable[5] > 0:
                    marray = empty(nbndry)
                    marray.fill(masstable[5])
                    bndry.create_dataset("Mass",data=marray)                
                    
            try:    
                del marray      #If there are no masses in the mass table, marray will have never been defined
            except NameError:
                pass
                
            if struct.unpack('<I',f.read(4))[0] != mass_size:   #And read the size of the block again.
                raise IOError("The block size at the end of the mass block doesn't match that at the beginning: "+str(mass_size)+". There is something wrong with the file, or my reading of it.")
            continue
                
        elif label == 'U   ':     #Then it's internal energy, and I should do that
            u_size = struct.unpack('<I',f.read(4))[0]
            print "Reading gas internal energy for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            gas.create_dataset("InternalEnergy",data=fromstring(f.read(fprec*ngas),dtype='f'))
            if struct.unpack('<I',f.read(4))[0] != u_size:
                raise IOError("The block size at the end of the internal energy block doesn't match that at the beginning: "+str(u_size)+". There is something wrong with the file, or my reading of it.")
            continue
        
        elif label == "RHO ":
            print "Reading density for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            rho_size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset("Density",data=fromstring(f.read(fprec*ngas),dtype='f'))  
            if struct.unpack('<I',f.read(4))[0] != rho_size:
                raise IOError("The block size at the end of the density block doesn't match that at the beginning: "+str(rho_size)+". There is something wrong with the file, or my reading of it.")
            continue        #These continues are probably unnecessary, but better safe (and faster) then sorry
            
        elif label == "HSML":
            print "Reading smoothings lengths for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            hsml_size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset("SmoothingLength",data=fromstring(f.read(fprec*ngas),dtype='f'))       
            if struct.unpack('<I',f.read(4))[0] != hsml_size:
                raise IOError("The block size at the end of the HSML block doesn't match that at the beginning: "+str(hsml_size)+". There is something wrong with the file, or my reading of it.")
            continue
            
        #Now check if it's a standard Makefile enabled block:
        elif label == "POT ":
            print "Reading gravitational potentials for "+str(ngas+nhalo+ndisk+nbulge+nstar+nbndry)+" particles with "+str(fprec)+" bytes per float"
            phi_size = struct.unpack('<I',f.read(4))[0]
            
            if ngas > 0:
                gas.create_dataset("Potential",data=fromstring(f.read(fprec*ngas),dtype='f'))     
            if nhalo > 0:
                halo.create_dataset("Potential",data=fromstring(f.read(fprec*nhalo),dtype='f'))
            if ndisk > 0:
                disk.create_dataset("Potential",data=fromstring(f.read(fprec*ndisk),dtype='f'))
            #print "Read gravitational potentials for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
            if nbulge > 0:
                bulge.create_dataset("Potential",data=fromstring(f.read(fprec*nbulge),dtype='f'))
            if nstar > 0:
                star.create_dataset("Potential",data=fromstring(f.read(fprec*nstar),dtype='f'))
            if nbndry > 0:
                bndry.create_dataset("Potential",data=fromstring(f.read(fprec*nbndry),dtype='f'))
            
            if struct.unpack('<I',f.read(4))[0] != phi_size:
                raise IOError("The block size at the end of the gravitational potential block doesn't match that at the beginning: "+str(phi_size)+".  There is something wrong with the file, or my reading of it.")
            continue
            
        elif label == "ACCE":
            print "Reading accelerations for "+str(ngas+nhalo+ndisk+nbulge+nstar+nbndry)+" particles with "+str(fprec)+" bytes per integer"
            accel_size = struct.unpack('<I',f.read(4))[0]
            
            if ngas > 0:
                gas.create_dataset("Acceleration",data=fromstring(f.read(3*fprec*ngas),dtype='f').reshape((-1,3)))
            if nhalo > 0:
                halo.create_dataset("Acceleration",data=fromstring(f.read(3*fprec*nhalo),dtype='f').reshape((-1,3)))
            if ndisk > 0:
                disk.create_dataset("Acceleration",data=fromstring(f.read(3*fprec*ndisk),dtype='f').reshape((-1,3)))
            #print "Read accelerations for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
            if nbulge > 0:
                bulge.create_dataset("Acceleration",data=fromstring(f.read(3*fprec*nbulge),dtype='f').reshape((-1,3)))
            if nstar > 0:
                star.create_dataset("Acceleration",data=fromstring(f.read(3*fprec*nstar),dtype='f').reshape((-1,3)))
            if nbndry > 0:
                bndry.create_dataset("Acceleration",data=fromstring(f.read(3*frec*nbndry),dtype='f').reshape((-1,3)))
            
            if struct.unpack('<I',f.read(4))[0] != accel_size:
                raise IOError("The block size at the end of the acceleration block doesn't match that at the beginning: "+str(accel_size)+". There is something wrong with the file, or my reading of it.")
            continue
            
        elif label == "ENDT":
            print "Reading rate of change of entropy for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            dsdt_size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset("RateofChangeofEntropy",data=fromstring(f.read(fprec*ngas),dtype='f'))
            if struct.unpack('<I',f.read(4))[0] != dsdt_size:
                raise IOError("The block size at the end of the entropy block doesn't match that at the beginning: "+str(dsdt_size)+". There is something wrong with the file, or my reading of it.")
            continue
            
        elif label == "TSTP":
            print "Reading timesteps for "+str(ngas+nhalo+ndisk+nbulge+nstar+nbndry)+" particles with "+str(fprec)+" bytes per integer"
            timestep_size = struct.unpack('<I',f.read(4))[0]
        
            if ngas > 0:
                gas.create_dataset("TimeStep",data=fromstring(f.read(fprec*ngas),dtype='f'))
            if nhalo > 0:
                halo.create_dataset("TimeStep",data=fromstring(f.read(fprec*nhalo),dtype='f'))
            if ndisk > 0:
                disk.create_dataset("TimeStep",data=fromstring(f.read(fprec*ndisk),dtype='f'))     
            #print "Read timesteps for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))    
            if nbulge > 0:
                bulge.create_dataset("TimeStep",data=fromstring(f.read(fprec*nbulge),dtype='f')) 
            if nstar > 0:
                star.create_dataset("TimeStep",data=fromstring(f.read(fprec*nstar),dtype='f'))
            if nbndry > 0:
                bndry.create_dataset("TimeStep",data=fromstring(f.read(fprec*nbndry),dtype='f'))
                                
            if struct.unpack('<I',f.read(4))[0] != timestep_size:
                raise IOError("The block size at the end of the timestep block doesn't match that at the beginning: "+str(timesteps_size)+". There is something wrong with the file, or my reading of it.")
            continue
        
        #Ok, so if I'm here, it's a nonstandard block.  How to handle those....hmm
        #Well, I'm not going to bother renaming things; rather, I'll just use the strings as names
        #The primary difficulty lies in determining what type of particle the various blocks belong to...
        #Certain blocks I know belong to gas and/or stars, others it's not obvious
        #If a block is unrecognized, I'm going to give the option to skip or exit.
        elif label == "NE  ":   #A gas only thing
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
        
        elif label == "NH  ":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
                        
        elif label == "SFR ":   #A gas only block
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
                
        elif label == "AGE ":   #Stars only block
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            star.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*nstar),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
                        
        elif label == "Z   ":   #Gas and stars
            print "Reading "+label+" for "+str(ngas+nstar)+" (gas+stars) particles with "+str(fprec)+" bytes per float ("+str(zfields)+" fields per particle)"
            size = struct.unpack('<I',f.read(4))[0]
            if ngas > 0:  gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(zfields*fprec*ngas),dtype='f').reshape((-1,3)))
            if nstar > 0:  star.create_dataset(label.rstrip(' '),data=fromstring(f.read(zfields*fprec*nstar),dtype='f').reshape((-1,3)))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")    
            continue
        
        elif label == "BFLD":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
        
        elif label == "DBDT":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
        
        elif label == "DIVB":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")  
            continue
                
        elif label == "ABVC":  
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*gas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
                
        elif label == "COOR":    
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
        
        elif label == "CONR":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue
            
        elif label == "BFSM":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue

        elif label == "DENN":
            print "Reading "+label+" for "+str(ngas)+" (gas) particles with "+str(fprec)+" bytes per float"
            size = struct.unpack('<I',f.read(4))[0]
            gas.create_dataset(label.rstrip(' '),data=fromstring(f.read(fprec*ngas),dtype='f'))    
            if struct.unpack('<I',f.read(4))[0] != size:
                raise IOError("The block size at the end of the "+label.rstrip(' ')+" block doesn't match that at the beginning: "+str(size)+". There is something wrong with the file, or my reading of it.")
            continue

        else:
            cont = raw_input("Block labeled '"+label+"' isn't recognized by this script.  Press enter to skip it or type 'exit' to exit. ")
            if len(cont) > 0:
                if massflag == False:
                    #Then I need to fill in mass blocks for those particles that have their masses defined in the header before exiting because the mass block hasn't come up yet or doesn't exist.
                    if ngas > 0 and masstable[0] > 0:    #There are particles in the group, and their masses are in the header (so I have to fill in a variable block)
                            marray = empty(ngas)
                            marray.fill(masstable[0])
                            gas.create_dataset("Mass",data=marray)
                    if nhalo > 0 and masstable[1] > 0:
                            marray = empty(nhalo)
                            marray.fill(masstable[1])
                            halo.create_dataset("Mass",data=marray)
                    if ndisk > 0 and masstable[2] > 0:
                            marray = empty(ndisk)
                            marray.fill(masstable[2])
                            disk.create_dataset("Mass",data=marray)
                    if nbulge > 0 and masstable[3] > 0:
                            marray = empty(nbulge)
                            marray.fill(masstable[3])
                            bulge.create_dataset("Mass",data=marray)
                    if nstar > 0 and masstable[4] > 0:
                            marray = empty(nstar)
                            marray.fill(masstable[4])
                            star.create_dataset("Mass",data=marray)
                    if nbndry > 0 and masstable[5] > 0:
                            marray = empty(nbndry)
                            marray.fill(masstable[5])
                            bndry.create_dataset("Mass",data=marray)
                    try:    
                        del marray      #If there are no masses in the mass table, marray will have never been defined
                    except NameError:
                        pass
                f.close()
                irate.close()
                print "All blocks up to "+label+" from "+inname+" have been saved to "+outname+"."
                sys.exit(1)
            else:
                f.seek(btot,1)
                continue
    
    #Now I need to fill in mass blocks for those particles that have their masses defined in the header
    if massflag == False:   #Then mass wasn't anywhere in the file, so I need to create them for all the groups
        if ngas > 0 and masstable[0] > 0:    #There are particles in the group, and their masses are in the header (so I have to fill in a variable block)
                marray = empty(ngas)
                marray.fill(masstable[0])
                gas.create_dataset("Mass",data=marray)
        if nhalo > 0 and masstable[1] > 0:
                marray = empty(nhalo)
                marray.fill(masstable[1])
                halo.create_dataset("Mass",data=marray)
        if ndisk > 0 and masstable[2] > 0:
                marray = empty(ndisk)
                marray.fill(masstable[2])
                disk.create_dataset("Mass",data=marray)
        if nbulge > 0 and masstable[3] > 0:
                marray = empty(nbulge)
                marray.fill(masstable[3])
                bulge.create_dataset("Mass",data=marray)
        if nstar > 0 and masstable[4] > 0:
                marray = empty(nstar)
                marray.fill(masstable[4])
                star.create_dataset("Mass",data=marray)
        if nbndry > 0 and masstable[5] > 0:
                marray = empty(nbndry)
                marray.fill(masstable[5])
                bndry.create_dataset("Mass",data=marray)
        try:    
            del marray      #If there are no masses in the mass table, marray will have never been defined
        except NameError:
            pass

    f.close()
    irate.close()
    msg = ""
    if posflag == False:
        msg = msg + "There was no position block identified in the file; therefore, "+outname+" doesn't conform to the IRATE format specifications.\n"
    elif velflag == False:
        msg = msg +  "There was no velocity block identified in the file; therefore, "+outname+" doesn't conform to the IRATE format specifications.\n"
    elif idflag == False:
        msg = msg +  "There was no particle ID block identified in the file; therefore, "+outname+" doesn't conform to the IRATE format specifications.\n"
    else:
        msg = msg + "Finished reading type 2 file "+inname+" and saved it to "+outname
    print msg


def irate_to_gbin(inname,outname,snapnum):
    """
    Converts a single snapshot in an IRATE format file into a SnapFormat = 1
    Gadget binary file.  If there are groups with "Star" in the name, they're
    assumed to contain star particles and will be placed in Gadget group 4;
    likewise, if there are groups with "Gas" in the name, their particles will
    be placed in Gadget group 0. Dark matter particles with the lowest mass will
    be plced in group 1 and all other particles will be placed in Gadget group 5.
    
    :param inname:  
        The input IRATE format file.
    :param outname:  
        The name of the output Gadget file
    :param int snapnum:  
        The integer identifying where the data is stored; that is,
        particle data is expected to be under '/Snapshot{snapnum}/ParticleData.
    """
    import h5py,sys,os
    from numpy import empty,append
    from struct import pack
    from .core import check_version
    
    try:
        snapnum = int(snapnum)
    except ValueError:
        raise ValueError("The number used for the snapshot is not an integer.")
    
    irate = h5py.File(inname,'r')
    
    try:
        check_version(irate,update=False)
    except ValueError,e:
        print "\nWARNING:  ",e,"\n" 
    
    #Find all the data
    try:
        snap = irate['Snapshot{0:05}'.format(snapnum)]
        pdata = snap['ParticleData']
    except KeyError:
        print "There is no group /Snapshot{0:05}/ParticleData.  Please check your file format and the snapshot number you provided and try again.".format(snapnum)
        sys.exit(1337)
    
    dHiResPos = empty(0,dtype="float32").reshape((-1,3))
    dHiResVel = empty(0,dtype="float32").reshape((-1,3))
    dHiResID = empty(0,dtype="I")
    dHiResMass = empty(0,dtype="float32")    
    
    dpos = empty(0,dtype="float32").reshape((-1,3))
    dvel = empty(0,dtype="float32").reshape((-1,3))
    dID = empty(0,dtype="I")
    dmass = empty(0,dtype="float32")
    
    gpos = empty(0,dtype="float32").reshape((-1,3))
    gvel = empty(0,dtype="float32").reshape((-1,3))
    gID = empty(0,dtype="I")
    gmass = empty(0,dtype="float32")
    u = empty(0,dtype="float32")
    rho = empty(0,dtype="float32")
    hsml = empty(0,dtype="float32")    
        
    spos = empty(0,dtype="float32").reshape((-1,3))
    svel = empty(0,dtype="float32").reshape((-1,3))
    sID = empty(0,dtype="I")
    smass = empty(0,dtype="float32")

    keys = []
    masses = []
    for key in pdata.keys():
        keys = append(keys,key)
        masses = append(masses,(pdata[key]['Mass'][:]).min()) 
    hiResKey = keys[masses.argmin()]
    
    for key in pdata.keys():
        if "Gas" in key:
            print "Saving data under /Snapshot{0:05}/ParticleData/{1} in the gas group".format(snapnum,key)
            gpos = append(gpos,pdata[key]['Position'][:],axis=0)
            gvel = append(gvel,pdata[key]['Velocity'][:],axis=0)
            gID = append(gID,pdata[key]['ID'][:])
            gmass = append(gmass,pdata[key]['Mass'][:])
            u = append(u,pdata[key]['InternalEnergy'][:])
            rho = append(rho,pdata[key]['Density'][:])
            hsml = append(hsml,pdata[key]['SmoothingLength'][:])
        elif "Star" in key:
            print "Saving data under /Snapshot{0:05}/ParticleData/{1} in the star group".format(snapnum,key)
            spos = append(spos,pdata[key]['Position'][:],axis=0)
            svel = append(svel,pdata[key]['Velocity'][:],axis=0)
            sID = append(sID,pdata[key]['ID'][:])
            smass = append(smass,pdata[key]['Mass'][:])
        elif ("Dark" in key) & (key == hiResKey):
            print "Saving data under /Snapshot{0:05}/ParticleData/{1} in the Halo group".format(snapnum,key)
            dHiResPos = append(dHiResPos,pdata[key]['Position'][:],axis=0)
            dHiResVel = append(dHiResVel,pdata[key]['Velocity'][:],axis=0)
            dHiResID = append(dHiResID,pdata[key]['ID'][:])
            dHiResMass = append(dHiResMass,pdata[key]['Mass'][:])
        else:
            print "Saving data under /Snapshot{0:05}/ParticleData/{1} in the boundary group".format(snapnum,key)
            dpos = append(dpos,pdata[key]['Position'][:],axis=0)
            dvel = append(dvel,pdata[key]['Velocity'][:],axis=0)
            dID = append(dID,pdata[key]['ID'][:])
            dmass = append(dmass,pdata[key]['Mass'][:])
            
    
    ndHiRes = dHiResPos.shape[0]
    dHiResPos = dHiResPos.astype("f")
    dHiResVel = dHiResVel.astype("f")
    dHiResID = dHiResID.astype("I")
    dHiResMass = dHiResMass.astype("f")
    
    ndark = dpos.shape[0]
    dpos = dpos.astype("f")
    dvel = dvel.astype("f")
    dID = dID.astype("I")
    dmass = dmass.astype("f")

    ngas = gpos.shape[0]
    gpos = gpos.astype("f")
    gvel = gvel.astype("f")
    gID = gID.astype('I')
    gmass = gmass.astype("f")
    u = u.astype("f")
    rho = rho.astype("f")
    hsml = hsml.astype("f")
       
    nstar = spos.shape[0]
    spos = spos.astype("f")
    svel = svel.astype("f")
    sID = sID.astype("I")
    smass = smass.astype("f")
                
    #Ok, now I have all the data--just need to write it out appropriately, starting with the header
    if os.path.isfile(outname):
        print "The file you've chosen to output to exists already.  Please choose a different outname."
        sys.exit(1337)
        
    simprops = irate['SimulationProperties'].attrs
    cosmo = irate['Cosmology'].attrs
        
    out = open(outname,'w')
    out.write(pack('<I',256))
    
    out.write(pack('<6I',ngas,ndHiRes,0,0,nstar,ndark))
    out.write(pack('<6d',0,0,0,0,0,0))  #just going to include the mass for each particle individually
    out.write(pack('<d',snap.attrs['ScaleFactor']))
    out.write(pack('<d',snap.attrs['Redshift']))
    out.write(pack('<i',simprops['FlagSFR']))
    out.write(pack('<i',simprops['FlagFeedback']))
    out.write(pack('<6i',ngas,0,0,0,nstar,ndark))
    out.write(pack('<i',simprops['FlagCooling']))
    out.write(pack('<i',1))
    out.write(pack('<d',simprops['Boxsize']))
    out.write(pack('<d',cosmo['OmegaMatter']))
    out.write(pack('<d',cosmo['OmegaLambda']))
    out.write(pack('<d',cosmo['HubbleParam']))
    out.write(pack('<i',simprops['FlagAge']))
    out.write(pack('<i',simprops['FlagMetals']))
    out.write(pack('<6i',0,0,0,0,0,0))
    out.write(pack('<i',simprops['FlagEntropyICs']))
    
    left = 260 - out.tell()
    assert left > 0
    for i in range(left):
        out.write(pack('<x'))
    out.write(pack('<I',256))
    assert out.tell() == 264
    
    irate.close()
   
    #Now to do the rest of the file
    ntot = sum([ngas,ndHiRes,nstar,ndark])
    vec_size = pack('<I',12*ntot)
    scal_size = pack('<I',4*ntot)
    gas_size = pack('<I',4*ngas)
    
    out.write(vec_size)
    if ngas: out.write(gpos.tostring())
    if ndHiRes: out.write(dHiResPos.tostring())
    if nstar: out.write(spos.tostring())
    if ndark: out.write(dpos.tostring())
    out.write(vec_size)
    
    out.write(vec_size)
    if ngas: out.write(gvel.tostring())
    if ndHiRes: out.write(dHiResVel.tostring())
    if nstar: out.write(svel.tostring())
    if ndark: out.write(dvel.tostring())
    out.write(vec_size)
    
    out.write(scal_size)
    if ngas: out.write(gID.tostring())
    if ndHiRes: out.write(dHiResID.tostring())
    if nstar: out.write(sID.tostring())
    if ndark: out.write(dID.tostring())
    out.write(scal_size)
    
    out.write(scal_size)
    if ngas: out.write(gmass.tostring())
    if ndHiRes: out.write(dHiResMass.tostring())
    if nstar: out.write(smass.tostring())
    if ndark: out.write(dmass.tostring())
    out.write(scal_size)
    
    if ngas:
        out.write(gas_size)
        out.write(u.tostring())
        out.write(gas_size)
        
        out.write(gas_size)
        out.write(rho.tostring())
        out.write(gas_size)
        
        out.write(gas_size)
        out.write(hsml.tostring())
        out.write(gas_size)
        
    out.close()
    
    print "Finished writing {0} gas particles, {1} hi-res dark particles, {2} star particles, and {3} dark particles to {4}".format(ngas,ndHiRes,nstar,ndark,outname)
    

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
    from .core import add_cosmology,check_version
    print "Opening GADGET format HDF5 file"+inname
    
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
        
        try:
            check_version(irate,update=False)
        except ValueError,e:
            print "\nWARNING:  ",e,"\n" 
        
        if "Cosmology" not in irate.keys():
            add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)
        else:
            try:
                add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=False)
            except KeyError:    #Then ns, s8, or omegaB was passed, but wasn't in the file before (in which case I want to add it)
                add_cosmology(outname,s8=s8,ns=ns,omegaB=omegaB,update=True)
    else:
        new = True
        print "Creating new IRATE file {0}".format(outname)
        irate = h5py.File(outname,'w')
        check_version(irate,update=True)
        add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)


    #Create the simulation properties group in the output file
    try:
        simprops = irate.create_group("SimulationProperties")
        addprops = True
    except ValueError:
        #Then it already exists, so let's check that what I was going to put there matches
        simprops = irate['SimulationProperties']
        addprops = False
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

class GadgetSimulationPropertiesValidator(Validator):
    groupname = 'SimulationProperties'
    def validate(self,grp):
        """
        Validates that the SimulationProperties group in an IRATE file has the
        necessary properties that are typically contained in a Gadget file
        header:  Boxsize, FlagSFR, FlagFeedback, FlagCooling, FlagAge,
        FlagMetals, and FlagEntropyICs.
        
        :param grp:
            The 'SimulationProperties' :class:`~h5py.Group` to be validated
        """ 
        #Ok, so the simulation properties need the boxsize and some flags
        required = ['Boxsize','FlagSFR','FlagFeedback','FlagCooling','FlagAge','FlagMetals','FlagEntropyICs']
        for prop in required:
            if prop not in grp.attrs:
                self.invalid("Simulation properties is missing required attribute {0}".format(prop))
        
class GadgetGasValidator(Validator):
    groupname = 'Gas'
    def validate(self,grp):
        """
        Validates that a Gas group contains the required datasets for gas and
        that they're all the proper length and dimensions.  Since the standard
        Position, Velocity, Mass, and ID are checked for by the default
        validator, this checks only for InternalEnergy, Density, and
        SmoothingLength.
        
        :param grp:
            The :class:`~h5py.Group` that will be validated
        """
        #Gas group needs to have the usual datasets, plus SmoothingLength, InternalEnergy, and Density
        if 'Position' not in grp:
            self.invalid('Position is missing')
            length = 0
        else:
            length = grp['Position'].shape[0]
            self.check_units(grp['Position'])
        if 'InternalEnergy' not in grp:
            self.invalid("InternalEnergy is missing")
        else:
            self.check_units(grp['InternalEnergy'])
            #The internal energy group needs to be 1D and the same length as position
            if len(grp['InternalEnergy'].shape) != 1:
                self.invalid('InternalEnergy dataset is not 1d')
            if length != 0:     #Then it found position fine, and I want to make sure that this is as long as that
                if grp['InternalEnergy'].shape[0] != length:
                    self.invalid("InternalEnergy dataset is not the same length as Position dataset")
            else:
                #Then position wasn't found, but I want to keep going, so set this as the new length
                length = grp['InternalEnergy'].shape[0]
        generallyrequired = ['SmoothingLength','Density']
        for key in generallyrequired:
            if key not in grp:
                self.invalid('{0} is missing (this is acceptable for initial conditions only).'.format(key))
            else:
                self.check_units(grp[key])
                #They need to be 1D and the same length as those above
                if len(grp[key].shape) != 1:
                    self.invalid('{0} dataset is not 1d'.format(key))
                if length != 0:
                    #Then the length of these needs to match that length
                    if grp[key].shape[0] != length:
                        self.invalid("{0} dataset is not the same length as other datasets in that group".format(key))
                else:
                    #Then I want to set this as the new length, because apparently neither of the ones before it had a length
                    length = grp[key].shape[0]
        
def _add_validators():
    add_validator_type('gadget',[GadgetSimulationPropertiesValidator,GadgetGasValidator])  

    
class GadgetSimulationPropertiesValidator(Validator):
    groupname = 'SimulationProperties'
    def validate(self,grp):
        """
        Validates that the SimulationProperties group in an IRATE file has the
        necessary properties that are typically contained in a Gadget file
        header:  Boxsize, FlagSFR, FlagFeedback, FlagCooling, FlagAge,
        FlagMetals, and FlagEntropyICs.
        
        :param grp:
            The 'SimulationProperties' :class:`~h5py.Group` to be validated
        """ 
        #Ok, so the simulation properties need the boxsize and some flags
        required = ['Boxsize','FlagSFR','FlagFeedback','FlagCooling','FlagAge','FlagMetals','FlagEntropyICs']
        for prop in required:
            if prop not in grp.attrs:
                self.invalid("Simulation properties is missing required attribute {0}".format(prop))
        
class GadgetGasValidator(Validator):
    groupname = 'Gas'
    def validate(self,grp):
        """
        Validates that a Gas group contains the required datasets for gas and
        that they're all the proper length and dimensions.  Since the standard
        Position, Velocity, Mass, and ID are checked for by the default
        validator, this checks only for InternalEnergy, Density, and
        SmoothingLength.
        
        :param grp:
            The :class:`~h5py.Group` that will be validated
        """
        #Gas group needs to have the usual datasets, plus SmoothingLength, InternalEnergy, and Density
        if 'Position' not in grp:
            self.invalid('Position is missing')
            length = 0
        else:
            length = grp['Position'].shape[0]
            self.check_units(grp['Position'])
        if 'InternalEnergy' not in grp:
            self.invalid("InternalEnergy is missing")
        else:
            self.check_units(grp['InternalEnergy'])
            #The internal energy group needs to be 1D and the same length as position
            if len(grp['InternalEnergy'].shape) != 1:
                self.invalid('InternalEnergy dataset is not 1d')
            if length != 0:     #Then it found position fine, and I want to make sure that this is as long as that
                if grp['InternalEnergy'].shape[0] != length:
                    self.invalid("InternalEnergy dataset is not the same length as Position dataset")
            else:
                #Then position wasn't found, but I want to keep going, so set this as the new length
                length = grp['InternalEnergy'].shape[0]
        generallyrequired = ['SmoothingLength','Density']
        for key in generallyrequired:
            if key not in grp:
                self.invalid('{0} is missing (this is acceptable for initial conditions only).'.format(key))
            else:
                self.check_units(grp[key])
                #They need to be 1D and the same length as those above
                if len(grp[key].shape) != 1:
                    self.invalid('{0} dataset is not 1d'.format(key))
                if length != 0:
                    #Then the length of these needs to match that length
                    if grp[key].shape[0] != length:
                        self.invalid("{0} dataset is not the same length as other datasets in that group".format(key))
                else:
                    #Then I want to set this as the new length, because apparently neither of the ones before it had a length
                    length = grp[key].shape[0]
        
def _add_validators():
    add_validator_type('gadget',[GadgetSimulationPropertiesValidator,GadgetGasValidator])
