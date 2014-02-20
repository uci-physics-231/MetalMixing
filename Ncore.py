"""
The :mod:`irate.core` module contains the core utilities used by IRATE. These
include a variety of convinience functions for working with IRATE format files.
While they are technically part of the :mod:`irate.core` module, they are all 
also included in the base irate namespace, so the best way to import these is to 
do::

    from irate import create_irate
    
rather than importing from :mod:`irate.core`.
"""
from __future__ import division,with_statement

try:
    from collections import namedtuple
except ImportError: #py<2.6
    namedtuple = lambda *args,**kwargs:tuple

#Names of variables that should actually be imported into the root irate package
__all__ = ['add_cosmology','add_standard_cosmology','get_irate_particle_nums',
           'get_irate_catalog_halo_nums','create_irate','get_all_particles',
           'scatter_files','gather_files','get_metadata','set_metadata',
           'get_units','set_units','get_cgs_factor']



class HDF5FileOrFilenameContext(object):
    """
    This is a convinience class intended to be used as follows::
    
        with HDF5FileOrFilenameContext(fnorhdf5f,'a') as f:
            f.dosomething(whatever)
            ...
        
        #f is now closed if fnorhdf5f was a filename
            
    if `fnorhdf5f` is an hdf5 file, it will simply be bound to `f`. If it's a
    filename, it will be opened as a new :class:`h5py.File` bound to `f`, and
    will close when leaving the with block.
    
    The second argument sets the mode with which the file is opened, or if an
    :class:`h5py.File` is given, the mode is tested to see if it matches. If
    left off, the default ('a') is used for a new file, or no checking is done
    if an existing :class:`h5py.File` is used.
    
    """
    def __init__(self,forfn,mode=None):
        self.forfn = forfn
        self.mode = mode
        self.fobj = None
        
    def __enter__(self):
        import h5py
        
        if isinstance(self.forfn,basestring):
            mode = 'a' if self.mode is None else self.mode
            self.fobj = h5py.File(self.forfn,mode)
        elif isinstance(self.forfn,h5py.File):
            self.fobj = self.forfn
            if self.mode is not None and self.mode!='r' and self.fobj.mode=='r':
                raise IOError('Tried to get write access for an h5py file that is read-only')
        else:
            raise TypeError('Tried to open a non-filename/h5py.File as an hdf5 file') 
            
        return self.fobj
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.fobj is not None and self.fobj is not self.forfn:
            self.fobj.close()
            

def create_irate(fn,dark,star,gas,headername,headerdict=None,compression=None):
    """ A Convinience function to create an IRATE format file data as numpy 
    arrays. This will overwrite a currently existing file.
    
    :param str fn: The file name to use for the new file
    :param dark: 
        The data for the dark matter particles. Should be a dict with keys
        'Position', 'Velocity', 'Mass', (and possibly others) mapping to
        :class:`numpy.ndarray` arrays. Alternatively, it can be a structured
        array (i.e. dtype with names and formats) with the names of the dtype
        giving the names of the resulting datasets.
    :param star: The data for the star particles. See `dark` for format.
    :param gas: The data for the gas particles. See `dark` for format.
    :param str headername: 
        The group name to use for the header entry or None to create no header.
    :param headerdict: 
        A dictionary mapping strings to attribute values. These will be used for
        the attributes of the header entry. If None, an attribute will be added
        to indicate header content is missing.
    :param compression:
        The type of compression to use for the datasets.  Can be any of the 
        compression that h5py supports on your system ('gzip','lzf', or 'szip')
        or None to do no compression.
    
    :returns: The :class:`h5py.File` of the newly-created file (still open).
    
    :except IRATEFormatError: 
        If any of the data arrays are missing the necessary fields.
    :except ValueError: 
        If the input data is invalid.
    
    """
    import  h5py
    
    raise NotImplementedError('this function needs to be updated to the latest IRATE format')
    from .validate import validate_file
    
    hdf = h5py.File(fn,'w')
    
    #create and populate the Header
    hdr = hdf.create_group('Header')
    if headername is not None:
        hdrsgrp = hdr.create_group(headername)
        if headerdict is None:
            hdrsgrp.attrs['noheader'] = 'Header data was not provided when this file was created.'
        else:
            attrs = hdrsgrp.attrs
            for k,v in headerdict.iteritems():
                attrs[k] = v
                
    hdf.create_group('Analysis') #leave empty
    
    #now populate the data sets
    for n,d in [('Dark',dark),('Star',star),('Gas',gas)]:
        g = hdf.create_group(n)
        
        if hasattr(d,'dtype'):
            #should be a structured array
            if d.dtype.kind != 'V':
                raise ValueError('%s data is an array, but not structured'%n)
            for dn in d.dtype.names:
                if d[dn].shape[0]>0:
                    g.create_dataset(dn,data=d[dn],compression=compression)
        else:
            for k,v in d.iteritems():
                if v.shape[0]>0:
                    g.create_dataset(k,data=v,compression=compression)
    

    validate_file(hdf) #check this in case some of the entries are invalid
    
    return hdf


IrateParticleNums = namedtuple('IrateParticleNums','ndark,nstar,ngas') #:nodoc
def get_irate_particle_nums(fn,validate=True):
    """ Convinience function to determine the number of particles of each type
    in an IRATE format file.
    
    :param str fn: 
        The filename of an IRATE formatted file or a :class:`h5py.File` object.
    :param bool validate: 
        If True, the file will be validate as IRATE format. If 'strict', strict
        validation will be used.  Otherwise no validation will be performed.
    
    :returns: a 3-namedtuple (ndark,nstar,ngas)
    """
    from .validate import validate_file
    
    with HDF5FileOrFilenameContext(fn,'r') as f:
        if validate:
            validate_file(f)
            
        ns = {}
        for pname in ('Dark','Star','Gas'):
            pgrp = f[pname]
            while 'Mass' not in pgrp:
                pgrp = iter(pgrp).next()
            ns[pname] = pgrp['Mass'].shape[0]
    
    return IrateParticleNums(ns['Dark'],ns['Star'],ns['Gas'])
    
def check_version(fn,update=False):
    from . import formatversion as v
    with HDF5FileOrFilenameContext(fn) as f:
        if not update:
            if 'IRATEVersion' in f.attrs.keys():
                if f.attrs['IRATEVersion'] != v:
                    raise ValueError("The version number in {0} doesn't match the version number of these tools.".format(f.filename))
        else:
            if 'IRATEVersion' not in f.attrs.keys() or f.attrs['IRATEVersion'] < v:
                f.attrs['IRATEVersion'] = v
            
                

def get_irate_catalog_halo_nums(fn,validate=True):
    """ Convinience function to determine the number of entries in an IRATE halo
    catalog file.
    
    :param str fn: 
        The filename of an IRATE halo catalog file or a :class:`h5py.File`
        object.
    :param bool validate: 
        If True, the file will be validate as IRATE halo catalog format.
    
    :returns: an integer with the total number of halos
    """
    from .validate import validate_file
    
    with HDF5FileOrFilenameContext(fn,'r') as f:
        if validate:
            validate_file(f)
            
        return f['Catalog']['Center'].shape[0]


def get_all_particles(fn,ptype,dataname,subsample=None,validate=True):
    """ Convinience function to build an array with *all* the particles of a 
    requested type or types (i.e. all subgroups, if present, will be visited).
    
    :param str fn: 
        The filename of an IRATE formatted file or a :class:`h5py.File` object.
    :param str ptype: 
        The particle type: 'Dark','Star', 'Gas', a sequence of a mix of those
        three, or None for all particles of all types.
    :param str dataname: The dataset name to extract (e.g. 'Position', 'Mass')
    :param int subsample: 
        Gives the stride of the data - e.g. each array will be spliced to only
        return every `subsample` data point. If None, all the data will be
        returned.
    :param bool validate: 
        If True, the file will be validate as IRATE format. If 'strict', strict
        validation will be used.  Otherwise no validation will be performed.
    
    :returns: 
        A 3-tuple (data,grpidx,grpmap) where `data` is the data arrays,
        `grpidx` is an int array with the same first dimension as `data`, and
        `grpmap` is a dictionary mapping the values in `grpidx` to group names.
    
    
    :except ValueError: If a group does not have the requested `dataname`.
    
    """
    import h5py
    import numpy as np
    
    from .validate import validate_file
    
    if ptype is None:
        ptype = ['Dark','Star','Gas']
    elif isinstance(ptype,basestring):
        ptype = [ptype]
    
    with HDF5FileOrFilenameContext(fn,'r') as f:
        if validate:
            validate_file(f)
            
        grps = {}
        for pt in ptype:
            if isinstance(iter(f[pt]).next(),h5py.Group):
                for n in f[pt].keys():
                    grps[pt+'/'+n] = f[pt][n]
            else:
                grps[pt] = f[pt]
                
        datarrs = []
        grpnames = []
        for n,g in grps.iteritems():
            if dataname not in g:
                msg = 'Group {0} in file {1} has no dataset {2}'
                raise ValueError(msg.format(n,fn,dataname))
            if subsample is None:
                datarrs.append(g[dataname][:])
            else:
                datarrs.append(g[dataname][::subsample])
            grpnames.append(n)
            
        data = np.concatenate(datarrs)
        grpidx = []
        for i,a in enumerate(datarrs):
            idx = np.empty(a.shape[0])
            idx.fill(i)
            grpidx.append(idx)
        grpidx = np.concatenate(grpidx)
        grpmap = dict(enumerate(grpnames)) 
        
        return data,grpidx,grpmap
        

def add_cosmology(iratefile,omegaM=None,omegaL=None,h=None,s8=None,ns=None,omegaB=None,
                  cosmoname=None,update=True,verbose=False):
    """ Creates or updates a cosmology section of an IRATE file.
    
    This function creates an IRATE file or opens it if it already exists, and
    either adds the given cosmology-related values to it, checks that they match
    with what's in the existing file, or updates the cosmology with the given
    values.
        
    :param iratefile:  
        The name of the file to open or a :class:`h5py.File` object.
    :param float omegaM:  
        The present day matter density or None to skip this parameter.
    :param float omegaL:  
        The present day dark energy density or None to skip this parameter.
    :param float h:  
        The reduced Hubble Parameter :math:`h=H_0/100` or None to skip this 
        parameter.
    :param float s8:  
        :math:`\sigma_8`, the amplitude of the power spectrum at 8 Mpc/h, or None
        to skip this parameter.
    :param float ns:  
        The index of the primordial power spectrum or None to skip this 
        parameter.
    :param float omegaB:
        The present day baryon density or None to skip this parameter
    :param str cosmoname: 
        The name of the cosmology as a string or None to skip this value. This
        attribute is optional in the IRATE format.
    :param bool update: 
        If an existing file is given and this is True, the parameters will be
        updated with those given to this function.  Otherwise, the file will be
        *checked* to ensure the file matches. If this is False and the file does
        not exist, an exception will be raised.
    :param bool verbose:
        If True, informational messages are printed as the function does various
        actions.
        
    :raises TypeError: If the `iratefile` input is invalid
    :raises ValueError: 
        If update is False and the file's values don't match those passed into
        the function, or the file does not exist.
    :raises KeyError:
        If update is False and there are values passed that aren't in the file.
    """
    from os import path
    
    import h5py
    
    closefile = True
    
    if isinstance(iratefile,h5py.File):
        filename = iratefile.filename
        closefile = False
    elif isinstance(iratefile,basestring):
        filename = iratefile
        if path.isfile(iratefile):
            iratefile = h5py.File(filename,'a')
        else:
            if not update:
                raise ValueError("{0} is not a file - can't check that it" +\
                                 " matches supplied parameters".format(filename))
            if verbose:
                print  "Creating IRATE file {0}".format(filename)
            iratefile = h5py.File(filename,'w')
    else:
        msg = '{0} is not a valid file name nor an h5py.File'
        raise TypeError(msg.format(iratefile))
        
    try:
        if 'Cosmology' in iratefile:
            cgrp = iratefile['Cosmology']
        else:
            if not update:
                msg = "{0} does not have a Cosmology section - can't " +\
                        "check that it matches supplied parameters"
                raise ValueError(msg.format(filename))
            if verbose:
                print "Adding cosmology section to {0}.".format(filename)   
            cgrp = iratefile.create_group('Cosmology')
        
        attrnametovar = {'OmegaMatter':None if omegaM is None else float(omegaM),
                         'OmegaLambda':None if omegaL is None else float(omegaL),
                         'HubbleParam':None if h is None else float(h) if h > 1.0 else h*100.0,
                         'PowerSpectrumIndex':None if ns is None else float(ns),
                         'sigma_8':None if s8 is None else float(s8),
                         'OmegaBaryon':None if omegaB is None else float(omegaB),
                         'Name':None if cosmoname is None else str(cosmoname)}
        
        if update:
            if verbose:
                msg = 'Updating {0} with provided cosmological parameters'
                print msg.format(filename)
            for k,v in attrnametovar.iteritems():
                if v is not None:
                    cgrp.attrs[k] = v
        else:
            if verbose:
                print 'Checking cosmological parameters in {0}'.format(filename)
            for k,v in attrnametovar.iteritems():
                if v is not None:
                    if k not in cgrp.attrs:
                        msg = 'File {0} is missing {1} cosmological parameter'
                        raise KeyError(msg.format(filename,k))
                    if cgrp.attrs[k] != v:
                        msg = 'Cosmological parameter {0} ({2}) does not match ' +\
                              'provided value in file {1} ({3})'
                        raise ValueError(msg.format(k,filename,v,cgrp.attrs[k]))
    finally:
        if closefile:
            iratefile.close()
            
_standard_cosmos= {
'WMAP7':dict(omegaM=.27,omegaL=.73,h=.71,s8=.8,ns=.96),
'WMAP5':dict(omegaM=.26,omegaL=.74,h=.72,s8=.8,ns=.96),
'WMAP3':dict(omegaM=.24,omegaL=.76,h=.73,s8=.76,ns=.96),
'WMAP1':dict(omegaM=.27,omegaL=.73,h=.72,s8=.9,ns=.99)}

def add_standard_cosmology(iratefile,cname,update=True,verbose=False):
    """ Adds or updates a cosmology from the name of a standard cosmology.
    
    The `cname` parameter gives the name of the cosmology to be used from the 
    list below, while all other parameters are the same as for the
    :func:`add_cosmology` function.
    
    * 'WMAP7'
    
        LCDM Cosmological parameters for the 7-year WMAP data (Komatsu et al. 
        2011). For details, see 
        http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm
        
    * 'WMAP5'
    
        LCDM Cosmological parameters for the 5-year WMAP data (Komatsu et al. 
        2009). For details, see 
        http://lambda.gsfc.nasa.gov/product/map/dr3/params/lcdm_sz_lens_wmap7.cfm
        
    * 'WMAP3'
    
        LCDM Cosmological parameters for the 3-year WMAP data (Spergel et al. 
        2007). For details, see 
        http://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_sz_lens_wmap7.cfm
        
    * 'WMAP1'
    
        LCDM Cosmological parameters for the 1-year WMAP data. For details, see 
        Spergel et al. 2003.
    
    """
    args = [iratefile]
    kwargs = _standard_cosmos[cname].copy()
    kwargs['cosmoname'] = cname
    kwargs['update'] = update
    kwargs['verbose'] = verbose
    
    return add_cosmology(*args,**kwargs)
    
    
def get_metadata(dataset,metadataname):
    """ Gets a metadata value from an IRATE file dataset.
    
    :param dataset: An :class:`~h5py.Dataset` object in an IRATE file.
    :param str metadataname:  The name of the meteadata field requested.
    :returns: 
        The value of the metadata requested, or None if the metadata with the
        requested name does not exist.
    
    :raises TypeError: If the `dataset` is not a dataset
    """
    import h5py
    
    if not isinstance(dataset,h5py.Dataset):
        raise TypeError('input is not a dataset')
        
    if metadataname in dataset.attrs:
        return dataset.attrs[metadataname]
        
    nm = dataset.name + metadataname
    curr = dataset
    while curr!=curr.parent:
        if nm in curr.attrs:
            return curr.attrs[nm]
        curr = curr.parent
    else:
        return None
    
class MetadataHiddenWarning(Warning):
    """
    This warning indicates that metadata is being set if 
    """
    
def set_metadata(dataset,metadataname,value,group=None):
    """ Sets a metadata value for an IRATE file dataset.
    
    :param dataset: An :class:`~h5py.Dataset` object in an IRATE file.
    :param str metadataname:  The name of the meteadata field requested.
    :param value:  The name of the meteadata field requested.
    :param group:  
        The group this metadata should be set for, or None.  
        
    .. note:: 
        Any other datasets with the same as this one below this group in the
        heirarchy will be assigned the same units (unless that dataset overrides
        them itself).
    
    :raises TypeError: If `dataset` is not a dataset.
    :raises ValueError: If `group` is not a parent of `dataset`.
    
    .. note::
        A :class:`MetadataHiddenWarning` will be issued as a warning if the
        metadata for this dataset has already been set somewhere lower down in
        the heirarchy than the requested `group`.  The result of this is that 
        this function will not alter the resulting units for the `dataset` 
        because it will be overshadowed by the metadata setting further down.
        (see :mod:`warnings` for an explanation of warnings)
    """
    import h5py
    
    from warnings import warn
    
    if not isinstance(dataset,h5py.Dataset):
        raise TypeError('input is not a dataset')
    
    if group is None:
        mdgrp = dataset
    else:
        mdgrp = group
        if metadataname in dataset.attrs:
            msg = 'Metadata {0} already set in {1}'
            warn(msg.format(metadataname,dataset),MetadataHiddenWarning)
        metadataname = dataset.name+metadataname
    
        curr = dataset
        while (curr!=curr.parent):
            if curr==group:
                break
            elif metadataname in curr.attrs:
                msg = 'Metadata {0} already set in {1}'
                warn(msg.format(metadataname,curr),MetadataHiddenWarning)
                
            curr = curr.parent
        else:
            msg = 'The given unitgroup {0} is not a parent of the dataset {1}'
            raise ValueError(msg.format(unitgroup,dataset))
        
    mdgrp.attrs[metadataname] = value

def get_units(dataset):
    """ Gets the units for a dataset in an IRATE file.
    
    Because units are typically expressed in factors of the hubble constant and
    may or may not be comoving or otherwise varying with scale factor, two 
    additional factors are needed beyond the raw unit conversion.  To use these
    together, the appropriate factor to multiple the dataset by to get to 
    physical units for an assumed hubble parameter `h` and scale factor `a` is:
    `tocgs`*h^`hubbleexponent`*a^`scalefactorexponent`  . Hence, for example, 
    comoving Mpc/h would have `tocgs` =3.0857e24 , `hubbleexponent` =-1, and 
    `scalefactorexponent` =1 .
    
    :param dataset: A :class:`~h5py.Dataset` object in an IRATE file
    :returns: 
        `name,tocgs,hubbleexponent,scalefactorexponent`. `name` is the 
        name of the unit for this dataset, `tocgs` is a factor to multiply the
        dataset by to get cgs units, `hubleexponent` is the exponent of the 
        reduced hubble parameter (h=H0/100) for this unit, and 
        `scalefactorexponent` is the exponent describing how this unit varies
        with the scale factor.  Alternatively, if no unit is found, ``None``
        is returned.
        
    .. note:: 
        I f you are not working with a file where you know for sure what the
        units are, is important to check whether or not the unit is None, 
        because the IRATE format does not require units for all datasets, as 
        some quantities are dimensionless.  The easiest way to do this is::
        
            res = get_units(mydataset)
            if res is None:
                ... do something if it is unitless ...
            else:
                nm,tocgs,hexp,sfexp = res
                ... do something with the units ...
    
    :raises TypeError: If `dataset` is not a dataset
    """
    nm = get_metadata(dataset,'unitname')
    if nm is None:
        return None
    else:
        tocgs,hexp,sfexp = get_metadata(dataset,'unitcgs')
        return nm,tocgs,hexp,sfexp

def get_cgs_factor(dataset,h=.7,a=1):
    """ Retrieves the factor to multiply by the dataset to convert to physical 
    CGS units at a given scale factor and hubble constant.
    
    
    :param dataset: An :class:`~h5py.Dataset` for which to read the units.
    :param h: Reduced hubble parameter to assume for the conversion.
    :param a: Scale factor to assume for the conversion = 1/(z+1).
    
    :raises TypeError: If `dataset` is not a dataset
    :raises ValueError: If the dataset is dimensionless
    """
    res = get_metadata(dataset,'unitcgs')
    if res is None:
        raise ValueError('dataset {0} does not have units'.format(dataset))
    tocgs,hexp,sfexp = res
    return tocgs * h**hexp * a**sfexp
    

def set_units(dataset,unitname,tocgs,hubbleexponent,scalefactorexponent,
              unitgroup=None):
    """
    Sets the units of a dataset to the provided values.
    
    Because units are typically expressed in factors of the hubble constant and
    may or may not be comoving or otherwise varying with scale factor, two 
    additional factors are needed beyond the raw unit conversion.  To use these
    together, the appropriate factor to multiple the dataset by to get to 
    physical units for an assumed hubble parameter `h` and scale factor `a` is:
    `tocgs`*h^`hubbleexponent`*a^`scalefactorexponent`  . Hence, for example, 
    comoving Mpc/h would have `tocgs` =3.0857e24 , `hubbleexponent` =-1, and 
    `scalefactorexponent` =1 .
    
    :param dataset: The dataset for which to set the units.
    :param str unitname: A human-readable name of the unit.
    :param str tocgs: 
        A numerical factor to multiply the dataset by to convert to cgs units.
    :param str hubbleexponent: 
        The exponent of the hubble parameter for this unit.
    :param str scalefactorexponent: 
        The exponent of the scale factor for this unit.
    
    :raise TypeError: If an input is not the correct type
    """
    
    if not isinstance(unitname,basestring):
        raise TypeError('unitname is not a string')
    try:
        tocgs = float(tocgs)
        he = float(hubbleexponent)
        sfe = float(scalefactorexponent)
    except ValueError:
        raise TypeError('tocgs, hubbleexponet, or scalefactorexponent are not float-compatible')
        
    set_metadata(dataset,'unitname',unitname,unitgroup)
    set_metadata(dataset,'unitcgs',(tocgs,he,sfe),unitgroup)
    
                
def scatter_files(infile,outfile=None,splitgroups=False):
    """ Splits a single IRATE file into separate files for each dataset.
    
    :param infile: A file name or an `h5py.File` object to be scattered.
    :param outfile: 
        The filename or an `h5py.File` object to use as the new base IRATE format
        file or None to overwrite the old name. Other output files will be placed
        in the same directory as this file.
    :param bool splitgroups: If True, also split groups into separate files.
    
    :returns: 
        A list of filenames that were created with the base name as the first. 
    """
    raise NotImplementedError

def _recursive_copy(ingrp,outgrp):
    """
    This function is used to recursively copy over all parts of one HDF5 file
    to another, returning a list of external files.
    """
    extfiles = []
    for k,v in ingrp.attrs.iteritems():
        outgrp.attrs[k] = v 
    for k in ingrp:
        v = ingrp.get(k,getlink=True)
        if hasattr(v,'path') and hasattr(v,'filename'):
            extfiles.append(v.filename)
            v = ingrp.get(k)
        if hasattr(v,'keys'):
            #group
            newv = outgrp.create_group(k)
            extfiles.extend(_recursive_copy(v,newv))
        else:
            #dataset
            kwds = dict(name=k,data=v.value)
            #copy over dataset properties
            for kkw in ('chunks','maxshape','compression','compression_opts',
                        'shuffle','fletcher32'):
                kwds[kkw] = getattr(v,kkw)
            #copy over dataset
            newv = outgrp.create_dataset(**kwds)
            #copy over dataset attributes
            for ki,vi in v.attrs.iteritems():
                newv.attrs[ki] = vi
    return extfiles

def gather_files(infile,outfile=None):
    """ Combine all external file that are part of a master IRATE file into a
    single monolithic file.
    
    .. note::
        This function can be driven from the command line via the ``irategather``
        script.  See :doc:`scripts` for details.
    
    :param infile: A file name or an `h5py.File` object to be gathered.
    :param outfile: 
        The filename or an `h5py.File` object to use as the new IRATE format file
        or None to overwrite the old name.
    
    :returns: 
        (outfilename,filelist) where `outfilename` is the name of the output
        file, and `filelist` is a list of all the files gathered together to
        create this file.
    """
    from os.path import exists
    
    if outfile is None:
        if not isinstance(infile,basestring):
            raise IOError('cannot overwrite open hdf5 file provided as infile')
        outfile = infile+'.tmp'
        i = 1
        while exists(outfile):
            outfile = infile+'.tmp'+str(i)
        moveouttoin = True
    else:
        moveouttoin = False
    
    gatherlist = []    
    with HDF5FileOrFilenameContext(infile) as fin:
        infilename = fin.filename
        gatherlist.append(infilename)
        with HDF5FileOrFilenameContext(outfile) as fout:
            outfilename = fout.filename
            gatherlist.extend(set(_recursive_copy(fin,fout)))
    if moveouttoin:
        from os import remove
        from shutil import move
        remove(infilename)
        move(outfilename,infilename)
        resfile = infilename
    else:
        resfile = outfilename
        
    return resfile,gatherlist
