import numpy as np
import xarray as xr

def hdf5tods(hfile, variables=None, labels=None):
    """
    Reads data from an HDF5 into an xarray dataset
    Assumes data is organized as: group(variable) / group(snapshot) / data

    Can be used for movie.h5
    """
    import h5py as h5
    from .utils import parse_dims

    if type(hfile) is str:
        hfile = h5.File(hfile, 'r')

    DAs = dict()
    for key in hfile.keys():
        if variables is not None and key not in variables: continue
        print(key)
        dims = parse_dims(key)

        snaps, label_list = [], []
        #-----
        # Iterate through each label (snapshot) and append it
        for label in hfile[key].keys():
            if labels is not None and float(label) not in labels: continue
            snaps.append(hfile[key][label][:].T)
            label_list.append(float(label))
        #-----

        #-----
        # Put all labels together into one DataArray that holds one variable
        DA = xr.DataArray(np.array(snaps), dims=['label', *dims], 
                          coords=dict(label=label_list))
        DAs[key] = DA
        #-----

    ds = xr.Dataset(DAs)
    return ds



def read_grid(rundir='.'):
    """
    Reads important information about runs from key files:
    input.dat
    grid.h5
    grid_def
    """
    import h5py as h5
    from os.path import join
    from .utils import get_grid_centers
    run = read_input(fname=join(rundir, 'input.dat'))
    Ngrid = read_griddef(gfile=join(rundir,'grid_def'))
    gridfile = h5.File(join(rundir, 'grid.h5'), 'r')
    run = Ngrid.merge(run)

    NX, NZ = Ngrid.NX.item(), Ngrid.NZ.item()
    grids = gridfile['grids']
    X_f = np.linspace(0, run.LX, NX+1)
    Z_f = np.linspace(0, run.LZ, NZ+1)
    Y_f = grids['y'][:]

    X, Z = get_grid_centers(X_f, Z_f)
    Y = Y_f[:-1]
    run['X'] = X; run['Y'] = Y; run['Z'] = Z
    run['X_f'] = X_f; run['Y_f'] = Y_f; run['Z_f'] = Z_f

    return run



def read_input(fname):
    """ Reads some lines from input.dat file """
    from os.path import join

    ds = xr.Dataset()
    f = open(fname,'r')
    for i, line in enumerate(f):
        if i==8:
            ds['ν'], ds['β'], ds['LX'], ds['LY'], ds['LZ'] = map(float, line.split())
        elif i==23:
            ds['Ri'] = float(line.split()[0])
            ds['Pr'] = float(line.split()[1])
    return ds

def read_griddef(gfile='grid_def'):
    """ Reads files of the kind grid_def which define parameters """
    gfile = open(gfile, 'r')
    ds = xr.Dataset()
    for s in gfile:
        if s.startswith("!") or s.lower().startswith("c"): continue
        statement = s[s.find("(")+1:s.find(")")]
        param, val = statement.strip().split('=')
        ds[param] = int(val)
    return ds



def read_flowfiles(outfiles=None, rundir=None, tstep_name='Timestep', variables=None):
    """ Iterates through flow files and reads them """
    from glob import glob
    import h5py as h5
    from os.path import join

    if outfiles is None:
        outfiles = sorted(glob(join(rundir, 'out.*.h5')))
    snaps = []
    for outf in outfiles:
        print(outf)
        try:
            hfile = h5.File(outf, 'r')
        except:
            print('Skipping', outf)
            continue
        dsaux = xr.Dataset()
        for key in hfile[tstep_name].keys():
            if variables is not None and key not in variables: continue
            print(key)
            dsaux[key] = xr.DataArray(hfile[tstep_name][key][:].T)
    
        snaps.append(dsaux)
    dsout = xr.concat(snaps, dim='itime').rename(dim_0='X', dim_1='Y', dim_2='Z')
    return dsout
    

