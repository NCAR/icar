#!/usr/bin/env python
import os
import glob
import time
# import multiprocessing as mp

import numpy as np
import xarray as xr

import sys

pool = None

# This should be an input, this is the search string that is assumed to match
# the output files to be aggregated.
file_search = "icar_out_{ens}_*"


def load_file(file_name):
    '''Load a netcdf dataset into memory'''
    return xr.open_dataset(file_name).load()


def get_dims(dataset, section="d"):
    '''Get the global attributes defining the domain, memory, or tile space'''
    results = []
    for axis in ["i","j","k"]:
        for position in ["s","e"]:
            results.append(int(dataset.attrs[axis + section + position]))
    return results

def get_dim_offset(dims):
    '''Return x_offset, y_offset
    For the staggered dims, offset=1, otherwise offset=0'''
    x_off = 0
    if 'lon_u' in dims: x_off = 1

    y_off = 0
    if 'lat_v' in dims: y_off = 1

    return x_off, y_off

def set_up_dataset(d):
    '''Create a dataset to cover the entire domain with the variables present in d

    d : an input dataset covering part of the domain
    d must have global attributes ids, ide, jds, jde, kds, kde that define the full domain

    A new dataset is created with all the variables+attributes in d covering the full domain
    '''
    ids, ide, jds, jde, kds, kde = get_dims(d, section='d')
    nx = ide - ids + 1
    ny = jde - jds + 1
    nz = kde - kds + 1

    data_vars = dict()

    for v in d.variables:
        coords = [c for c in d[v].coords]
        dims   = d[v].dims
        name   = d[v].name
        attrs  = d[v].attrs

        x_off, y_off = get_dim_offset(dims)

        if len(dims) == 1:
            nt = d.sizes[dims[0]]
            data = np.zeros((nt))
        if len(dims) == 2:
            data = np.zeros((ny + y_off, nx + x_off))
        if len(dims) == 3:
            data = np.zeros((d.sizes[dims[0]], ny + y_off, nx + x_off))
        if len(dims) == 4:
            nt = d.sizes[dims[0]]
            nz = d.sizes[dims[1]]
            data = np.zeros((nt, nz, ny + y_off, nx + x_off))

        # print(name, data.shape, dims, attrs)
        data_vars[v] = xr.DataArray(data.astype(d[v].dtype), dims=dims, name=name, attrs=attrs)#, coords=coords)

    ds = xr.Dataset(data_vars, attrs=d.attrs)
    ds.encoding = d.encoding
    ds["time"] = d["time"]
    return ds.set_coords([c for c in d.coords])


def agg_file(first_file, verbose=True):
    '''Aggregated all files that come from the same time step as first_file

    first_file should have _001_ in the filename somewhere.  This will be replaced
    with * to search for all matching files from this date. Once files are found, a
    dataset containing the entire domain is created and the data from each file are
    added to the master dataset.

    Result: aggregated dataset is written to a netcdf file'''

    if verbose:print(first_file)
    date_search = first_file.replace("000001_","*")
    outputfile = first_file.replace("000001_","_").replace("__","_")
    if os.path.isfile(outputfile):
        return

    this_date_files = glob.glob(date_search)
    this_date_files.sort()

    # Run this in serial instead of using the parallel map functionality.
    # all_data = []
    # for f in this_date_files:
    #     all_data.append(load_file(f))

    if pool is None:
        all_data = []
        for f in this_date_files:
            all_data.append(load_file(f))
    else:
        results = pool.map_async(load_file, this_date_files)
        all_data = results.get()


    data_set = set_up_dataset(all_data[0])

    ids, ide, jds, jde, kds, kde = get_dims(all_data[0], section='d')
    for d in all_data:
        ims, ime, jms, jme, kms, kme = get_dims(d, section='m')
        its, ite, jts, jte, kts, kte = get_dims(d, section='t')

        if ims==ids:
            its = ids
        if ime==ide:
            ite = ide

        if jms==jds:
            jts = jds
        if jme==jde:
            jte = jde

        xts, xte = its - ims, ite - ims + 1
        yts, yte = jts - jms, jte - jms + 1
        zts, zte = kts - kms, kte - kms + 1

        xs, xe = its - ids, ite - ids + 1
        ys, ye = jts - jds, jte - jds + 1
        zs, ze = kts - kds, kte - kds + 1

        for v in d.variables:
            dims   = d[v].dims
            x_off, y_off = get_dim_offset(dims)

            if len(dims) == 2:
                data_set[v].values[ys:ye, xs:xe] = d[v].values[yts:yte, xts:xte]
            if len(dims) == 3:
                if dims[0] == "time":
                    data_set[v].values[:, ys:ye+y_off, xs:xe+x_off] = d[v].values[:, yts:yte+y_off, xts:xte+x_off]
                else:
                    data_set[v].values[zs:ze, ys:ye+y_off, xs:xe+x_off] = d[v].values[zts:zte, yts:yte+y_off, xts:xte+x_off]
            if len(dims) == 4:
                data_set[v].values[:,zs:ze, ys:ye+y_off, xs:xe+x_off] = d[v].values[:,zts:zte, yts:yte+y_off, xts:xte+x_off]

    print(outputfile)
    data_set.to_netcdf(outputfile)

def main(file_search = "icar_out_{ens}_*"):
    first_files = glob.glob(file_search.format(ens="000001"))
    first_files.sort()

    # For some reason running the parallelization this far out seems to have far worse performance...
    #  would map_async be faster for some reason?  I assume map is still parallel.
    # pool.map(agg_file, first_files)

    for f in first_files:
        agg_file(f)

# number of processors to parallelize reading the files over
# n_processors = 1
# pool = mp.Pool(n_processors)


def continuous(file_search):
    print("Running continuous aggregation, Ctrl-C to stop")
    while True:
        first_files = glob.glob(file_search.format(ens="000001"))
        first_files.sort()

        # skip the last file in the list as ICAR might still be running
        for f in first_files[:-1]:
            agg_file(f, verbose=False)

        time.sleep(10)


if __name__ == '__main__':
    if len(sys.argv) > 2:
        try:
            continuous(sys.argv[1])
        except KeyboardInterrupt:
            pass
    elif len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()
