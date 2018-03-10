#!/usr/bin/env python

"""
SYNOPSIS

    template_argparse.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION

    
"""
from __future__ import print_function
import sys
import os
import traceback
import argparse
import glob

import numpy as np
from mpl_toolkits.basemap import Basemap
# import matplotlib.pyplot as plt
from netCDF4 import Dataset
import mygis

from bunch import Bunch

# note this is only the default that can be overwritten by CLI args
global topo_dataset
topo_dataset="/glade/u/home/gutmann/work/topo_1km_global.nc"

# convert km to meters
km = 1000.0

def create_mapset(center_lat, center_lon, width, height):
    """create a basemap instance"""
    width_in_meters=width*km
    height_in_meters = height*km

    lat1 = center_lat+height/166.
    lat2 = center_lat-height/166.

    # setup Lambert Conformal basemap.
    print("Setting up map projection")
    m = Basemap(width=width_in_meters,height=height_in_meters,projection='lcc',
                resolution='l',lat_1=lat1,lat_2=lat2,lat_0=center_lat,lon_0=center_lon)
    return m


def load_topo_data(latmin,latmax,lonmin,lonmax):
    data = Dataset(topo_dataset)
    lons = data.variables["X"][:]
    lats = data.variables["Y"][:]
    
    left  = np.where(lons<lonmin)[0][-1]
    right = np.where(lons>lonmax)[0][0]

    # note lats are reversed from what you might expect
    # e.g. 90N would be at position 0, not -1
    bottom = np.where(lats<latmin)[0][0]
    top    = np.where(lats>latmax)[0][-1]
    
    print("Reading topo data")
    topoin = data.variables["topo"][bottom:top-1:-1,left:right+1]
    
    lat=lats[bottom:top-1:-1]
    lon=lons[left:right+1]
    
    lon,lat=np.meshgrid(lon,lat)
    
    if type(topoin)==np.ma.core.MaskedArray:
        mask=topoin.mask.copy()
        topoin.mask[:]=False
        topoin[mask]=0.0
        
    return Bunch(data=topoin,lon=lon,lat=lat)


def setup_grid(lat,lon,width,height,dx,m):
    """setup the latitude and longitude grids over the domain"""
    x,y = m(lon,lat)
    xll = x-width/2
    yll = y-height/2
    
    allx = np.arange(xll, xll+width+dx, dx)
    ally = np.arange(yll, yll+height+dx, dx)
    nx=len(allx)
    ny=len(ally)
    latgrid=np.zeros((ny,nx))
    longrid=np.zeros((ny,nx))

    ulat=np.zeros((ny,nx+1))
    ulon=np.zeros((ny,nx+1))
    vlat=np.zeros((ny+1,nx))
    vlon=np.zeros((ny+1,nx))
    
    for i,y in enumerate(ally):
        for j,x in enumerate(allx):
           longrid[i,j],latgrid[i,j] = m(x,y,inverse=True)
           ulon[i,j],ulat[i,j] = m(x-dx/2,y,inverse=True)
           vlon[i,j],vlat[i,j] = m(x,y-dx/2,inverse=True)
           if i==0:
               vlon[ny,j],vlat[ny,j] = m(x,ally[-1]+dx/2,inverse=True)
           
        ulon[i,nx],ulat[i,nx] = m(x+dx/2,y,inverse=True)
        
    
    return latgrid,longrid,ulat,ulon,vlat,vlon
    

def topo2grid(topo,lat,lon):
    """convert a topodata set to match a set of lat,lon grids"""
    
    output_data=np.zeros(lat.shape)
    output_n=np.zeros(lat.shape)
    
    max_dlat=np.max(lat[1:,:]-lat[:-1,:])/2.0
    max_dlon=np.max(lon[:,1:]-lon[:,:-1])/2.0

    hi_dlat=np.max(topo.lat[1:,:]-topo.lat[:-1,:])/2.0
    hi_dlon=np.max(topo.lon[:,1:]-topo.lon[:,:-1])/2.0
    ratio=max(hi_dlat/max_dlat, hi_dlon/max_dlon)
    window=max(2,2*int(np.floor(ratio)))
    
    ny,nx=topo.data.shape
    nyl,nxl=lat.shape

    print("Converting topographic grid")
    if (nyl>100):
        print("for large grids, this can take a while")
    for i in range(ny):
        if (nyl>100):
            if (i%20)==0:
                print("\r{0:4.0f}%".format(100.0*i/float(ny)), end="")
                sys.stdout.flush()
        dists=(lat-topo.lat[i,0])**2 + (lon-topo.lon[i,0])**2
        cury,curx=np.unravel_index(np.argmin(dists),dists.shape)
        for j in range(nx):
            xmin=max(0,curx-window)
            xmax=min(nxl,curx+window)
            ymin=max(0,cury-window)
            ymax=max(nyl,cury+window)
            
            dists=(lat[ymin:ymax,xmin:xmax]-topo.lat[i,j])**2 + (lon[ymin:ymax,xmin:xmax]-topo.lon[i,j])**2
            cury,curx=np.unravel_index(np.argmin(dists),dists.shape)
            cury+=ymin
            curx+=xmin
            if ((np.abs(lat[cury,curx]-topo.lat[i,j]) < max_dlat) 
                and (np.abs(lon[cury,curx]-topo.lon[i,j]) < max_dlon)) :
                
                output_data[cury,curx]+=topo.data[i,j]
                output_n[cury,curx]+=1
            
    
    output_n[output_n==0]=1
    
    print("  Finished")
    return Bunch(topo=output_data/output_n, lat=lat, lon=lon)
    

def write_outputfile(filename,dataset,mapset,dx):
    """docstring for write_outputfile"""
    globalatts=mapset.projparams
    globalatts.pop("x_0")
    globalatts.pop("y_0")
    globalatts.pop("units")
    globalatts["dx"]=dx
    
    latvar=Bunch(data=dataset.lat,name="lat",dims=('lat','lon'),dtype='f',attributes=Bunch(long_name="latitude",units="degrees"))
    lonvar=Bunch(data=dataset.lon,name="lon",dims=('lat','lon'),dtype='f',attributes=Bunch(long_name="longitude",units="degrees"))
    ulatvar=Bunch(data=dataset.ulat,name="lat_u",dims=('lat','lon_u'),dtype='f',attributes=Bunch(long_name="latitude_ugrid",units="degrees"))
    ulonvar=Bunch(data=dataset.ulon,name="lon_u",dims=('lat','lon_u'),dtype='f',attributes=Bunch(long_name="longitude_ugrid",units="degrees"))
    vlatvar=Bunch(data=dataset.vlat,name="lat_v",dims=('lat_v','lon'),dtype='f',attributes=Bunch(long_name="latitude_vgrid",units="degrees"))
    vlonvar=Bunch(data=dataset.vlon,name="lon_v",dims=('lat_v','lon'),dtype='f',attributes=Bunch(long_name="longitude_vgrid",units="degrees"))
    landvar=Bunch(data=dataset.xland,name="xland",dims=('lat','lon'),dtype='i',attributes=Bunch(long_name="land-sea-mask",units="[0,1]"))
    evars=[latvar,lonvar,ulatvar,ulonvar,vlatvar,vlonvar,landvar]
    # evars=[latvar,lonvar]
    
    print("Writing: "+filename)
    mygis.write(filename,
                dataset.topo,varname="HGT",dims=('lat','lon'),dtype='f', attributes=Bunch(long_name="topography",units="m"),
                extravars=evars, global_attributes=globalatts, history="data from make_domain.py")
    

def main (lat, lon, width, height, dx, outputfile):

    m = create_mapset(lat, lon, width, height)
    
    latgrid,longrid,ulat,ulon,vlat,vlon = setup_grid(lat, lon, width*km, height*km, dx, m)
    
    latmin,latmax,lonmin,lonmax = latgrid.min(),latgrid.max(), longrid.min(), longrid.max()
    topo = load_topo_data(latmin,latmax,lonmin,lonmax)
    
    topo=topo2grid(topo,latgrid,longrid)
    topo.ulat=ulat
    topo.ulon=ulon
    topo.vlat=vlat
    topo.vlon=vlon
    topo.xland=np.ones(latgrid.shape)
    
    write_outputfile(outputfile, topo, m, dx)

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Create an ICAR high-res domain')
        parser.add_argument('center_lat',action='store', type=float,
                            help="Center latitude of domain")
        parser.add_argument('center_lon',action='store', type=float,
                            help="Center longitude of domain")
        parser.add_argument('width',action='store', type=float,
                            help="Width of domain (km)")
        parser.add_argument('height',action='store', type=float,
                            help="Height of domain (km)")
        parser.add_argument('-d','--dx',dest="dx",type=float, default=4000,
                            help="Grid size [m] (default=4000)")
        parser.add_argument('-o','--outputfile',dest="outputfile",default="domain.nc",
                            help="Output file name (default=domain.nc)")
        parser.add_argument('-z','--dem',dest="dem",default=topo_dataset,
                            help="Name of base topography file")
        parser.add_argument('-v', '--version',action='version',
                version='make_domain 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        topo_dataset = args.dem
        
        exit_code = main(args.center_lat,args.center_lon,args.width,args.height,args.dx,args.outputfile)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
