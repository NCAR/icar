#!/usr/bin/env python
import glob,os,sys
from copy import copy

import numpy as np

from bunch import Bunch
import mygis

g=9.81

def adjust_p(p,hin,hout):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals
    # h,dz in meters
    slp = p/(1 - 2.25577E-5*hin)**5.25588
    p[:]=slp*(1 - 2.25577E-5*hout)**5.25588


def main(inputfile,sounding_file=None):
    filename="icar_"+inputfile
    print(filename)
    yaxis=2
    yaxis2d=1
    u   =mygis.read_nc(inputfile,"U").data.repeat(2,axis=yaxis)
    # w   =mygis.read_nc(inputfile,"W").data.repeat(2,axis=yaxis)
    v   =mygis.read_nc(inputfile,"V").data.repeat(2,axis=yaxis)[:,:,:-1,:] # v has one extra cell in y, so when doubling ydim we have to remove gridcell
    qv  =mygis.read_nc(inputfile,"QVAPOR").data.repeat(2,axis=yaxis)
    qc  =mygis.read_nc(inputfile,"QCLOUD").data.repeat(2,axis=yaxis)
    qi  =mygis.read_nc(inputfile,"QICE").data.repeat(2,axis=yaxis)
    th  =mygis.read_nc(inputfile,"T").data.repeat(2,axis=yaxis)
    pb  =mygis.read_nc(inputfile,"PB").data.repeat(2,axis=yaxis)
    p   =mygis.read_nc(inputfile,"P").data.repeat(2,axis=yaxis) + pb
    phb =mygis.read_nc(inputfile,"PHB").data.repeat(2,axis=yaxis)
    ph  =mygis.read_nc(inputfile,"PH").data.repeat(2,axis=yaxis) + phb

    hgt =mygis.read_nc(inputfile,"HGT").data.repeat(2,axis=yaxis2d)
    land=mygis.read_nc(inputfile,"XLAND").data.repeat(2,axis=yaxis2d)

    nt,nz,ny,nx=qv.shape
    print(nx,ny,nz)
    dims=np.array(qv.shape)

    z=(ph)/g
    dz=np.diff(z,axis=1)
    # dz shape = (time,nz-1,ny,nx)
    # ph/phb are defined between model levels, we want z in the middle of each model level
    # e.g. z[0]=0, but wrfz[0]=dz[0]/2 where dz[0]=z[1]-z[0]
    wrfz=(z[:,:-1,:,:]+z[:,1:,:,:])/2
    # p=units.z2p(wrfz,100000)
    # wrfz=np.zeros(dz.shape)
    # wrfz[:,0,...]=dz[:,0,...]/2+hgt
    # for i in range(1,nz):
        # wrfz[:,i,:,:]=(dz[:,i,:,:]+dz[:,i-1,:,:])/2+wrfz[:,i-1,:,:]

    mean_dz=dz[0,...].mean(axis=1).mean(axis=1)
    print("MEAN LEVELS:")
    print("dz_levels=[")
    for i in range(0,nz,10):
        curlist=[str(cur) for cur in mean_dz[i:i+10]]
        print(",".join(curlist)+",")

    if i+10<mean_dz.shape[0]:
        curlist=[str(cur) for cur in mean_dz[i+10:]]
        print(",".join(curlist)+"]")
    else:
        print("]")


    print("FIRST LEVELS:")
    print("dz_levels=[")
    for i in range(0,nz,10):
        curlist=[str(cur) for cur in dz[0,i:i+10,0,0]]
        print(",".join(curlist)+",")

    if i+10<dz.shape[1]:
        curlist=[str(cur) for cur in dz[0,i+10:,0,0]]
        print(",".join(curlist)+"]")
    else:
        print("]")


    # dz=np.zeros(dz.shape)+mean_dz[np.newaxis,:,np.newaxis,np.newaxis]
    dz=np.zeros(dz.shape)+dz[:,:,:,np.newaxis,0]
    z=np.zeros(dz.shape)
    z[:,0,...]=dz[:,0,...]/2+hgt
    for i in range(1,nz):
        z[:,i,:,:]=(dz[:,i,:,:]+dz[:,i-1,:,:])/2+z[:,i-1,:,:]

    # adjust_p(p,wrfz,z)

    dx=mygis.read_attr(inputfile,"DX")
    if type(dx)==np.ndarray:
        dx=dx[0]
    dlon=dx/111.1
    dlat=dx/111.1
    lonmin=-110.0; lonmax=lonmin+nx*dlon
    latmin=40.0; latmax=latmin+ny*dlat

    udims=copy(dims)
    udims[-1]+=1
    vdims=copy(dims)
    vdims[-2]+=1

    lon=np.arange(lonmin,lonmax,dlon)[:nx]
    lat=np.arange(latmin,latmax,dlat)[:ny]
    lon,lat=np.meshgrid(lon,lat)

    ulon=np.arange(lonmin-dlon/2,lonmax+dlon/2,dlon)[:nx+1]
    ulat=np.arange(latmin,latmax,dlat)[:ny]
    ulon,ulat=np.meshgrid(ulon,ulat)

    vlon=np.arange(lonmin,lonmax,dlon)[:nx]
    vlat=np.arange(latmin-dlat/2,latmax+dlat/2,dlat)[:ny+1]
    vlon,vlat=np.meshgrid(vlon,vlat)

    lat=lat.reshape((1,ny,nx))
    lon=lon.reshape((1,ny,nx))
    hgt=hgt.reshape((1,ny,nx))

    d3dname=("t","z","y","x")
    ud3dname=("t","z","y","xu")
    ud2dname=("t","y","xu")
    vd3dname=("t","z","yv","x")
    vd2dname=("t","yv","x")
    d2dname=("t","y","x")

    othervars=[Bunch(data=v,   name="V",      dims=vd3dname,dtype="f",attributes=dict(units="m/s",  description="Horizontal (y) wind speed")),
               # Bunch(data=w,   name="W",      dims=d3dname, dtype="f",attributes=dict(units="m/s",  description="Vertical wind speed")),
               Bunch(data=qv,  name="QVAPOR", dims=d3dname, dtype="f",attributes=dict(units="kg/kg",description="Water vapor mixing ratio")),
               Bunch(data=qc,  name="QCLOUD", dims=d3dname, dtype="f",attributes=dict(units="kg/kg",description="Cloud water mixing ratio")),
               Bunch(data=qi,  name="QICE",   dims=d3dname, dtype="f",attributes=dict(units="kg/kg",description="Cloud ice mixing ratio")),
               Bunch(data=p*0, name="P",      dims=d3dname, dtype="f",attributes=dict(units="Pa",   description="Pressure (perturbation)")),
               Bunch(data=p,   name="PB",     dims=d3dname, dtype="f",attributes=dict(units="Pa",   description="Pressure (base)")),
               Bunch(data=th,  name="T",      dims=d3dname, dtype="f",attributes=dict(units="K",    description="Potential temperature")),
               Bunch(data=dz,  name="dz",     dims=d3dname, dtype="f",attributes=dict(units="m",    description="Layer thickness")),
               Bunch(data=z,   name="Z",      dims=d3dname, dtype="f",attributes=dict(units="m",    description="Layer Height AGL (also ASL here)")),
               Bunch(data=wrfz,name="WRFZ",   dims=d3dname, dtype="f",attributes=dict(units="m",    description="Layer Height AGL (also ASL here) in the WRF input")),
               # Bunch(data=z*0, name="PH",     dims=d3dname, dtype="f",attributes=dict(units="m2/s2",description="Geopotential Height ASL (perturbation)")),
               # Bunch(data=z*g, name="PHB",    dims=d3dname, dtype="f",attributes=dict(units="m2/s2",description="Geopotential Height ASL (base)")),
               Bunch(data=lat, name="XLAT",   dims=d2dname, dtype="f",attributes=dict(units="deg",  description="Latitude")),
               Bunch(data=lon, name="XLONG",  dims=d2dname, dtype="f",attributes=dict(units="deg",  description="Longitude")),
               Bunch(data=ulat,name="XLAT_U", dims=ud2dname,dtype="f",attributes=dict(units="deg",  description="Latitude on U stagger")),
               Bunch(data=ulon,name="XLONG_U",dims=ud2dname,dtype="f",attributes=dict(units="deg",  description="Longitude on U stagger")),
               Bunch(data=vlat,name="XLAT_V", dims=vd2dname,dtype="f",attributes=dict(units="deg",  description="Latitude on V stagger")),
               Bunch(data=vlon,name="XLONG_V",dims=vd2dname,dtype="f",attributes=dict(units="deg",  description="Longitude on V stagger")),
               Bunch(data=hgt, name="HGT",    dims=d2dname, dtype="f",attributes=dict(units="m",    description="Terrain Elevation")),
               Bunch(data=land,name="XLAND",  dims=d2dname, dtype="f",attributes=dict(units="",     description="Land Mask [1=land,2=water]"))
               ]
    fileexists=glob.glob(filename) or glob.glob(filename+".nc")
    if fileexists:
        print("Removing : "+fileexists[0])
        os.remove(fileexists[0])

    mygis.write(filename,  u,varname="U", dims=ud3dname,dtype="f",attributes=dict(units="m/s",description="Horizontal (x) wind speed"),
            extravars=othervars)


if __name__ == '__main__':
    if len(sys.argv)>1:
        filename=sys.argv[1]
    else:
        filename="wrfinput_d01"
    if len(sys.argv)>2:
        sounding_file=sys.argv[2]
    else:
        sounding_file=None

    main(filename,sounding_file)
