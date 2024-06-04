#!/usr/bin/env python
import glob,os

import numpy as np

from bunch import Bunch
import mygis as io
import load_data

def adjust_p(p,h,dz):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals
    # h,dz in meters
    # slp = p/(1 - 2.25577E-5*h)**5.25588
    # p=slp*(1 - 2.25577E-5*(h+dz))**5.25588
    p*=(1 - 2.25577E-5*(h+dz))**5.25588

def update_base(base,filename,nz):
    data=load_data.cols(filename)
    nz=min(data.shape[0]-1,nz)
    base.z=data[:nz,0]
    base.dz=np.diff(data[:nz+1,0]).reshape((1,nz,1))
    base.th=data[:nz,1].reshape((1,nz,1))
    base.qv=data[:nz,2].reshape((1,nz,1))/1000.0

def main():
    filename="init"
    nx,nz,ny=(20.,10.,19)
    dims=[int(nx),int(nz),int(ny)]

    # po=np.log(100000.0)
    # p1=np.log(50000.0)
    # dp=(po-p1)/nz
    # p=np.exp(np.arange(po,p1,-dp)),
    lonmin=-110.0; lonmax=-100.0; dlon=(lonmax-lonmin)/nx
    latmin=35.0; latmax=45.0; dlat=(latmax-latmin)/ny

    base=Bunch(u=10.0,w=0.0,v=0.0,
               qv=0.0013,qc=0.0,
               p=100000.0,
               th=np.arange(273.0,300,(300-273.0)/nz),
               dz=400.0)
    base.z=np.arange(0,nz*base.dz,base.dz)
    if glob.glob("sounding.txt"):
        update_base(base,"sounding.txt",nz)
        nz=base.th.size
        dims=[nx,nz,ny]
    # base.p=base.p[:nz]

    u=np.zeros(dims,dtype="f")+base.u
    w=np.zeros(dims,dtype="f")+base.w
    v=np.zeros(dims,dtype="f")+base.v
    qv=np.zeros(dims,dtype="f")+base.qv
    # qv+=np.sin(np.arange(dims[1])/dims[1]*2*np.pi)+1
    qc=np.zeros(dims,dtype="f")+base.qc
    coscurve=np.cos(np.arange(dims[0])/dims[0]*2*np.pi+np.pi)+1
    # p-=coscurve.reshape((nx,1,1))*15000
    hgt=(coscurve*1000).reshape((int(nx),1)).repeat(ny,axis=1)

    lon=np.arange(lonmin,lonmax,dlon)
    lat=np.arange(latmin,latmax,dlat)
    lat,lon=np.meshgrid(lat,lon) #note that this appears "backwards" but these are C-style, fortran will be reversed
    dz=np.zeros(dims)+base.dz
    z=np.zeros(dims,dtype="f")+base.z.reshape((1,int(nz),1))+hgt.reshape((int(nx),1,int(ny)))

    layer1=(dz[:,0,:]/2)
    z[:,0,:]+=layer1
    for i in range(1,int(nz)):
        z[:,i,:]=z[:,i-1,:]+(dz[:,i-1,:]+dz[:,i,:])/2.0

    p=np.zeros(dims,dtype="f")+base.p# .reshape((1,nz,1))
    adjust_p(p,0.0,z)
    th=np.zeros(dims,dtype="f")+base.th.reshape((1,int(nz),1))


    d3dname=("x","z","y")
    d2dname=("x","y")
    othervars=[Bunch(data=lat,name="XLAT", dims=d2dname,dtype="f",attributes=dict(units="deg",  description="Latitude")),
               Bunch(data=lon,name="XLONG",dims=d2dname,dtype="f",attributes=dict(units="deg",  description="Longitude")),
               Bunch(data=v,  name="v",    dims=d3dname,dtype="f",attributes=dict(units="m/s",  description="Horizontal (y) wind speed")),
               Bunch(data=w,  name="w",    dims=d3dname,dtype="f",attributes=dict(units="m/s",  description="Vertical wind speed")),
               Bunch(data=qv, name="qv",   dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Water vapor mixing ratio")),
               Bunch(data=qc, name="qc",   dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Cloud water mixing ratio")),
               Bunch(data=p,  name="p",    dims=d3dname,dtype="f",attributes=dict(units="Pa",   description="Pressure")),
               Bunch(data=th, name="th",   dims=d3dname,dtype="f",attributes=dict(units="K",    description="Potential temperature")),
               Bunch(data=dz, name="dz",   dims=d3dname,dtype="f",attributes=dict(units="m",    description="Layer thickness")),
               Bunch(data=z,  name="z",    dims=d3dname,dtype="f",attributes=dict(units="m",    description="Layer Height AGL")),
               Bunch(data=hgt,name="hgt",  dims=d2dname,dtype="f",attributes=dict(units="m",    description="Terrain Elevation"))
               ]
    print(filename)
    fileexists=glob.glob(filename) or glob.glob(filename+".nc")
    if fileexists:
        print("Removing : "+fileexists[0])
        os.remove(fileexists[0])
    print(u.shape)
    print(lat.shape)
    print(hgt.shape)
    print(lon.shape)
    io.write(filename,  u,varname="u", dims=d3dname,dtype="f",attributes=dict(units="m/s",description="Horizontal (x) wind speed"),
            extravars=othervars)


if __name__ == '__main__':
    main()
