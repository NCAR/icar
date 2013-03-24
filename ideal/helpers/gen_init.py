#!/usr/bin/env python
import glob,os

import numpy as np

from bunch import Bunch
from stat_down import myio as io
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
    nx,nz,ny=(100.,20.,3)
    dims=[nx,nz,ny]
    
    # po=np.log(100000.0)
    # p1=np.log(50000.0)
    # dp=(po-p1)/nz
    # p=np.exp(np.arange(po,p1,-dp)),
    base=Bunch(u=10.0,w=0.0,v=0.0,
               qv=0.0013,qc=0.0,
               p=100000.0,
               th=np.arange(273.0,300,(300-273.0)/nz),
               dz=200.0)
    
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
    hgt=(coscurve*1000).reshape((nx,1)).repeat(ny,axis=1)
    dz=np.zeros(dims)+base.dz
    z=np.zeros(dims,dtype="f")+base.z.reshape((1,nz,1))+hgt.reshape((nx,1,ny))
    
    layer1=(dz[:,0,:]/2)
    print(layer1.shape)
    z[:,0,:]+=layer1
    for i in range(1,int(nz)):
        z[:,i,:]=z[:,i-1,:]+(dz[:,i-1,:]+dz[:,i,:])/2.0
    
    p=np.zeros(dims,dtype="f")+base.p# .reshape((1,nz,1))
    adjust_p(p,0.0,z)
    print(p[0,-1,0],p.min())
    th=np.zeros(dims,dtype="f")+base.th.reshape((1,nz,1))
    
    
    d3dname=("x","z","y")
    d2dname=("x","y")
    othervars=[Bunch(data=v,  name="v",  dims=d3dname,dtype="f",attributes=dict(units="m/s",  description="Horizontal (y) wind speed")),
               Bunch(data=w,  name="w",  dims=d3dname,dtype="f",attributes=dict(units="m/s",  description="Vertical wind speed")),
               Bunch(data=qv, name="qv", dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Water vapor mixing ratio")),
               Bunch(data=qc, name="qc", dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Cloud water mixing ratio")),
               Bunch(data=p,  name="p",  dims=d3dname,dtype="f",attributes=dict(units="Pa",   description="Pressure")),
               Bunch(data=th, name="th", dims=d3dname,dtype="f",attributes=dict(units="K",    description="Potential temperature")),
               Bunch(data=dz, name="dz", dims=d3dname,dtype="f",attributes=dict(units="m",    description="Layer thickness")),
               Bunch(data=z,  name="z",  dims=d3dname,dtype="f",attributes=dict(units="m",    description="Layer Height AGL")),
               Bunch(data=hgt,name="hgt",dims=d2dname,dtype="f",attributes=dict(units="m",    description="Terrain Elevation"))
               ]
    if glob.glob(filename):
        os.remove(filename)
    io.write(filename,  u,varname="u", dims=d3dname,dtype="f",attributes=dict(units="m/s",description="Horizontal (x) wind speed"),
            extravars=othervars)


if __name__ == '__main__':
    main()

