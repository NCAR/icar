#!/usr/bin/env python
import glob,os
from copy import copy

import numpy as np

from bunch import Bunch
import mygis as io
import load_data

case_study=2
wind_speed=10.0

xdomain=[800,500,100]
xdomain=[1500,1000,200]
xdomain=[1500,1000,500]
dx=2.0
# nx,nz,ny
master_dims=list([(xdomain[0]/dx*2,20,5),
                  (xdomain[1]/dx*2,20,5),
                  (xdomain[2]/dx*2,20,5)])

def adjust_p(p,h,dz):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals
    # h,dz in meters
    # slp = p/(1 - 2.25577E-5*h)**5.25588
    # p=slp*(1 - 2.25577E-5*(h+dz))**5.25588
    p*=(1 - 2.25577E-5*(h+dz))**5.25588

def update_base(base,filename,nz):
    """update the base information using data from a sounding file

    filename should be a space delimited text file with 3 columns
        height [m], potential temperature [K], and specific humidity [g/kg]"""
    print("Using Sounding from : "+filename)
    data=load_data.cols(filename)
    nz=min(data.shape[0]-1,nz)
    base.z=data[:nz,0]
    base.dz=np.diff(data[:nz+1,0]).reshape((1,nz,1,1))
    base.th=data[:nz,1].reshape((1,nz,1,1))
    base.qv=data[:nz,2].reshape((1,nz,1,1))/1000.0

def build_topography(experiment,dims):
    """create the topography to be used for a given case study"""

    # experiment                        D1      D2      D3
    # h         height of the hill    1800    1400    1040   [meters]
    # sigma     half-width              60      40       3.1 [grid cells]
    # z0        base of the hill      1700    2000    2200   [meters]
    # G         number of grids        420     250      52   [grid cells]

    Nx    = dims[3]                               # % length of domain  (grid cells)
    hm    = [1800.0,1400.0,1040.0][experiment]    # % mnt height (m)
    xm    = Nx/2.0                                # % mountain location in domain (grid cell)
    am    = [60.0,40.0,3.1][experiment]           # % mountain half-width (grid cells)

    dx    = 2000.0                                # % grid spacing (m)
    Lx    = Nx*dx                                 # % length of domain (m)
    x     = np.linspace(0,Lx,Nx)                  # % distance array (m)
    zo    = [1700.0,2000.0,2200.0][experiment]    # mountain base height (m) NOT REALLY USED CORRECTLY YET, NEED to truncate the sounding...
    zo    = 0.0


    zs=hm/(1.0+((x/dx-xm)/am)**2.)
    # zs-=zs[0]
    zs=zs.reshape((1,dims[3])).repeat(dims[2],axis=0)
    return zs

def main():
    filename="ideal_{}_{}".format(case_study,int(wind_speed))
    print(filename)
    nx,nz,ny=master_dims[case_study]
    dims=[1,int(nz),int(ny),int(nx)]

    # this is just arbitrary for now
    dlon=dx/111.1
    dlat=dx/111.1
    lonmin=-110.0; lonmax=lonmin+nx*dlon
    latmin=40.0; latmax=latmin+ny*dlat

    base=Bunch(u=wind_speed,w=0.0,v=0.0,
               qv=0.0013,qc=0.0,
               p=100000.0,
               th=np.arange(273.0,300,(300-273.0)/nz).reshape((nz,1,1)),
               dz=100.0)
    base.z=np.arange(0,nz*base.dz,base.dz)
    if glob.glob("sounding.txt"):
        update_base(base,"sounding.txt",nz)
        nz=base.th.size
        dims=[1,nz,ny,nx]

    udims=copy(dims)
    udims[-1]+=1
    vdims=copy(dims)
    vdims[-2]+=1
    u=np.zeros(udims,dtype="f")+base.u
    w=np.zeros(dims,dtype="f")+base.w
    v=np.zeros(vdims,dtype="f")+base.v
    qv=np.zeros(dims,dtype="f")+base.qv
    qc=np.zeros(dims,dtype="f")+base.qc

    # simple topography = a cosine
    # coscurve=np.cos(np.arange(dims[3])/dims[3]*2*np.pi+np.pi)+1
    # hgt=(coscurve*1000).reshape((1,nx)).repeat(ny,axis=0)
    hgt=build_topography(case_study,dims)

    lon=np.arange(lonmin,lonmax,dlon)[:int(nx)]
    lat=np.arange(latmin,latmax,dlat)[:ny]
    lon,lat=np.meshgrid(lon,lat)

    ulon=np.arange(lonmin-dlon/2,lonmax+dlon/2,dlon)[:int(nx)+1]
    ulat=np.arange(latmin,latmax,dlat)[:ny]
    ulon,ulat=np.meshgrid(ulon,ulat)

    vlon=np.arange(lonmin,lonmax,dlon)[:int(nx)]
    vlat=np.arange(latmin-dlat/2,latmax+dlat/2,dlat)[:ny+1]
    vlon,vlat=np.meshgrid(vlon,vlat)

    dz=np.zeros(dims)+base.dz
    z=np.zeros(dims,dtype="f")+base.z.reshape((1,nz,1,1))+hgt.reshape((1,1,ny,int(nx)))

    layer1=(dz[0,:,:]/2)
    z[0,:,:]+=layer1
    for i in range(1,int(nz)):
        z[:,i,:,:]=z[:,i-1,:,:]+(dz[:,i-1,:,:]+dz[:,i,:,:])/2.0

    p=np.zeros(dims,dtype="f")+base.p
    adjust_p(p,0.0,z)
    th=np.zeros(dims,dtype="f")+base.th

    lat=lat.reshape((1,ny,int(nx)))
    lon=lon.reshape((1,ny,int(nx)))
    hgt=hgt.reshape((1,ny,int(nx)))

    d3dname=("t","z","y","x")
    ud3dname=("t","z","y","xu")
    ud2dname=("t","y","xu")
    vd3dname=("t","z","yv","x")
    vd2dname=("t","yv","x")
    d2dname=("t","y","x")
    g=9.81
    othervars=[Bunch(data=v,  name="V",    dims=vd3dname,dtype="f",attributes=dict(units="m/s",  description="Horizontal (y) wind speed")),
               Bunch(data=w,  name="W",     dims=d3dname,dtype="f",attributes=dict(units="m/s",  description="Vertical wind speed")),
               Bunch(data=qv, name="QVAPOR",dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Water vapor mixing ratio")),
               Bunch(data=qc, name="QCLOUD",dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Cloud water mixing ratio")),
               Bunch(data=qc, name="QICE",  dims=d3dname,dtype="f",attributes=dict(units="kg/kg",description="Cloud ice mixing ratio")),
               Bunch(data=p*0,name="P",     dims=d3dname,dtype="f",attributes=dict(units="Pa",   description="Pressure (perturbation)")),
               Bunch(data=p,  name="PB",    dims=d3dname,dtype="f",attributes=dict(units="Pa",   description="Pressure (base)")),
               Bunch(data=th-300, name="T", dims=d3dname,dtype="f",attributes=dict(units="K",    description="Potential temperature")),
               Bunch(data=dz, name="dz",    dims=d3dname,dtype="f",attributes=dict(units="m",    description="Layer thickness")),
               Bunch(data=z,  name="Z",     dims=d3dname,dtype="f",attributes=dict(units="m",    description="Layer Height AGL (also ASL here)")),
               Bunch(data=z*0,name="PH",    dims=d3dname,dtype="f",attributes=dict(units="m2/s2",description="Geopotential Height ASL (perturbation)")),
               Bunch(data=z*g,name="PHB",   dims=d3dname,dtype="f",attributes=dict(units="m2/s2",description="Geopotential Height ASL (base)")),
               Bunch(data=lat,name="XLAT",  dims=d2dname,dtype="f",attributes=dict(units="deg",  description="Latitude")),
               Bunch(data=lon,name="XLONG", dims=d2dname,dtype="f",attributes=dict(units="deg",  description="Longitude")),
               Bunch(data=ulat,name="XLAT_U",dims=ud2dname,dtype="f",attributes=dict(units="deg",  description="Latitude on U stagger")),
               Bunch(data=ulon,name="XLONG_U",dims=ud2dname,dtype="f",attributes=dict(units="deg",  description="Longitude on U stagger")),
               Bunch(data=vlat,name="XLAT_V",dims=vd2dname,dtype="f",attributes=dict(units="deg",  description="Latitude on V stagger")),
               Bunch(data=vlon,name="XLONG_V",dims=vd2dname,dtype="f",attributes=dict(units="deg",  description="Longitude on V stagger")),
               Bunch(data=hgt,name="HGT",   dims=d2dname,dtype="f",attributes=dict(units="m",    description="Terrain Elevation"))
               ]
    fileexists=glob.glob(filename) or glob.glob(filename+".nc")
    if fileexists:
        print("Removing : "+fileexists[0])
        os.remove(fileexists[0])

    io.write(filename,  u,varname="U", dims=ud3dname,dtype="f",attributes=dict(units="m/s",description="Horizontal (x) wind speed"),
            extravars=othervars)


if __name__ == '__main__':
    for case in range(3):
        for ws in [5,10,15,25]:
            wind_speed=float(ws)
            case_study=case
            main()
