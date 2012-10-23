#!/usr/bin/env python

"""
SYNOPSIS

    swim_narr.py [-h] [--verbose] [-v, --version] [NARR_NC_file_search_pattern] [topofile]

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    swim_narr.py

EXIT STATUS

    None

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    0.2 - the hey it sort of works version!
    
"""
# Python standard modules
import sys
import os
import traceback
#import argparse
from multiprocessing import Pool
import glob
import time
# Standard 3rd party modules
import numpy.fft as fft
from scipy.signal import convolve
import numpy as np
# my common modules
# import netCDF4
# import mygis
# from nc import NC_writer
import units
# import lop_model_2d as lop
from bunch import Bunch
# swim_narr specific modules
import swim_io # this class combines the functionality of the three file IO classes below
import swim as swim_lib
import swim_fast as swim
import fast_mean

NPROCESSORS=2
R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

# this is the half kernel window size for various smoothin applied to the input data
halfk=8
# Note, this defines the subdomain within the WRF 2km model run to use
# sub_y0=125-halfk
# sub_y1=400+halfk
# sub_x0=220-halfk
# sub_x1=400+halfk
# io_ratio=8.0 #this is the ratio of the input grid (NARR 32km) to the output grid (WRF 2km)
# wrfres=2000.0
# to use the entire WRF domain use this instead, but be careful of the smoothing above...: 
sub_y0=0
sub_y1=None
sub_x0=0
sub_x1=None
io_ratio=8.0 #this is the ratio of the input grid (NARR 32km) to the output grid (WRF 4km)
wrfres=4000.0

def match_xy(lat1,lon1,lat2,lon2):
    N=lat2.shape
    x=np.zeros(N)
    y=np.zeros(N)
    winhalfsize=2

    for i in range(N[1]):
        j=0
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        latwin=lat1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
        lonwin=lon1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
        for j in range(N[0]):
            if (prevx!=lastx) or (prevy!=lasty):
                latwin=lat1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
                lonwin=lon1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
                (prevx,prevy)=(lastx,lasty)
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            
            # print(j,i,newx,newy)
            lastx=(newx-winhalfsize)+lastx
            lasty=(newy-winhalfsize)+lasty
            x[j,i]=lastx
            y[j,i]=lasty
        
    return (x,y)
 
def bilin_weights(yi,y,xi,x):
    x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1-x0 # np.abs((xi-x[1])/(x[1]-x[0]))
    x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1-x2 # np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    f1=(yi-y5)/(y6-y5)
    f2=1-f1# (y6-yi)/(y6-y5)
    return np.array([x1*f2,x0*f2,x3*f1,x2*f1])
    
    
def match_xy_bilin(lat1,lon1,lat2,lon2):
    N=lat2.shape
    out=np.zeros((N[0],N[1],4,3))
    x=np.zeros(4).astype('i')
    y=np.zeros(4).astype('i')
    winhalfsize=5
    
    dxinc=np.sign(lon1[1,1]-lon1[0,0]).astype('i')
    dyinc=np.sign(lat1[1,1]-lat1[0,0]).astype('i')

    for i in range(N[1]):
        j=0
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        latwin=lat1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
        lonwin=lon1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
        for j in range(N[0]):
            if (prevx!=lastx) or (prevy!=lasty):
                latwin=lat1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
                lonwin=lon1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
                (prevx,prevy)=(lastx,lasty)
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            lastx=(newx-winhalfsize)+lastx
            lasty=(newy-winhalfsize)+lasty
            x[0]=newx
            y[0]=newy
            if latwin[newy,newx]<lat2[j,i]:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=newx+dxinc
                    x[2]=newx
                    x[3]=newx+dxinc
                    y[1]=newy
                    y[2]=newy+dyinc
                    y[3]=newy+dyinc
                else:
                    x[1]=newx-dxinc
                    x[2]=newx
                    x[3]=newx-dxinc
                    y[1]=newy
                    y[2]=newy+dyinc
                    y[3]=newy+dyinc
            else:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=newx+dxinc
                    x[2]=newx
                    x[3]=newx+dxinc
                    y[1]=newy
                    y[2]=newy-dyinc
                    y[3]=newy-dyinc
                else:
                    x[1]=newx-dxinc
                    x[2]=newx
                    x[3]=newx-dxinc
                    y[1]=newy
                    y[2]=newy-dyinc
                    y[3]=newy-dyinc
                    
            # bilinear interpolation for an arbitrary grid spacing 
            # (must be a grid, but this will handle a lot of irregularities)
            weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])

            out[j,i,:,0]=(y-winhalfsize)+prevy
            out[j,i,:,1]=(x-winhalfsize)+prevx
            out[j,i,:,2]=weights
        
    return out
              
                
def convert_p(p,h,dz):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals
    # h,dz in meters
    slp = p/(1 - 2.25577E-5*h)**5.25588
    pout=slp*(1 - 2.25577E-5*(h+dz))**5.25588
    return pout

class NARR_Reader(object):
    _filenames=None
    _sfc_files=None
    _curfile=0
    x=-99
    y=-99
    _bilin=False
    _nn=True
    _geoLUT=None
    topo=None
    # atm variables in NARR 3d files
    tvar='TMP_221_ISBL'
    hvar='HGT_221_ISBL'
    qcvar='CLWMR_221_ISBL'
    qvvar='SPF_H_221_ISBL'
    qivar='ICMR_221_ISBL'
    wvar='V_VEL_221_ISBL'
    uvar='U_GRD_221_ISBL'
    vvar='V_GRD_221_ISBL'
    # land surface variables in NARR SFC files
    lhvar="LHTFL_221_SFC_ave3h"
    shvar="SHTFL_221_SFC_ave3h"
    LWuvar="ULWRF_221_SFC_ave3h"
    albvar="ALBDO_221_SFC"
    
    def init_xy(self):
        narrfilename='/glade/scratch/gutmann/usbr/narr_data/baseline_info/HGT_2006010100.nc'
        if wrfres==2000:
            wrffilename='/glade/scratch/gutmann/usbr/narr_data/baseline_info/2km_wrf_input_d01'
        elif wrfres==4000:
            wrffilename='/glade/scratch/gutmann/usbr/narr_data/baseline_info/4km_wrf_output.nc'
        nlat=swim_io.read_nc(narrfilename,var='latitude').data
        nlon=swim_io.read_nc(narrfilename,var='longitude').data
        wlat=swim_io.read_nc(wrffilename,var='XLAT').data[0,sub_y0:sub_y1,sub_x0:sub_x1]
        wlon=swim_io.read_nc(wrffilename,var='XLONG').data[0,sub_y0:sub_y1,sub_x0:sub_x1]
        hgt=swim_io.read_nc(narrfilename,var='HGT_1-HYBL').data[:,:]
        
        if self._nn:
            (x,y)=match_xy(nlat,nlon,wlat,wlon)
            self.x=x.astype('i')
            self.y=y.astype('i')
            self.topo=hgt[y,x]
        elif self._bilin:
            self._geoLUT=match_xy_bilin(nlat,nlon,wlat,wlon)
            curx=(self._geoLUT[:,:,:,1]).astype('i')
            cury=(self._geoLUT[:,:,:,0]).astype('i')
            w=self._geoLUT[:,:,:,2]
            self.topo=np.zeros(curx.shape[0:2],dtype=np.float32)
            # for i in range(4):self.topo+=hgt[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            for i in range(4):self.topo+=hgt[cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            self.topo=self.topo[halfk:-halfk,halfk:-halfk]
            
        
    def __init__(self, file_search,sfc=None,nn=True, bilin=False, *args, **kwargs):
        super(NARR_Reader,self).__init__(*args, **kwargs)
        self._filenames=np.sort(glob.glob(file_search))
        print('Calculating XY match lookup table, this may take a while.')
        if sfc!=None:
            self._sfc_files=glob.glob(sfc)
        if nn:
            self._nn=True
            self._bilin=False
        elif bilin:
            self._nn=False
            self._bilin=True
        self.init_xy()
    
    
    # we are our own iterator...
    def __iter__(self):
        return self

    def nextNN(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        minx=self.x.min()
        maxx=self.x.max()+1
        miny=self.y.min()
        maxy=self.y.max()+1
        curx=self.x-minx
        cury=self.y-miny
        
        d=swim_io.Dataset(self._filenames[curfile], 'r')
        temperature=d.variables[self.tvar][:,miny:maxy,minx:maxx][:,curx,cury]
        pressure=d.variables[self.hvar][:,miny:maxy,minx:maxx][:,curx,cury]
        specific_humidity=d.variables[self.qcvar][:,miny:maxy,minx:maxx][:,curx,cury]
        relative_humidity=d.variables[self.qvvar][:,miny:maxy,minx:maxx][:,curx,cury]
        wind_u=d.variables[self.uvar][:,miny:maxy,minx:maxx][:,curx,cury]
        wind_v=d.variables[self.vvar][:,miny:maxy,minx:maxx][:,curx,cury]
        d.close()

        datestr=self._filenames[curfile].split('.')[1]

        self._curfile+=1
        
        return Bunch(ta=temperature, p=pressure,
                     sh=specific_humidity, 
                     rh=relative_humidity, 
                     u=wind_u, v=wind_v,
                     date=datestr)
    
    def convolve2to3(self,data,kern,mode="same"):
        # n=data.shape
        # print(kern.shape)
        if mode!="same":
            i=0
            datatmp=convolve(data[i,:,:],kern,mode=mode)
            n=datatmp.shape
            outputdata=np.zeros((data.shape[0],n[0],n[1]))
            outputdata[0,:,:]=datatmp
            for i in range(data.shape[0]-1):
                outputdata[i+1,:,:]=convolve(data[i+1,:,:],kern,mode=mode)
        else:
            for i in range(data.shape[0]):
                data[i,:,:]=convolve(data[i,:,:],kern,mode=mode)
        return data
    
    # def smooth_wind(self,wind):
    #     kernelsize=20
    #     gaussian=gauss_kern(kernelsize)
    #     windsm=self.convolve2to3(wind,gaussian,mode="full")
    #     halfk=kernelsize/2
    #     fillfrom=np.arange(halfk)
    #     refill=(halfk-fillfrom)+halfk
    #     windsm[:,:,refill]+=windsm[:,:,fillfrom]
    #     windsm[:,refill,:]+=windsm[:,fillfrom,:]
    #     windsm[:,:,0-refill]+=windsm[:,:,0-fillfrom]
    #     windsm[:,0-refill,:]+=windsm[:,0-fillfrom,:]
    #     wind=windsm[:,halfk:-halfk,halfk:-halfk]
    #     return wind
    # def smooth_wind2(self,wind,chalfk):
    #     windout=wind.copy()
    #     for i in range(wind.shape[0]):
    #         for j in range(chalfk,wind.shape[1]-chalfk):
    #             for k in range(chalfk,wind.shape[2]-chalfk):
    #                 windout[i,j,k]=np.mean(wind[i,j-chalfk:j+chalfk,k-chalfk:k+chalfk])
    #     return windout

    def _next_surface(self):
        geoLUT=self._geoLUT
        x=geoLUT[:,:,:,1]
        y=geoLUT[:,:,:,0]
        w=geoLUT[:,:,:,2]
        # io_ratio = 16 = ratio of NARR grid size (32km) to WRF grid size (2km)
        offset=np.int(np.ceil(halfk/io_ratio))
        minx=x.min().astype('i')-offset 
        maxx=x.max().astype('i')+offset
        miny=y.min().astype('i')-offset
        maxy=y.max().astype('i')+offset
        curx=(x-minx).astype('i')
        cury=(y-miny).astype('i')
        d=swim_io.Dataset(self._sfc_files[self._curfile], 'r')
        Nxy=curx.shape[0:2]
        sensible_heat=np.zeros(Nxy,dtype=np.float32)
        latent_heat=np.zeros(Nxy,dtype=np.float32)
        longwave_up=np.zeros(Nxy,dtype=np.float32)
        albedo=np.zeros(Nxy,dtype=np.float32)
        
        gaussian=gauss_kern(2)
        gaussian2=gauss_kern(np.floor(2*offset))
        
        curdata=convolve(d.variables[self.shvar][miny:maxy,minx:maxx],gaussian,mode="same")
        for i in range(4):sensible_heat+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=convolve(d.variables[self.lhvar][miny:maxy,minx:maxx],gaussian,mode="same")
        for i in range(4):latent_heat+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=convolve(d.variables[self.LWuvar][miny:maxy,minx:maxx],gaussian,mode="same")
        for i in range(4):longwave_up+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=convolve(d.variables[self.albvar][miny:maxy,minx:maxx],gaussian,mode="same")
        for i in range(4):albedo+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        d.close()
        sensible_heat=sensible_heat[halfk:-halfk,halfk:-halfk]
        latent_heat=latent_heat[halfk:-halfk,halfk:-halfk]
        longwave_up=longwave_up[halfk:-halfk,halfk:-halfk]
        albedo=albedo[halfk:-halfk,halfk:-halfk]
        return Bunch(sensible_heat=-sensible_heat,latent_heat=-latent_heat,
                     longwave_up=longwave_up,albedo=albedo/100.0)

        
    
    def nextBilin(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        geoLUT=self._geoLUT
        x=geoLUT[:,:,:,1]
        y=geoLUT[:,:,:,0]
        w=geoLUT[:,:,:,2]
        offset=np.int(np.ceil(halfk/io_ratio))
        minx=x.min().astype('i')-offset # 16 = ratio of NARR grid size (32km) to WRF grid size (2km)
        maxx=x.max().astype('i')+offset
        miny=y.min().astype('i')-offset
        maxy=y.max().astype('i')+offset
        curx=(x-minx).astype('i')
        cury=(y-miny).astype('i')

        d=swim_io.Dataset(self._filenames[curfile], 'r')
        pressure=d.variables['lv_ISBL0'][:]*100 # convert hPa to Pa
#        evenpressures=np.where(pressure%5000 ==0)[0][4:-2]
        evenpressures=np.where(pressure>000)[0]
        print(evenpressures.shape)
        pressure=pressure[evenpressures]
        Nz=evenpressures.size
        Nxy=curx.shape[0:2]
        N=[Nz,Nxy[0],Nxy[1]]
        
        temperature=np.zeros(N,dtype=np.float32)
        hgt=np.zeros(N,dtype=np.float32)
        qc=np.zeros(N,dtype=np.float32)
        # qi=np.zeros(N,dtype=np.float32) #for now qi is added to qc so we don't have to worry about ni
        specific_humidity=np.zeros(N,dtype=np.float32)
        wind_u=np.zeros(N,dtype=np.float32)
        wind_v=np.zeros(N,dtype=np.float32)
        gaussian=gauss_kern(2)
        gaussian2=gauss_kern(np.floor(2*offset))
        # Note, this currently reads all pressure levels from disk then subsets to "evenpressures"
        curdata=self.convolve2to3(d.variables[self.tvar][:,miny:maxy,minx:maxx][evenpressures,:,:],gaussian,mode="same")
        for i in range(4):temperature+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=self.convolve2to3(d.variables[self.hvar][:,miny:maxy,minx:maxx][evenpressures,:,:],gaussian,mode="same")
        for i in range(4):hgt+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=self.convolve2to3(d.variables[self.qcvar][:,miny:maxy,minx:maxx][evenpressures,:,:],gaussian,mode="same")
        for i in range(4):qc+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=self.convolve2to3(d.variables[self.qivar][:,miny:maxy,minx:maxx][evenpressures,:,:],gaussian,mode="same")
        for i in range(4):qc+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=self.convolve2to3(d.variables[self.qvvar][:,miny:maxy,minx:maxx][evenpressures,:,:],gaussian,mode="same")
        for i in range(4):specific_humidity+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
#        curdata=fast_mean.fast_smooth(d.variables[self.uvar][:,miny:maxy,minx:maxx][evenpressures,:,:],offset)
        curdata=d.variables[self.uvar][:,miny:maxy,minx:maxx][evenpressures,:,:]
        for i in range(4):wind_u+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
#        curdata=fast_mean.fast_smooth(d.variables[self.vvar][:,miny:maxy,minx:maxx][evenpressures,:,:],offset)
        curdata=d.variables[self.vvar][:,miny:maxy,minx:maxx][evenpressures,:,:]
        for i in range(4):wind_v+=np.float32(curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        d.close()
        pressure=pressure[:,np.newaxis,np.newaxis].repeat(N[1],axis=1).repeat(N[2],axis=2).astype(np.float32)
        datestr=self._filenames[curfile].split('.')[1]
        kernelsize=np.int(np.floor(halfk))
        gaussian=gauss_kern(kernelsize)
#        wind_u=fast_mean.fast_smooth(wind_u,halfk)[:,halfk:-halfk,halfk:-halfk]
#        wind_v=fast_mean.fast_smooth(wind_v,halfk)[:,halfk:-halfk,halfk:-halfk]
        wind_u=wind_u[:,halfk:-halfk,halfk:-halfk]
        wind_v=wind_v[:,halfk:-halfk,halfk:-halfk]
        # print(time.time()-t1)
        temperature=temperature[:,halfk:-halfk,halfk:-halfk]
        pressure=pressure[:,halfk:-halfk,halfk:-halfk]
        specific_humidity=specific_humidity[:,halfk:-halfk,halfk:-halfk]
        hgt=hgt[:,halfk:-halfk,halfk:-halfk]
        qc=qc[:,halfk:-halfk,halfk:-halfk]
        potential_temperature=temperature*(100000.0/pressure)**(R/cp)

        if self._sfc_files!=None:
            sfc=self._next_surface()
        else:
            sfc=None

        self._curfile+=1
        return Bunch(ta=temperature, p=pressure,th=potential_temperature,
                     sh=specific_humidity, hgt=hgt,qc=qc,
                     qv=specific_humidity/(1-specific_humidity), 
                     u=wind_u, v=wind_v,w=np.zeros(N,dtype=np.float32,order="F"),
                     date=datestr,sfc=sfc,
                     qr=np.zeros(N,dtype=np.float32,order="F"),nr=np.zeros(N,dtype=np.float32,order="F"),
                     qs=np.zeros(N,dtype=np.float32,order="F"),qi=np.zeros(N,dtype=np.float32,order="F"),
                     ni=np.zeros(N,dtype=np.float32,order="F"),qg=np.zeros(N,dtype=np.float32,order="F"),
                     )
                     # dp=units.sh2dp(temperature,pressure,specific_humidity), 
                     # rh=relative_humidity, 

    def next(self):
        if self._geoLUT==None:
            return self.nextNN()
        else:
            return self.nextBilin()
            
    def close(self):
        pass
    def __enter__(self):
        return self
    
    def __exit__(self):
        self.close()
    
    def __del__(self):
        self.close()
        
def topo_adjust_weather(hitopo,lowtopo,weather):
    # first find the first pressure level that is above the coarse topography
    zlevel=np.zeros(weather.hgt.shape[1:],dtype='i')
    weights=np.zeros(weather.hgt.shape[1:])
    dhgt= -np.gradient(weather.hgt)[0]
    # dz3d=weather.hgt-hitopo[np.newaxis,:,:].repeat(weather.hgt.shape[0],axis=0)
    for i in range(weather.hgt.shape[0]):
        curdz=weather.hgt[i,:,:]-lowtopo
        curlevel=np.where(curdz>0)
        if len(curlevel[0])>0:
            zlevel[curlevel]=i
            weights[curlevel]=curdz[curlevel]/(dhgt[i,:,:][curlevel])
    N=hitopo.shape
    # and subset all variables to the bottom n layers above the topography
    nlevels=8
    print(weather.hgt.shape)
    x=np.arange(N[1])[np.newaxis,:].repeat(N[0],axis=0)
    y=np.arange(N[0])[:,np.newaxis].repeat(N[1],axis=1)
    for key in weather.keys():
        if (key!="date") and (key!="sfc"):
            tmpoutput=np.zeros((nlevels,N[0],N[1]),dtype=np.float32)
            for level in range(nlevels):
                # tmpoutput[level,...]=weather[key][zlevel-level,y,x]
                tmpoutput[level,...]=(weather[key][zlevel-level-1,y,x]*(1-weights)
                                     +weather[key][zlevel-level,y,x]*(weights))
            weather[key]=tmpoutput
    # then adjust all pressures by dz (difference between coarse topo and high res topo)
    dz=hitopo-lowtopo
    weather.p=convert_p(weather.p,weather.hgt,dz[np.newaxis,:,:].repeat(weather.p.shape[0],axis=0))


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)    
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)               
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    g=g[size/2:-size/2,sizey/2:-sizey/2]
    return g / g.sum()


def update_weather(newatm,oldatm):
    for k in newatm.keys():
        if (k!="hgt") and (k!="date") and (k!="sfc"):
            newatm[k][:,1:-1,1:-1]=oldatm[k][:,1:-1,1:-1]

# Force vertical convergence/divergence to balance horizontal convergence/divergence
def adjust_winds(u,v,w):
    ustag=(u[:,:,1:]+u[:,:,:-1])/2#[4:-2]
    vstag=(v[:,1:,:]+v[:,:-1,:])/2
    conv=(ustag[:,1:-1,1:]-ustag[:,1:-1,:-1])+(vstag[:,1:,1:-1]-vstag[:,:-1,1:-1])
    w[0,1:-1,1:-1]=-conv[0,:,:]
    for i in range(w.shape[0]-1):
        w[i+1,1:-1,1:-1]=w[i,1:-1,1:-1]-conv[i+1,:,:]
    return(ustag,vstag,w)

file_search="nc3d/merged*.nc"
sfc_file_search="ncsfc/merged*.nc"
if wrfres==2000:
    topofile='/glade/scratch/gutmann/usbr/narr_data/baseline_info/2km_wrf_input_d01'
elif wrfres==4000:
    topofile='/glade/scratch/gutmann/usbr/narr_data/baseline_info/4km_wrf_output.nc'
outputdir='output/'
def main(): # (file_search="nc3d/merged*.nc",topofile='/d2/gutmann/usbr/narr_data/baseline_info/2km_wrf_input_d01'):
    topoinfo=swim_io.read_nc(topofile,var='HGT')
    if sub_x1==None:
        xendpt=-halfk
    else:
        xendpt=sub_x1-halfk
    if sub_y1==None:
        yendpt=-halfk
    else:
        yendpt=sub_y1-halfk
    topo=topoinfo.data[0,sub_y0+halfk:yendpt,sub_x0+halfk:xendpt]
    (Ny,Nx)=topo.shape
    # for use with the linear theory winds eventually
    # topo=np.choose(topo<2000,(topo,2000))
    # topo[:,:32]=2000
    # Fzs=fft.fftshift(fft.fft2(topo))/(Nx*Ny)

    gkern=gauss_kern(40)
    narr=NARR_Reader(file_search,sfc=sfc_file_search,bilin=True,nn=False)
    timestep=3.0*60.0*60.0 #3hrs
    oldweather=None
    print("Off and running...")
    weather=narr.next()
    topo_adjust_weather(topo, narr.topo, weather)
    (weather.u,weather.v,weather.w)=adjust_winds(weather.u.astype(np.float64),
                                                weather.v.astype(np.float64),
                                                weather.w.astype(np.float64))
    (Nz,Ny,Nx)=weather.qv.shape
    oldt=weather.th.copy()
    # print("convolving...")
    # for i in range(weather.u.shape[0]):
    #     tmpu=convolve(weather.u[i,:,:],gkern,mode="same")
    #     weather.u[i,:,:]=tmpu
    #     tmpv=convolve(weather.v[i,:,:],gkern,mode="same")
    #     weather.v[i,:,:]=tmpv
    oldweather=weather
    print("Initializing microphysics...")
    swim_lib.swim_step.init() #calls thompson_init internally

    # processPool=Pool(NPROCESSORS)
    processPool=None
    print("Reading NARR")
    for weather in narr:
        # for i in range(weather.u.shape[0]):
        #     tmpu=convolve(weather.u[i,:,:],gkern,mode="same")
        #     weather.u[i,:,:]=tmpu
        #     tmpv=convolve(weather.v[i,:,:],gkern,mode="same")
        #     weather.v[i,:,:]=tmpv
        topo_adjust_weather(topo, narr.topo, weather)
        newt=weather.th.copy()
        dTdt=(newt-oldt)/timestep
        oldt=newt
        newp=weather.p.copy()
        (weather.u,weather.v,weather.w)=adjust_winds(weather.u.astype(np.float64),weather.v.astype(np.float64),weather.w.astype(np.float64))
        newu=weather.u.copy()
        newv=weather.v.copy()
        neww=weather.w.copy()
        update_weather(weather,oldweather)
        weather.u=newu[4:-2]
        weather.v=newv
        weather.w=neww
        # P=swim.swim2d(topo,weather,mp, timestep=timestep,fname="outputtest_"+weather.date+".nc")
        print(" ")
        print("-------------------------------------------------")
        print(weather.date)
        print("-------------------------------------------------")
        # print("pressure range:",weather.p.max(),weather.p.min())
        # print("hgt range:",weather.hgt.max(),weather.hgt.min())
        P=swim.swim2d(topo,weather,oldweather,newp,newu,newv,neww,swim_lib,dTdt,processPool,timestep=timestep,dx=wrfres)
        #,fname="outputtest_"+weather.date+".nc")
        curdate=np.int64(weather.date)
        ncpout=swim_io.NC_writer(outputdir+'swim_p_'+str(curdate),Nx,Ny,var='PRECIP')
        nctout=swim_io.NC_writer(outputdir+'swim_t_'+str(curdate),Nx,Ny,Nz,var='TEMP')
        ncqout=swim_io.NC_writer(outputdir+'swim_qv_'+str(curdate),Nx,Ny,Nz,var='qv')
        ncpout.appendToVar(P, date=curdate)
        nctout.appendToVar(weather.th, date=curdate)
        ncqout.appendToVar(weather.qv, date=curdate)
        nctout.close()
        ncpout.close()
        ncqout.close()
        oldweather=weather
    narr.close()
    # processPool.close()




if __name__ == '__main__':
    main()
