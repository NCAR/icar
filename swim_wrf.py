#!/usr/bin/env python

"""
SYNOPSIS

    swim_wrf.py [-h] [--verbose] [-v, --version] [WRF_filename] [hi_res_topofile]

DESCRIPTION
    Reads forcing data to drive the: 
        Simple Weather Interpolation Model  (SWIM)
        Physics based Orographic Precipitation model (POP)
        Hybrid Orographic Precipitation model (HOP)
        or similar silly acronyms / abbreviations (SSAA)

EXAMPLES

    swim_wrf.py

EXIT STATUS

    None

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    0.6   - adding WRF high res winds
    0.5.1 - SPEED (esp. with new microphysics ~50x faster)
    0.5   - added LT winds & fancy 3D grid(?)
    0.4   - adding PBL mixing and (stupid) radiative cooling 1.5K/day
    0.3   - WRF based (previously from swim_narr.py)
    0.2   - the hey it sort of works version!
    0.1   - skeletal, pieces work
"""
# Python standard modules
import sys
import os
import traceback
#import argparse # requires python >=2.7
from multiprocessing import Pool,Queue,Process
import glob
import time
# Standard 3rd party modules
import numpy.fft as fft
from scipy.signal import convolve
import numpy as np
# my common modules
import units
from bunch import Bunch
# swim specific modules
import swim_io # this class combines the functionality of the three file IO classes below
import swim as swim_lib #this is the meet of the fortran code
import swim_fast as swim # a faster version of the hand off between python and fortran
import fast_mean # inline C to do a running spatial mean
import lt_winds #computes Linear Theory 3D winds


NPROCESSORS=2 #not used currently
R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2
physics=int(1) #1=thompson other=simple
# use_linear_winds=False #flag to use linear theory winds or not
use_linear_winds=True  #comment out one or the other
use_wrf_winds=False # ditto for wrf_based winds
# use_wrf_winds=True

# this is the half kernel window size for smoothing applied to the input data (wind)
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
sub_y0=30
sub_y1=-30
sub_x0=30
sub_x1=-30
wrfres=4000.0
io_ratio=36000.0/wrfres #this is the ratio of the input grid (36km) to the output grid (4km)

file_search="forcing/wrfout_d01_200*00"
wind_files="winds/wrfout_d01_200*00"
sfc_file_search=file_search
base_dir="baseline_info/"
if wrfres==2000:
    # topofile='/d2/gutmann/usbr/narr_data/baseline_info/2km_wrf_input_d01'
    topofile=base_dir+'2km_wrf_input_d01'
elif wrfres==4000:
    # topofile='/d2/gutmann/usbr/narr_data/baseline_info/4km_wrf_output.nc'
    topofile=base_dir+'4km_wrf_output.nc'
outputdir='output/'


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
    N1=lat1.shape
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
        ymin=lasty-winhalfsize
        if ymin<0:ymin=0
        ymax=lasty+winhalfsize
        if ymax>=N1[0]:ymax=N1[0]-1
        xmin=lastx-winhalfsize
        if xmin<0:xmin=0
        xmax=lastx+winhalfsize
        if xmax>=N1[1]:xmax=N1[1]-1
        latwin=lat1[ymin:ymax,xmin:xmax]
        lonwin=lon1[ymin:ymax,xmin:xmax]
        for j in range(N[0]):
            if (prevx!=lastx) or (prevy!=lasty):
                ymin=lasty-winhalfsize
                if ymin<0:ymin=0
                ymax=lasty+winhalfsize
                if ymax>=N1[0]:ymax=N1[0]-1
                xmin=lastx-winhalfsize
                if xmin<0:xmin=0
                xmax=lastx+winhalfsize
                if xmax>=N1[1]:xmax=N1[1]-1
                latwin=lat1[ymin:ymax,xmin:xmax]
                lonwin=lon1[ymin:ymax,xmin:xmax]
                (prevx,prevy)=(lastx,lasty)
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            lastx=xmin+newx
            lasty=ymin+newy
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
            x[x<0]=0
            y[y<0]=0
            x[x>=(xmax-xmin)]=xmax-xmin-1
            y[y>=(ymax-ymin)]=ymax-ymin-1
                    
            # bilinear interpolation for an arbitrary grid spacing 
            # (must be a grid, but this will handle a lot of irregularities)
            weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])

            out[j,i,:,0]=y+ymin
            out[j,i,:,1]=x+xmin
            out[j,i,:,2]=weights
        
    return out
              
                
def convert_p(p,h,dz):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals
    # h,dz in meters
    slp = p/(1 - 2.25577E-5*h)**5.25588
    pout=slp*(1 - 2.25577E-5*(h+dz))**5.25588
    return pout

def fix_top_bottom(topo,dz):
    # find the worst dz between the top and bottom of the dataset
    worst_offset=np.max(np.abs(topo[0,:]-topo[-1,:]))
    # find the minimum distance required to smoothly transition between those
    # two endpoints
    mindist=np.round(worst_offset/(dz/50.0))
    # calculate the mid point / average topography that will be the new borders
    aves=(topo[0,:]+topo[-1,:])/2.0
    # create the output topography array
    sz=topo.shape
    if mindist%2==1:
        mindist+=1
    outputtopo=np.zeros((sz[0]+mindist,sz[1]))
    # fill in everything inside the padding zone with the original topography
    padding=np.round(mindist/2.0)
    outputtopo[padding:sz[0]+padding,:]=topo
    # now calculate a fading factor (0-1) between the actual topography and the new boundary
    fade=(np.arange(padding)/np.float(padding)).reshape((padding,1)).repeat(sz[1],axis=1)
    # finally calculate the topography over the padding zones by linearly combining the 
    # aves / mid-point with the actual topography fading from 100% ave to 100% real topo
    border=aves.reshape((1,-1)).repeat(padding,axis=0)
    topo_top=topo[0,:].reshape((1,-1)).repeat(padding,axis=0)
    topo_bottom=topo[-1,:].reshape((1,-1)).repeat(padding,axis=0)
    
    outputtopo[:padding,:]=(1-fade)*border + fade*topo_top
    outputtopo[-padding:,:]=fade*border + (1-fade)*topo_bottom
    return outputtopo,padding

    
def topo_preprocess(topo,test=False):
    acceptable_topo_diff=np.max(np.diff(topo))/3.0
    vpad=0.0
    hpad=0.0
    
    if test:
        topo=topo.copy() #dont change the topo that is outside this routine
        t1=np.arange(16.0,-1,-1.0)[np.newaxis,:].repeat(160,axis=0)/18.0
        zt=topo[91:251,:17]
        topo[91:251,:17]=1500*t1+zt*(1-t1)
    
        t1=np.arange(13.0)[:,np.newaxis].repeat(40,axis=1)/12.0
        zt=topo[250:,78:118]
        topo[250:,78:118]=2200*t1+zt*(1-t1)
    
        topo[245:,:50]=1500
        topo[topo<1500]=1500
    
    if np.max(np.abs(topo[0,:]-topo[-1,:]))>acceptable_topo_diff: # more a less a given
        (topo,vpad)=fix_top_bottom(topo,acceptable_topo_diff)
    if np.max(np.abs(topo[:,0]-topo[:,-1]))>acceptable_topo_diff: # ditto
        (topo,hpad)=fix_top_bottom(topo.T,acceptable_topo_diff)
        topo=topo.T
    #this might matter more than I thought... and maybe all edges should get exactly to 0?
    topo-=topo.min() #1500# (topo[0,:].mean()+topo[-1,:].mean()+topo[:,0].mean()+topo[:,-1].mean())/4.0
    return topo,vpad,hpad
    

class WRF_Reader(object):
    _filenames=None
    _sfc_files=None
    _curfile=0
    curpos=0
    x=-99
    y=-99
    _bilin=False
    _nn=True
    _geoLUT=None
    topo=None
    # atm variables in WRF 3d files
    timevar="Times"      # time...                                  1D
    hvar='HGT'           #terran height(?)                  [m]     2D
    tvar='T'             #potential temperature (-300)      [K]     3D
    gphvar_base="PHB"    #base geopotential height          [m]     3D
    gphvar_pert="PH"     #geopotential height perturbation  [m]     3D
    pvar_base="PB"       #base pressure                     [pa]    3D
    pvar_pert="P"        # pressure perturbation            [pa]    3D
    qcvar='QCLOUD'       # cloud water mixing ratio         [kg/kg] 3D
    qvvar='QVAPOR'       # water vapor mixing ratio         [kg/kg] 3D
    qivar='QICE'         # cloud ice mixing ratio           [kg/kg] 3D
    wvar='W'             # vertical winds                   [pa/s]  3D NOT USED
    uvar='U'             # U wind component                 [m/s]   3D 
    vvar='V'             # V wind component                 [m/s]   3D 
    # land surface variables in WRF files
    pblvar='PBLH'        # PBL height                       [m]     2D
    lhvar="LH"           # Latent heat flux from surface    [W/m^2] 2D
    shvar="HFX"          # Sensible heat flux from surface  [W/m^2] 2D

    def make_model_domain(self,nlevels,maxheight,topo):
        sz=topo.shape
        
        slopes=(maxheight-topo)/nlevels
        if use_linear_winds:
            slopes=slopes.mean()
            self.hgt3d=(topo[np.newaxis,...].repeat(nlevels,axis=0)+
                    slopes[np.newaxis,...]*np.arange(nlevels)[:,np.newaxis,np.newaxis])
        elif use_wrf_winds:
            self.hgt3d=(swim_io.read_nc(self._filenames[0],var=self.gphvar_base)
                    .data[0,self.usepressures,self.sub_y0:sub_y1,sub_x0:sub_x1]/9.81)
        else:
            self.hgt3d=(topo[np.newaxis,...].repeat(nlevels,axis=0)+
                    slopes[np.newaxis,...].repeat(nlevels,axis=0)*
                    np.arange(nlevels)[:,np.newaxis].repeat(sz[0]*sz[1],axis=1).reshape((nlevels,sz[0],sz[1])))
        
    
    def init_xy(self,driverfilename='forcing/wrfout_d01_2000-10-01_00:00:00'):
        if wrfres==2000:
            wrffilename=bas_dir+'2km_wrf_input_d01'
        elif wrfres==4000:
            wrffilename=base_dir+'4km_wrf_output.nc'
        d=swim_io.Dataset(driverfilename, 'r')
        nlat=d.variables['XLAT'][0,:,:]
        nlon=d.variables['XLONG'][0,:,:]
        nulat=d.variables['XLAT_U'][0,:,:]
        nulon=d.variables['XLONG_U'][0,:,:]
        nvlat=d.variables['XLAT_V'][0,:,:]
        nvlon=d.variables['XLONG_V'][0,:,:]
        hgt=d.variables['HGT'][0,:,:]
        self.nlevels=20
        self.usepressures=np.arange(self.nlevels)
        self.base_pressure=d.variables[self.pvar_base][0,...][self.usepressures,...]
        self.base_gph=d.variables[self.gphvar_base][0,...][self.usepressures,...]
        maxheight=self.base_gph[-1,:,:]/9.8
        d.close()
        wlat=swim_io.read_nc(wrffilename,var='XLAT').data[0,sub_y0:sub_y1,sub_x0:sub_x1]
        wlon=swim_io.read_nc(wrffilename,var='XLONG').data[0,sub_y0:sub_y1,sub_x0:sub_x1]
        hires_topo=swim_io.read_nc(wrffilename,var='HGT').data[0,sub_y0:sub_y1,sub_x0:sub_x1][halfk:-halfk,halfk:-halfk]
        print('Calculating XY match lookup table, this may take a while.')
        if self._nn:
            wlat=wlat[halfk:-halfk,halfk:-halfk]
            wlon=wlon[halfk:-halfk,halfk:-halfk]
            (x,y)=match_xy(nlat,nlon,wlat,wlon)
            (xu,yu)=match_xy(nulat,nulon,wlat,wlon)
            (xv,yv)=match_xy(nvlat,nvlon,wlat,wlon)
            self.x=x.astype('i')
            self.y=y.astype('i')
            self.xu=xu.astype('i')
            self.yu=yu.astype('i')
            self.xv=xv.astype('i')
            self.yv=yv.astype('i')
            self.topo=hgt[self.y,self.x]
            maxheight=maxheight[self.y,self.x]
            self.base_gph=self.base_gph[:,self.y,self.x]
            self.base_pressure=self.base_pressure[:,self.y,self.x]
            self.make_model_domain(self.nlevels,maxheight,hires_topo)
            
        elif self._bilin:
            self._geoLUT=match_xy_bilin(nlat,nlon,wlat,wlon)
            self._geoLUTu=match_xy_bilin(nulat,nulon,wlat,wlon)
            self._geoLUTv=match_xy_bilin(nvlat,nvlon,wlat,wlon)
            curx=(self._geoLUT[:,:,:,1]).astype('i')
            cury=(self._geoLUT[:,:,:,0]).astype('i')
            w=self._geoLUT[:,:,:,2]
            self.topo=np.zeros(curx.shape[0:2],dtype=np.float32)
            for i in range(4):self.topo+=hgt[cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            self.topo=self.topo[halfk:-halfk,halfk:-halfk]
            MH=np.zeros(curx.shape[0:2],dtype=np.float32)
            for i in range(4):MH+=maxheight[cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            MH=MH[halfk:-halfk,halfk:-halfk]
            self.make_model_domain(self.nlevels,MH,hires_topo)
            
        
    def __init__(self, file_search,sfc=None,nn=True, bilin=False, windonly=False, *args, **kwargs):
        super(WRF_Reader,self).__init__(*args, **kwargs)
        self._filenames=np.sort(glob.glob(file_search))
        if sfc!=None:
            self._sfc_files=glob.glob(sfc)
        if nn:
            self._nn=True
            self._bilin=False
        elif bilin:
            self._nn=False
            self._bilin=True
        self.windonly=windonly
        self.init_xy(self._filenames[0])
        self.curpos=19*24
        if windonly:
            self.curpos=0
        # d=swim_io.Dataset(self._filenames[self._curfile], 'r')
        # npos=d.variables[self.vvar].shape[0]-1
        # self.curpos=npos-10
        # d.close()
    
    
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

        minxu=self.xu.min()
        maxxu=self.xu.max()+1
        minyu=self.yu.min()
        maxyu=self.yu.max()+1
        curxu=self.xu-minxu
        curyu=self.yu-minyu

        minxv=self.x.min()
        maxxv=self.x.max()+1
        minyv=self.y.min()
        maxyv=self.y.max()+1
        curxv=self.x-minxv
        curyv=self.y-minyv
        plevels=self.usepressures
        d=swim_io.Dataset(self._filenames[curfile], 'r')
        if self.windonly:
            wind_u=d.variables[self.uvar][self.curpos,:,minyu:maxyu,minxu:maxxu][plevels,...]
            wind_u=wind_u[:,curyu,curxu]
            wind_v=d.variables[self.vvar][self.curpos,:,minyv:maxyv,minxv:maxxv][plevels,...]
            wind_v=wind_v[:,curyv,curxv]
            hgt=d.variables[self.gphvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            hgt=(hgt[:,cury,curx]+self.base_gph)/9.8
            temperature=d.variables[self.tvar][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            temperature=temperature[:,cury,curx]+300
            specific_humidity=d.variables[self.qvvar][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            specific_humidity=specific_humidity[:,cury,curx]
            pressure=d.variables[self.pvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            pressure=pressure[:,cury,curx]+self.base_pressure
            # temperature=None
            # pressure=None
            # specific_humidity=None
            relative_humidity=None
        else:
            temperature=d.variables[self.tvar][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]
            pressure=d.variables[self.pvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]
            specific_humidity=d.variables[self.qvvar][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]
            relative_humidity=specific_humidity.copy()
            wind_u=d.variables[self.uvar][self.curpos,:,minyu:maxyu,minxu:maxxu][plevels,curyu,curxu]
            wind_v=d.variables[self.vvar][self.curpos,:,minyv:maxyv,minxv:maxxv][plevels,curyv,curxv]
            hgt=(d.variables[self.gphvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]+self.base_gph)/9.8

        datestr=''.join(d.variables[self.timevar][self.curpos])
        d.close()


        self.time_inc()
        N=list(wind_v.shape)
        N[2]+=1 #v is staggered in y direction
        return Bunch(th=temperature, p=pressure,hgt=hgt,
                     qv=specific_humidity/(1-specific_humidity), 
                     rh=relative_humidity, 
                     u=wind_u, v=wind_v,w=np.zeros(N,dtype=np.float32,order="F"),
                     date=datestr)
    
    def convolve2to3(self,data,kern,mode="same"):
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
    

    def _next_surface(self,curpos):
        geoLUT=self._geoLUT
        x=geoLUT[:,:,:,1]
        y=geoLUT[:,:,:,0]
        w=geoLUT[:,:,:,2]
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
        pblh=np.zeros(Nxy,dtype=np.float32)
        # longwave_up=np.zeros(Nxy,dtype=np.float32)
        # albedo=np.zeros(Nxy,dtype=np.float32)
        
        gaussian=gauss_kern(2)
        
        curdata=d.variables[self.shvar][curpos,miny:maxy,minx:maxx]
        for i in range(4):sensible_heat+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=d.variables[self.lhvar][curpos,miny:maxy,minx:maxx]
        for i in range(4):latent_heat+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=d.variables[self.pblvar][curpos,miny:maxy,minx:maxx]
        for i in range(4):pblh+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        # curdata=convolve(d.variables[self.LWuvar][miny:maxy,minx:maxx],gaussian,mode="same")
        # for i in range(4):longwave_up+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        # curdata=convolve(d.variables[self.albvar][miny:maxy,minx:maxx],gaussian,mode="same")
        # for i in range(4):albedo+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        d.close()
        sensible_heat=sensible_heat[halfk:-halfk,halfk:-halfk]/6
        latent_heat=latent_heat[halfk:-halfk,halfk:-halfk]/6
        pblh=pblh[halfk:-halfk,halfk:-halfk]
        # print(sensible_heat.max())
        # print(latent_heat.max(),sensible_heat.max())
        return Bunch(sensible_heat=sensible_heat,latent_heat=latent_heat,pblh=pblh)
        # longwave_up=longwave_up[halfk:-halfk,halfk:-halfk]
        # albedo=albedo[halfk:-halfk,halfk:-halfk]
        # return Bunch(sensible_heat=-sensible_heat,latent_heat=-latent_heat,
        #              longwave_up=longwave_up,albedo=albedo/100.0)

        
    
    def nextBilin(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        x=self._geoLUT[:,:,:,1]
        y=self._geoLUT[:,:,:,0]
        w=self._geoLUT[:,:,:,2]
        xu=self._geoLUTu[:,:,:,1]
        yu=self._geoLUTu[:,:,:,0]
        wu=self._geoLUTu[:,:,:,2]
        xv=self._geoLUTv[:,:,:,1]
        yv=self._geoLUTv[:,:,:,0]
        wv=self._geoLUTv[:,:,:,2]
        offset=np.int(np.ceil(halfk/io_ratio))
        minx=x.min().astype('i')-offset
        maxx=x.max().astype('i')+offset
        miny=y.min().astype('i')-offset
        maxy=y.max().astype('i')+offset
        curx=(x-minx).astype('i')
        cury=(y-miny).astype('i')
        minux= xu.min().astype('i')-offset
        maxux= xu.max().astype('i')+offset
        minuy= yu.min().astype('i')-offset
        maxuy= yu.max().astype('i')+offset
        curux=(xu-minux).astype('i')
        curuy=(yu-minuy).astype('i')
        minvx= xv.min().astype('i')-offset
        maxvx= xv.max().astype('i')+offset
        minvy= yv.min().astype('i')-offset
        maxvy= yv.max().astype('i')+offset
        curvx=(xv-minvx).astype('i')
        curvy=(yv-minvy).astype('i')

        d=swim_io.Dataset(self._filenames[curfile], 'r')
        
        datestr=''.join(d.variables[self.timevar][self.curpos])
        usepressures=self.usepressures
        
        Nz=usepressures.size
        Nxy=curx.shape[0:2]
        N=[Nz,Nxy[0],Nxy[1]]
        
        wind_u=np.zeros(N,dtype=np.float32,order="F")
        wind_v=np.zeros(N,dtype=np.float32,order="F")
        
        if self.windonly:
            curdata=fast_mean.fast_smooth(d.variables[self.uvar][self.curpos,...][:,minuy:maxuy,minux:maxux][usepressures,:,:],offset)
            for i in range(4):wind_u+=curdata[:,curuy[:,:,i],curux[:,:,i]]*wu[:,:,i]
            # curdata=d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:]
            curdata=fast_mean.fast_smooth(d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:],offset)
            for i in range(4):wind_v+=curdata[:,curvy[:,:,i],curvx[:,:,i]]*wv[:,:,i]
            # we also want geopotential height so we can properly interpolate vertically if necessary
            curdata=(d.variables[self.gphvar_pert][self.curpos,...][usepressures,:,:]+self.base_gph)[:,miny:maxy,minx:maxx]
            for i in range(4):hgt+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            hgt=hgt[:,halfk:-halfk,halfk:-halfk]/9.8
            potential_temperature=None
            qc=None
            specific_humidity=None
            pressure=None
            qr=None
            nr=None
            qs=None
            qi=None
            ni=None
            qg=None
        else:
            potential_temperature=np.zeros(N,dtype=np.float32,order="F")
            hgt=np.zeros(N,dtype=np.float32,order="F")
            qc=np.zeros(N,dtype=np.float32,order="F")
            # qi=np.zeros(N,dtype=np.float32) #for now qi is added to qc so we don't have to worry about ni
            specific_humidity=np.zeros(N,dtype=np.float32,order="F")
            pressure=np.zeros(N,dtype=np.float32,order="F")
            
            # Note, this currently reads all pressure levels from disk then subsets to "usepressures"
            curdata=d.variables[self.pvar_pert][self.curpos,...]
            curdata=(curdata[usepressures,:,:]+self.base_pressure)[:,miny:maxy,minx:maxx]
            for i in range(4):pressure+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.tvar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]+300
            for i in range(4):potential_temperature+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=(d.variables[self.gphvar_pert][self.curpos,...][usepressures,:,:]+self.base_gph)[:,miny:maxy,minx:maxx]
            for i in range(4):hgt+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.qcvar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]
            for i in range(4):qc+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.qivar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]
            for i in range(4):qc+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.qvvar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]
            for i in range(4):specific_humidity+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            # curdata=d.variables[self.uvar][self.curpos,...][:,minuy:maxuy,minux:maxux][usepressures,:,:]
            curdata=fast_mean.fast_smooth(d.variables[self.uvar][self.curpos,...][:,minuy:maxuy,minux:maxux][usepressures,:,:],offset)
            for i in range(4):wind_u+=curdata[:,curuy[:,:,i],curux[:,:,i]]*wu[:,:,i]
            # curdata=d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:]
            curdata=fast_mean.fast_smooth(d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:],offset)
            for i in range(4):wind_v+=curdata[:,curvy[:,:,i],curvx[:,:,i]]*wv[:,:,i]

            potential_temperature=potential_temperature[:,halfk:-halfk,halfk:-halfk]
            pressure=pressure[:,halfk:-halfk,halfk:-halfk]
            specific_humidity=specific_humidity[:,halfk:-halfk,halfk:-halfk]
            hgt=hgt[:,halfk:-halfk,halfk:-halfk]/9.8
            qc=qc[:,halfk:-halfk,halfk:-halfk]
            
            N2=[Nz,Nxy[0]-halfk*2,Nxy[1]-halfk*2]
            qr=np.zeros(N2,dtype=np.float32,order="F")
            nr=np.zeros(N2,dtype=np.float32,order="F")
            qs=np.zeros(N2,dtype=np.float32,order="F")
            qi=np.zeros(N2,dtype=np.float32,order="F")
            ni=np.zeros(N2,dtype=np.float32,order="F")
            qg=np.zeros(N2,dtype=np.float32,order="F")
            
        
        self.curpos+=1
        kernelsize=np.int(np.floor(halfk))
        #        wind_u=wind_u[:,halfk:-halfk,halfk:-halfk]
        #        wind_v=wind_v[:,halfk:-halfk,halfk:-halfk]
        wind_u=fast_mean.fast_smooth(wind_u,halfk)[:,halfk:-halfk,halfk:-halfk]
        wind_v=fast_mean.fast_smooth(wind_v,halfk)[:,halfk:-halfk,halfk:-halfk]
        N=qc.shape
        #        print(wind_u.max(),wind_v.max())
        if self._sfc_files!=None:
            sfc=self._next_surface(self.curpos)
            sfc.pblh+=hgt[0,:,:]
            # convert the pblh into a vertical index into the 3D arrays
            pblindex=np.argmin(np.abs(sfc.pblh[np.newaxis,:,:]-hgt),axis=0)
            pblindex[pblindex<2]=2
            pblindex[pblindex>=17]=16
            sfc.pblh=pblindex+1
        else:
            sfc=None
        
        d.close()
        self.time_inc()
        return Bunch(p=pressure,th=potential_temperature,
                     sh=specific_humidity, hgt=hgt,qc=qc,
                     qv=specific_humidity/(1-specific_humidity), 
                     u=wind_u, v=wind_v,w=np.zeros(N2,dtype=np.float32,order="F"),
                     date=datestr,sfc=sfc,
                     qr=qr,qi=qi,qs=qs,qg=qg,ni=ni,nr=nr
                     )

    def time_inc(self):
        self.curpos+=1
        d=swim_io.Dataset(self._filenames[self._curfile], 'r')
        npos=d.variables[self.vvar].shape[0]-1
        if self.curpos>=npos:
            self._curfile+=1
            self.curpos=0
        d.close()
        
        
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
    N=hitopo.shape
    # just adjust all pressures by dz (difference between coarse topo and high res topo)
    dz=hitopo-lowtopo
    if len(dz.shape)==2:
        weather.p=convert_p(weather.p,weather.hgt,dz[np.newaxis,:,:].repeat(weather.p.shape[0],axis=0))
    else:
        weather.p=convert_p(weather.p,weather.hgt,dz)


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
    """Copies internal weather information from the old hi res atmosphere to the new one
        while maintaining the updated boundary conditions of the new atmosphere"""
    for k in newatm.keys():
        if (k!="hgt") and (k!="date") and (k!="sfc") and (k!="p") and (k!="u") and (k!="v") and (k!="w"):
            #newatm[k][:,1:-1,1:-1]=oldatm[k][:,1:-1,1:-1]
            tmp=np.where(np.isfinite(oldatm[k][:,1:-1,1:-1]))
            newatm[k][:,1:-1,1:-1][tmp]=oldatm[k][:,1:-1,1:-1][tmp]
            # tmp=np.where(~np.isfinite(oldatm[k][:,1:-1,1:-1]))
            # if len(tmp[0])>0:
            #     print("BROKEN:"+str(oldatm.date))
            #     oldatm[k][:,1:-1,1:-1][tmp]=newatm[k][:,1:-1,1:-1][tmp]

# Force vertical convergence/divergence to balance horizontal convergence/divergence
def adjust_winds(u,v,w,prestaggered=False):
    if not prestaggered:
        ustag=(u[:,:,1:]+u[:,:,:-1])/2
        vstag=(v[:,1:,:]+v[:,:-1,:])/2
    else:
        ustag=u
        vstag=v
    divergance=(ustag[:,1:-1,1:]-ustag[:,1:-1,:-1])+(vstag[:,1:,1:-1]-vstag[:,:-1,1:-1])
    w[0,1:-1,1:-1]=-divergance[0,:,:]
    for i in range(w.shape[0]-1):
        w[i+1,1:-1,1:-1]=w[i,1:-1,1:-1]-divergance[i+1,:,:]
    return(ustag,vstag,w)
    
def rotate_winds(u,v,hgt,dx=4000.0,rotation=None):
    if rotation==None:
        urot=np.cos(np.arctan(np.abs(hgt[:,:,1:]-hgt[:,:,:-1])/dx))
        vrot=np.cos(np.arctan(np.abs(hgt[:,1:,:]-hgt[:,:-1,:])/dx))
        rotation=[urot,vrot]
    else:
        urot,vrot=rotation
    u/=urot
    v/=vrot
    return rotation

def calc_ndsq(weather,base):
    R  = 287.0
    Rv = 461.0
    cp = 1004.0
    L   = 2.5e6
    g  = 9.81
    ratio = 18.015/28.964

    t0 = 273.15
    p0=weather.p[0,...]
    pii=1.0/((100000.0/weather.p)**(R/cp))
    T2m=weather.th[0,...]*pii[0,...]
    
    es = 611.21*np.exp(17.502*(T2m-t0)/(T2m-32.19))
    qs0 = ratio * es/(p0-es)

    cap_gamma = -(g * (1.+(L*qs0)/(R*T2m)) / (cp + (L**2 * qs0*ratio) / (R*T2m**2)))
    env_gamma = np.mean(np.diff(weather.th*pii,axis=0)/np.diff(base.hgt3d,axis=0),axis=0)
    
    dry_gamma=np.mean(env_gamma-cap_gamma)
    print(T2m.mean(),env_gamma.mean(),cap_gamma.mean(),dry_gamma)
    ndsq=(g/T2m.mean())*(dry_gamma)
    ndsq=max(min(1e-5,ndsq),1e-10)
    return ndsq
    

def simul_next(base,wrfwinds,q,r_matrix,Fzs,padx,pady):
    weather=base.next()
    topo_adjust_weather(base.hgt3d, weather.hgt, weather)

    if use_linear_winds:
        (junku,junkv,weather.w)=adjust_winds(weather.u,weather.v,weather.w)
        # calculate the (squared) brunt vaisala frequency, better off with the dry value (subtract cap_gamma)
        ndsq=calc_ndsq(weather,base)
        print("Ndsq=",ndsq)
        # print("raw winds Umax=",weather.u.max())
        # print("raw winds Vmax=",weather.v.max())
        # print("raw winds Wmax=",weather.w.max())
        lt_winds.update_winds(base.hgt3d,Fzs,weather.u,weather.v,weather.w,dx=wrfres,
                              Ndsq=ndsq,r_matrix=r_matrix,padx=padx,pady=pady,rotation=False)
        print("LT winds Umax=",weather.u.max())
        print("LT winds Vmax=",weather.v.max())
        print("LT winds Wmax=",weather.w.max())
        # weather.u=fast_mean.fast_smooth(weather.u,5)
        # weather.v=fast_mean.fast_smooth(weather.v,5)
        # weather.u=(weather.u[:,:,1:]+weather.u[:,:,:-1])/2
        # weather.v=(weather.v[:,1:,:]+weather.v[:,:-1,:])/2
        # if r_matrix!=None:
        #     rotate_winds(weather.u,weather.v,base.hgt3d,dx=wrfres,rotation=r_matrix)
        # (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)
        # print("adjusted LT winds Umax=",weather.u.max())
        # print("adjusted LT winds Vmax=",weather.v.max())
        # print("adjusted LT winds Wmax=",weather.w.max())
    elif use_wrf_winds:
        wind=wrfwinds.next()
        wind.u=fast_mean.fast_smooth(wind.u,2)
        wind.v=fast_mean.fast_smooth(wind.v,2)
        # wind.u=(wind.u[:,:,1:]+wind.u[:,:,:-1])/2
        # wind.v=(wind.v[:,1:,:]+wind.v[:,:-1,:])/2
        # if r_matrix !=None:
        #     r_matrix=rotate_winds(wind.u,wind.v,base.hgt3d,dx=wrfres,rotation=r_matrix)
        weather.u=wind.u
        weather.v=wind.v
        # (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)
        weather.p=wind.p
        weather.th=wind.th
        weather.qv=wind.qv
    # else:
    #     weather.u=(weather.u[:,:,1:]+weather.u[:,:,:-1])/2
    #     weather.v=(weather.v[:,1:,:]+weather.v[:,:-1,:])/2

    weather.u=(weather.u[:,:,1:]+weather.u[:,:,:-1])/2
    weather.v=(weather.v[:,1:,:]+weather.v[:,:-1,:])/2
    if r_matrix !=None:
        r_matrix=rotate_winds(weather.u,weather.v,base.hgt3d,dx=wrfres,rotation=r_matrix)
    (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)

    # (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w)
    
    q.put(weather)

def parallel_output(filename,Nx,Ny,Nz,varname,data,curdate):
    outputcurdate=''.join(''.join(''.join(curdate.split('-')).split('_')).split(':'))

    if Nz>0:
        ncout=swim_io.NC_writer(filename,Nx,Ny,Nz,var=varname)
    else:
        ncout=swim_io.NC_writer(filename,Nx,Ny,var=varname)
    ncout.appendToVar(data, date=outputcurdate)
    ncout.close()
    

def main(): # (file_search="nc3d/merged*.nc",topofile='/d2/gutmann/usbr/narr_data/baseline_info/2km_wrf_input_d01'):
    oldfiles=glob.glob(outputdir+"*.nc")
    for old in oldfiles:
        os.remove(old)
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

    # wrfwinds=WRF_Reader(wind_files,bilin=False,nn=True, windonly=True)
    # wind=wrfwinds.next()

    gkern=gauss_kern(40)
    wrf_base=WRF_Reader(file_search,sfc=sfc_file_search,bilin=True,nn=False)
    wrfwinds=None
    
    timestep=1.0*60.0*60.0 #3hrs
    oldweather=None
    print("Off and running...")
    weather=wrf_base.next()
    print("Mean PBL model level= "+str(weather.sfc.pblh.mean()))
    topo_adjust_weather(wrf_base.hgt3d, weather.hgt, weather)
    r_matrix=None
    if use_linear_winds:
        # for use with the linear theory winds
        (lt_topo,pady,padx)=topo_preprocess(topo)
        Fzs=fft.fftshift(fft.fft2(lt_topo))/((Nx+padx*2)*(Ny+pady*2))
        (junku,junkv,weather.w)=adjust_winds(weather.u,weather.v,weather.w)
        Ndsq=calc_ndsq(weather,wrf_base)
        r_matrix=lt_winds.update_winds(wrf_base.hgt3d,Fzs,weather.u,weather.v,weather.w,dx=wrfres,
                                      Ndsq=Ndsq,r_matrix=None,padx=padx,pady=pady,rotation=False)
        # (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,
        #                                 prestaggered=(weather.v.shape[1]!=weather.u.shape[1]))
        # weather.u=fast_mean.fast_smooth(weather.u,5)
        # weather.v=fast_mean.fast_smooth(weather.v,5)
        weather.u=(weather.u[:,:,1:]+weather.u[:,:,:-1])/2
        weather.v=(weather.v[:,1:,:]+weather.v[:,:-1,:])/2
        
        # r_matrix=rotate_winds(weather.u,weather.v,wrf_base.hgt3d,dx=wrfres,rotation=None)
        (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)
        
    elif use_wrf_winds:
        wrfwinds=WRF_Reader(wind_files,bilin=False,nn=True, windonly=True)
        wrf_base.hgt3d=wrfwinds.hgt3d #if we are using winds from wrf, we need to use the 3d domain from wrf too
        # to use high resolution WRF winds as a "TRUTH" run
        wind=wrfwinds.next()
        wind.u=(wind.u[:,:,1:]+wind.u[:,:,:-1])/2
        wind.v=(wind.v[:,1:,:]+wind.v[:,:-1,:])/2
        # r_matrix=rotate_winds(wind.u,wind.v,wrf_base.hgt3d,dx=wrfres,rotation=None)
        r_matrix=None
        weather.u=wind.u
        weather.v=wind.v
        weather.u=fast_mean.fast_smooth(weather.u,2)
        weather.v=fast_mean.fast_smooth(weather.v,2)
        weather.p=wind.p
        weather.th=wind.th
        weather.qv=wind.qv
        (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)
        padx=None
        pady=None
        Fzs=None
    else:
        weather.u=(weather.u[:,:,1:]+weather.u[:,:,:-1])/2
        weather.v=(weather.v[:,1:,:]+weather.v[:,:-1,:])/2
        # r_matrix=rotate_winds(weather.u,weather.v,wrf_base.hgt3d,dx=wrfres,rotation=None)
        r_matrix=None
        (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)
        padx=None
        pady=None
        Fzs=None
    (Nz,Ny,Nx)=weather.qv.shape
    oldt=weather.th.copy()
    oldweather=weather
    
    print("Reading WRF")
    q=Queue()
    p1=Process(target=simul_next,args=(wrf_base,wrfwinds,q,r_matrix,Fzs,padx,pady))
    p1.start()
    done=False
    print("Initializing microphysics...")
    swim_lib.swim_step.init(physics) #calls thompson_init internally if physics==1

    print(oldweather.date)
    processPool=None
    pIOlist=[None for i in range(11)]
    t1=t2=t3=t4=time.time()
    while not done:
        t1=time.time()
        weather=q.get()
        print("Mean PBL model level= "+str(weather.sfc.pblh.mean()))
        
        # because of the way the next time step is read now, the state of this object does not get updated
        # so update the time step for it here
        wrf_base.time_inc()
        if use_wrf_winds:
            wrfwinds.time_inc()
        p1.join()
        p1=Process(target=simul_next,args=(wrf_base,wrfwinds,q,r_matrix,Fzs,padx,pady))
        p1.start()

        update_weather(weather,oldweather)
        t2=time.time()
        print(" ")
        print("-------------------------------------------------")
        print(weather.date)
        print("-------------------------------------------------")
        P=swim.swim2d(topo,weather,oldweather,swim_lib,processPool,physics,timestep=timestep,dx=wrfres)
            
        t4=time.time()
        print("Finished Timestep")
        print("Total Time:"+str(t4-t3))
        print("     input :"+str(t2-t1))
        print("     phyics:"+str(t4-t2))
        print("     output:"+str(t1-t3))
        t3=time.time()
        curdate=weather.date
        
        outputnames=['swim_p_','swim_t_','swim_qv_','swim_qc_','swim_pres_','swim_qi_','swim_qs_','swim_qr_',"swim_u_","swim_v_","swim_w_"]
        varnames=['precip','temp','qv','qc','pressure','qi','qs','qr',"u","v","w"]
        outputvars=[P,weather.th,weather.qv,weather.qc,weather.p,weather.qi,weather.qs,weather.qr,weather.u,weather.v,weather.w]
        for i in range(len(varnames)):
            if i>=len(pIOlist):
                pIOlist.append(None)
            if pIOlist[i]!=None:
                pIOlist[i].join()
            if len(outputvars[i].shape)==2:
                curz=-1
                cury,curx=outputvars[i].shape
            else:
                curz,cury,curx=outputvars[i].shape
            pIOlist[i]=Process(target=parallel_output,args=(outputdir+outputnames[i]+str(curdate),
                                                            curx,cury,curz,
                                                            varnames[i],outputvars[i],curdate))
            pIOlist[i].start()
        
        oldweather=weather
    wrf_base.close()




if __name__ == '__main__':
    main()
