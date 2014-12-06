import datetime,os
import argparse
import glob

import numpy as np

from models import access,ccsm,fgoals,gfdl,ipsl
from bunch import Bunch

import io_routines as io

version="1.0"

GCM_NAMES=dict(
    access="ACCESS1-3",
    bcc="bcc-csm1-1-m",
    bnu="BNU-ESM",
    canesm="CanESM2",
    ccsm="CCSM4",
    cnrm_cm5="CNRM-CM5",
    fgoals="FGOALS-g2",
    gfdl_cm3="GFDL-CM3",
    gfdl_esm="GFDL-ESM2M",
    giss_e2h="GISS-E2-H",
    hadgem_es="HadGEM2-ES",
    ipsl_mr="IPSL-CM5A-MR",
    miroc5="MIROC5",
    miroc_esm="MIROC-ESM",
    mk3="CSIRO-Mk3-6-0",
    mri_cgcm3="MRI-CGCM3",
    noresm="NorESM1-M"
)

global_vert_coords=dict(
                        access=access.vcoord,
                        bcc = ccsm.vcoord,
                        bnu = ccsm.vcoord,
                        canesm = ipsl.vcoord,
                        ccsm = ccsm.vcoord,
                        cnrm_cm5 = ccsm.vcoord,
                        fgoals = fgoals.vcoord,
                        gfdl_cm3 = gfdl.vcoord,
                        gfdl_esm = gfdl.vcoord,
                        ipsl_mr = ipsl.vcoord,
                        miroc5 = ccsm.vcoord,
                        miroc_esm = ccsm.vcoord,
                        mk3 = ccsm.vcoord,
                        mri_cgcm3 = ccsm.vcoord,
                        noresm = ccsm.vcoord,
                        giss_e2h = ccsm.vcoord
                        )


def set_bounds(info):
    atm_file=info.atmdir+info.atmfile
    print(atm_file)
    cmip_file= atm_file.replace("_Y_",str(info.start_year))        \
                    .replace("_VAR_","hus")                        \
                    .replace("_ENS_",info.ensemble)                \
                    .replace("_EXP_",info.experiment)              \
                    .replace("_GCM_",info.gcm_name)
    print(cmip_file)
    cmip_file=glob.glob(cmip_file)[0]
    print(cmip_file)
    varlist=["lat","lon"]
    
    lat=io.read_nc(cmip_file,varlist[0]).data
    lon=io.read_nc(cmip_file,varlist[1]).data
    if lon.max()>180:
        lon[lon>180]=lon[lon>180]-360
    if len(lat.shape)>1 or len(lon.shape)>1:
        print(lat.shape,len(lat.shape))
        print(lon.shape,len(lon.shape))
        raise ValueError("config.py requires lat,lon to be 1D arrays")
    lonbounds=np.where((lon>=info.lon[0])&(lon<=info.lon[1]))[0]
    info.xmin=max(lonbounds[0]-1,0)
    info.xmax=min(lonbounds[-1]+1,len(lon))
    info.ymin=max(np.where(lat>=info.lat[0])[0][0]-1,0)
    info.ymax=min(np.where(lat<=info.lat[1])[0][-1]+1,len(lat))
    
    lon,lat=np.meshgrid(lon[info.xmin:info.xmax],lat[info.ymin:info.ymax])
    info.lat_data=lat
    info.lon_data=lon
    
def make_timelist(info,hrs=6.0):
    dt=datetime.timedelta(hrs/24.0)
    info.ntimes=np.int(np.round((info.end_date-info.start_date).total_seconds()/60./60./hrs))
    info.times=[info.start_date+dt*i for i in range(info.ntimes)]

def update_info(info):
    make_timelist(info)
    set_bounds(info)
    

def parse():
    parser= argparse.ArgumentParser(description='Convert cmip files to ICAR input forcing files')
    parser.add_argument('model',     nargs="?",action='store',help="GCM (ccsm,giss,miroc5,...)",         default="ccsm")
    parser.add_argument('experiment',nargs="?",action='store',help="Experiment (historical,rcp85)",      default="historical")
    parser.add_argument('start_year',nargs="?",action='store',help="Specify starting year (yyyy)",       default="1990")
    parser.add_argument('nyears',    nargs="?",action='store',help="Number of years to run (10)",        default="10")
    parser.add_argument('ensemble',  nargs="?",action='store',help="Ensemble member (r1i1p1)",           default="r1i1p1")
    parser.add_argument('lat_n',     nargs="?",action='store',help="northern latitude boundary",         default="60")
    parser.add_argument('lat_s',     nargs="?",action='store',help="southern latitude boundary",         default="20")
    parser.add_argument('lon_e',     nargs="?",action='store',help="eastern longitude boundary",         default="-60")
    parser.add_argument('lon_w',     nargs="?",action='store',help="western longitude boundary",         default="-130")
    parser.add_argument('output',    nargs="?",action='store',help="output file prefix (GCM replaced)",  default="GCM_")
    parser.add_argument('dir',       nargs="?",action='store',help="cmip file location",                 default="3D_6hrly_raw/")
    parser.add_argument('atmfile',   nargs="?",action='store',help="cmip atmospheric files",
                                     default="_VAR__6hrLev__GCM___EXP___ENS___Y_*.nc")

    parser.add_argument('-v', '--version',action='version',
            version='cmip2icar v'+version)
    parser.add_argument ('--verbose', action='store_true',
            default=False, help='verbose output', dest='verbose')
    parser.add_argument ('--print', action='store_true',
            default=False, help='print available models', dest='print_models')
    args = parser.parse_args()
    
    if args.print_models:
        for k in global_vert_coords.keys():
            print(k)
        sys.exit()
    
    start_year=int(args.start_year)
    nyears=int(args.nyears)
    start_date=datetime.datetime(start_year,1,1,0,0,0)
    end_date=datetime.datetime(start_year+nyears,1,1,0,0,0)
    
    info=Bunch(lat=[float(args.lat_s),float(args.lat_n)],
               lon=[float(args.lon_w),float(args.lon_e)],
               start_date=start_date,  end_date=end_date,
               start_year=int(args.start_year),experiment=args.experiment,
               ensemble=args.ensemble,
               model=args.model,
               atmdir=args.dir,
               atmfile=args.atmfile,
               gcm_name=GCM_NAMES[args.model],
               read_pressure=global_vert_coords[args.model],
               orog_file="orography.nc",
               output_file=args.output.replace("GCM",args.model)+args.ensemble+"_"+args.experiment+"_",
               version=version)
    
    return info
