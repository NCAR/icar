import datetime,os
import argparse
import glob

import numpy as np

from bunch import Bunch

import io_routines as io

version="1.0"

def set_bounds(info):
    atm_file=info.atmdir+info.atmfile
    ccsm_file=atm_file.replace("_Y_","2006").replace("_M_","01").replace("_D_","01").replace("_VAR_","hus")
    ccsm_file=glob.glob(ccsm_file)[0]
    varlist=["lat","lon"]

    lat=io.read_nc(ccsm_file,varlist[0]).data
    lon=io.read_nc(ccsm_file,varlist[1]).data-360

    info.xmin=np.where(lon>=info.lon[0])[0][0]
    info.xmax=np.where(lon<=info.lon[1])[0][-1]+1
    info.ymin=np.where(lat>=info.lat[0])[0][0]
    info.ymax=np.where(lat<=info.lat[1])[0][-1]+1

    lon,lat=np.meshgrid(lon[info.xmin:info.xmax],lat[info.ymin:info.ymax])
    info.lat_data=lat
    info.lon_data=lon

def make_timelist(info,hrs=6.0):
    dt=datetime.timedelta(hrs/24)
    info.ntimes=np.int64(np.round((info.end_date-info.start_date).total_seconds()/60./60./hrs))
    info.times=[info.start_date+dt*i for i in range(info.ntimes)]

def update_info(info):
    make_timelist(info)
    set_bounds(info)


def parse():
    parser= argparse.ArgumentParser(description='Convert CCSM files to ICAR input forcing files')
    parser.add_argument('start_date',nargs="?",action='store',help="Specify starting date (yyyy-mm-dd)", default="2026-01-01")
    parser.add_argument('end_date',  nargs="?",action='store',help="Specify end date (yyyy-mm-dd)",      default="2099-12-31")
    parser.add_argument('lat_n',     nargs="?",action='store',help="northern latitude boundary",         default="60")
    parser.add_argument('lat_s',     nargs="?",action='store',help="southern latitude boundary",         default="20")
    parser.add_argument('lon_e',     nargs="?",action='store',help="eastern longitude boundary",         default="-60")
    parser.add_argument('lon_w',     nargs="?",action='store',help="western longitude boundary",         default="-130")
    parser.add_argument('dir',       nargs="?",action='store',help="CCSM file location",                 default="/d5/gutmann/ccsm/")
    parser.add_argument('atmdir',    nargs="?",action='store',help="CCSM atmospheric data file location",default="ccsm_atm/")
    parser.add_argument('sfcdir',    nargs="?",action='store',help="CCSM surface data file location",    default="ccsm_sfc/")
    parser.add_argument('atmfile',   nargs="?",action='store',help="CCSM atmospheric files",             default="_VAR__6hrLev_CCSM4_rcp85_r6i1p1__Y__M__D_00-*.nc")
    parser.add_argument('sfcfile',   nargs="?",action='store',help="CCSM surface files",                 default="_VAR__3hr_CCSM4_rcp85_r6i1p1__Y__M__D_0000-*.nc")

    parser.add_argument('-v', '--version',action='version',
            version='CCSM2ICAR v'+version)
    parser.add_argument ('--verbose', action='store_true',
            default=False, help='verbose output', dest='verbose')
    args = parser.parse_args()

    date0=args.start_date.split("-")
    start_date=datetime.datetime(int(date0[0]),int(date0[1]),int(date0[2]))

    date0=args.end_date.split("-")
    end_date=datetime.datetime(int(date0[0]),int(date0[1]),int(date0[2]))

    info=Bunch(lat=[float(args.lat_s),float(args.lat_n)],
               lon=[float(args.lon_w),float(args.lon_e)],
               start_date=start_date,         end_date=end_date,
               atmdir=args.dir+args.atmdir,   sfcdir=args.dir+args.sfcdir,
               atmfile=args.atmfile,          sfcfile=args.sfcfile,
               version=version)

    return info
