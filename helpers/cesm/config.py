import datetime,os
import argparse
import glob
import sys

import numpy as np

from bunch import Bunch

import io_routines as io

version="1.0"

def set_bounds(info):
    atm_file=info.atmdir+info.atmfile
    cesm_file=atm_file.replace("_Y_",str(info.start_year)).replace("_VAR_","Q").replace("_ENS_",info.ensemble).replace("_EXP_",info.experiment)
    try:
        cesm_file=glob.glob(cesm_file)[0]
    except:
        print("ERROR searching for: "+cesm_file)
        sys.exit(1)
    varlist=["lat","lon"]
    
    lat=io.read_nc(cesm_file,varlist[0]).data
    lon=io.read_nc(cesm_file,varlist[1]).data#-360
    
    try:
        info.xmin=np.where(lon>=info.lon[0])[0][0]
        info.xmax=np.where(lon<=info.lon[1])[0][-1]+1
        info.ymin=np.where(lat>=info.lat[0])[0][0]
        info.ymax=np.where(lat<=info.lat[1])[0][-1]+1
    except:
        print("ERROR searching for lat/lon bounds:")
        print(info.lat, info.lon)
        print(lat.min(), lat.max())
        print(lon.min(), lon.max())
        sys.exit(1)
    
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
    parser= argparse.ArgumentParser(description='Convert cesm files to ICAR input forcing files',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('experiment',nargs="?",action='store',help="Experiment (20TR,RCP85)",            default="20TR")
    parser.add_argument('start_year',nargs="?",action='store',help="Specify starting year (yyyy)",       default="1990")
    parser.add_argument('ensemble',  nargs="?",action='store',help="Ensemble member (030)",              default="030")
    parser.add_argument('lat_n',     nargs="?",action='store',help="northern latitude boundary",         default="60")
    parser.add_argument('lat_s',     nargs="?",action='store',help="southern latitude boundary",         default="20")
    parser.add_argument('lon_e',     nargs="?",action='store',help="eastern longitude boundary",         default="-60")
    parser.add_argument('lon_w',     nargs="?",action='store',help="western longitude boundary",         default="-130")
    parser.add_argument('output',    nargs="?",action='store',help="output file prefix",                 default="cesm_")
    parser.add_argument('dir',       nargs="?",action='store',help="cesm file location",                 default="3D_6hrly_raw/")
    parser.add_argument('atmfile',   nargs="?",action='store',help="cesm atmospheric files",             default="b.e11.B_EXP_C5CNBDRD.f09_g16._ENS_.cam.h2._VAR_._Y_*.nc")
    parser.add_argument('nyears',    nargs="?",action='store',help="number of years to process",         default="10")

    parser.add_argument('-v', '--version',action='version',
            version='cesm2icar v'+version)
    parser.add_argument ('--verbose', action='store_true',
            default=False, help='verbose output', dest='verbose')
    args = parser.parse_args()
    
    nyears=int(args.nyears)
    start_year=int(args.start_year)
    start_date=datetime.datetime(start_year,1,1,0,0,0)
    end_date=datetime.datetime(start_year+nyears,1,1,0,0,0)
    
    if (float(args.lon_w)<0):
        args.lon_w = 360+float(args.lon_w)
    if (float(args.lon_e)<0):
        args.lon_e = 360+float(args.lon_e)
    info=Bunch(lat=[float(args.lat_s),float(args.lat_n)],
               lon=[float(args.lon_w),float(args.lon_e)],
               start_date=start_date,  end_date=end_date,
               start_year=int(args.start_year), nyears=nyears,
               experiment=args.experiment,
               ensemble=args.ensemble,
               atmdir=args.dir,
               atmfile=args.atmfile,
               output_file=args.output+args.ensemble+"_"+args.experiment+"_",
               version=version)
    
    return info
