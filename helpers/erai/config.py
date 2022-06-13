import datetime,os
import argparse

import numpy as np

from bunch import Bunch
import mygis

import io_routines as io

version="1.1"

def set_bounds(info):
    atm_file=info.atmdir+info.atmfile
    erai_file=atm_file.replace("_Y_","2000").replace("_M_","01").replace("_D_","01").replace("_h_","00")
    varlist=["g4_lat_0","g4_lon_1","Z_GDS4_HYBL","T_GDS4_HYBL","Q_GDS4_HYBL","LNSP_GDS4_HYBL","CLWC_GDS4_HYBL","CIWC_GDS4_HYBL","lv_HYBL2_a","lv_HYBL2_b","P0"]
    output_dir=info.nc_file_dir
    try:
        os.mkdir(output_dir)
    except:
        pass

    ncfile = io.grib2nc(erai_file,varlist,output_dir)
    lat = mygis.read_nc(ncfile,varlist[0]).data
    lon = mygis.read_nc(ncfile,varlist[1]).data - 360

    # print(lon, info.lon[0])
    lon[lon<-180]+=360
    info.xmin = max(0,np.argmin(np.abs(lon-info.lon[0]))-1)
    info.xmax = min(lon.size-1,np.argmin(np.abs(lon-info.lon[1]))+1)
    if (info.xmax < info.xmin):
        print("ERROR: attempting to wrap around the ERAi boundary lon="+str(lon[0])+str(lon[-1]))
        print("  Requested East lon = "+str(info.lon[0]))
        print("  Requested West lon = "+str(info.lon[1]))
        print("Requires Custom Solution!")
        raise IndexError

    # info.xmin = max(0,np.where(lon >= info.lon[0])[0][0]-1)
    # info.xmax = min(lon.size-1, np.where(lon[info.xmin:] >= info.lon[1])[0][0] + info.xmin + 1)
    #note lat is inverted from "expected"
    info.ymin=np.where(lat<=info.lat[1])[0][0]
    info.ymax=np.where(lat>=info.lat[0])[0][-1]+1
    lon,lat=np.meshgrid(lon[info.xmin:info.xmax],lat[info.ymin:info.ymax][::-1])
    info.lat_data=lat
    info.lon_data=lon

    # important to remove this file so that if this temporary date overlaps the actual time period we want to run,
    # later steps won't break when they try to use this temporary file by mistake
    os.remove(ncfile)

def make_timelist(info):
    hrs=6.0
    dt=datetime.timedelta(hrs/24)
    info.ntimes=np.int(np.round((info.end_date-info.start_date).total_seconds()/60./60./hrs))
    info.times=[info.start_date+dt*i for i in range(info.ntimes)]

def update_info(info):
    make_timelist(info)

    set_bounds(info)


def parse():
    parser= argparse.ArgumentParser(description='Convert ERAi files to ICAR input forcing files')
    parser.add_argument('start_date',nargs="?",    action='store',  help="Specify starting date (yyyy-mm-dd)", default="2000-10-01")
    parser.add_argument('end_date', nargs="?",     action='store',  help="Specify end date (yyyy-mm-dd)",      default="2000-10-02")
    parser.add_argument('lat_n',    nargs="?",     action='store',  help="northern latitude boundary",         default="60")
    parser.add_argument('lat_s',    nargs="?",     action='store',  help="southern latitude boundary",         default="20")
    parser.add_argument('lon_e',    nargs="?",     action='store',  help="eastern longitude boundary",         default="-50")
    parser.add_argument('lon_w',    nargs="?",     action='store',  help="western longitude boundary",         default="-140")
    parser.add_argument('dir',      nargs="?",     action='store',  help="ERAi file location",                 default="/glade/collections/rda/data/ds627.0/")
    parser.add_argument('atmdir',   nargs="?",     action='store',  help="ERAi atmospheric data file location",default="ei.oper.an.ml/_Y__M_/")
    parser.add_argument('sfcdir',   nargs="?",     action='store',  help="ERAi surface data file location",    default="ei.oper.fc.sfc/_Y__M_/")
    parser.add_argument('atmfile',  nargs="?",     action='store',  help="ERAi primary atmospheric file",      default="ei.oper.an.ml.regn128sc._Y__M__D__h_")
    parser.add_argument('atmuvfile',nargs="?",     action='store',  help="ERAi U/V atm file",                  default="ei.oper.an.ml.regn128uv._Y__M__D__h_")
    parser.add_argument('sfcfile',  nargs="?",     action='store',  help="ERAi surface file",                  default="ei.oper.fc.sfc.regn128sc._Y__M__D__h_")
    parser.add_argument('temp_nc_dir',nargs="?",   action='store',  help="temporary directory to store netCDF files in",default="temp_nc_dir")

    parser.add_argument('-v', '--version',action='version',
            version='ERAi2ICAR v'+version)
    parser.add_argument ('--verbose', action='store_true',
            default=False, help='verbose output', dest='verbose')
    args = parser.parse_args()

    date0=args.start_date.split("-")
    start_date=datetime.datetime(int(date0[0]),int(date0[1]),int(date0[2]))

    date0=args.end_date.split("-")
    end_date=datetime.datetime(int(date0[0]),int(date0[1]),int(date0[2]))

    if args.temp_nc_dir[-1]!="/":
        args.temp_nc_dir+="/"

    info=Bunch(lat=[float(args.lat_s),float(args.lat_n)],
               lon=[float(args.lon_w),float(args.lon_e)],
               start_date=start_date,end_date=end_date,
               atmdir=args.dir+args.atmdir,sfcdir=args.dir+args.sfcdir,
               atmfile=args.atmfile,uvfile=args.atmuvfile,sfcfile=args.sfcfile,
               nc_file_dir=args.temp_nc_dir,version=version)

    return info
