import datetime
import sys

import numpy as np

import mygis
from bunch import Bunch

def write_file(date,info,erai):
    """writes ERAi input data to a netcdf file"""

    filename=str(date).replace(" ","_")
    dims    = ("time", "level","lat","lon")
    dims2dt = ("time", "lat","lon")

    extra_vars=[]
    # 3D variables
    # cloud,ice,qv,u,v,t,p, z
    # 2D variables
    # hgt,latent_heat,PBL_height,sensible_heat,sfc_hgt (sfc_hgt not used currently should be ~the same as hgt)
    # 1D variables / coordinates
    # lat, lon
    atts=Bunch(long_name="Cloud liquid water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="cloud",data=erai["cloud"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Cloud ice water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="ice",data=erai["ice"],dims=dims,dtype="f",attributes=atts))

    # used as primary variable in io.write
    # atts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    # extra_vars.append(Bunch(name="qv",data=erai["qv"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="U (E/W) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="u",data=erai["u"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="V (N/S) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="v",data=erai["v"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Potential Temperature",units="K")
    extra_vars.append(Bunch(name="theta",data=erai["t"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Pressure",units="Pa")
    extra_vars.append(Bunch(name="p",data=erai["p"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Atmospheric Elevation",units="m",positive="up")
    extra_vars.append(Bunch(name="z",data=erai["z"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Topographic Height",units="m")
    extra_vars.append(Bunch(name="hgt",data=erai["hgt"],dims=dims[2:],dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface solar radiation (downwards)",units="W m**-2")
    extra_vars.append(Bunch(name="swdown",data=erai["sw"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface longwave radiation (downwards)",units="W m**-2")
    extra_vars.append(Bunch(name="lwdown",data=erai["lw"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Latent Heat flux (positive up)",units="W m**-2")
    extra_vars.append(Bunch(name="latent_heat",data=erai["latent_heat"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Sensible Heat flux (positive up)",units="W m**-2")
    extra_vars.append(Bunch(name="sensible_heat",data=erai["sensible_heat"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Planetary Boundary Layer Height",units="m")
    extra_vars.append(Bunch(name="PBL_height",data=erai["PBL_height"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Land fraction",units="")
    extra_vars.append(Bunch(name="landfraction",data=erai["landmask"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Skin Temperature",units="K")
    extra_vars.append(Bunch(name="tskin",data=erai["tskin"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Sea Surface Temperature",units="K")
    extra_vars.append(Bunch(name="sst",data=erai["sst"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Convective precipitation",units="mm")
    extra_vars.append(Bunch(name="cp",data=erai["cp"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="latitude",units="degrees_north")
    extra_vars.append(Bunch(name="lat",data=info.lat_data[:,0],dims=("lat",),dtype="f",attributes=atts))

    atts=Bunch(long_name="longitude",units="degrees_east")
    extra_vars.append(Bunch(name="lon",data=info.lon_data[0,:],dims=("lon",),dtype="f",attributes=atts))

    time_since_1900 = date - datetime.datetime(1900,1,1,0,0,0)
    time = time_since_1900.days + np.float64(time_since_1900.seconds/86400.0)
    atts=Bunch(long_name="time",units="days since 1900-01-01", calendar="gregorian")
    extra_vars.append(Bunch(name="time",data=time,dims=(dims[0],),dtype="d",attributes=atts))


    qvatts=Bunch(long_name="Specific Humidity",units="kg kg**-1")

    # write to output file
    mygis.write(filename=filename,varname="qv",data=erai.qv,attributes=qvatts,dtype="f",dims=dims,
                  extravars=extra_vars,history=" Produced by erai2icar v."+info.version+" "+" ".join(sys.argv))
