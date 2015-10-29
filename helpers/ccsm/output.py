import numpy as np
import mygis
from bunch import Bunch

def write_file(date,info,ccsm):
    """writes CCSM input data to a netcdf file"""
    
    filename=str(date).replace(" ","_")
    dims=("time","level","lat","lon")
    dims_3d=("time","lat","lon")
    dims_2d=("lat","lon")
    
    extra_vars=[]
    # 3D variables (+time)
    # cloud,ice,qv,u,v,t,p
    # 3D variables (constant in time)
    # z
    # 2D variables (+time)
    # latent_heat,PBL_height,sensible_heat
    # 2D variables (constant in time)
    # hgt, latitude, longitude
    
    atts=Bunch(long_name="Cloud liquid water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="cloud",data=ccsm["cloud"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Cloud ice water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="ice",data=ccsm["ice"],dims=dims,dtype="f",attributes=atts))

    # used as primary variable in io.write
    # atts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    # extra_vars.append(Bunch(name="qv",data=ccsm["qv"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="U (E/W) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="u",data=ccsm["u"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="V (N/S) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="v",data=ccsm["v"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Potential Temperature",units="kg kg**-1")
    extra_vars.append(Bunch(name="theta",data=ccsm["t"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Pressure",units="Pa")
    extra_vars.append(Bunch(name="p",data=ccsm["p"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Layer thicknesses",units="m")
    extra_vars.append(Bunch(name="dz",data=ccsm["dz"].astype("f"),dims=(dims[1],),dtype="f",attributes=atts))

    atts=Bunch(long_name="Topographic Height",units="m")
    extra_vars.append(Bunch(name="hgt",data=ccsm["hgt"],dims=dims_2d,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Latent Heat flux (positive up)",units="W m**-2")
    extra_vars.append(Bunch(name="latent_heat",data=ccsm["latent_heat"],dims=dims_3d,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Sensible Heat flux (positive up)",units="W m**-2")
    extra_vars.append(Bunch(name="sensible_heat",data=ccsm["sensible_heat"],dims=dims_3d,dtype="f",attributes=atts))

    atts=Bunch(long_name="latitude",units="degrees")
    extra_vars.append(Bunch(name="lat",data=info.lat_data,dims=dims_2d,dtype="f",attributes=atts))
    
    atts=Bunch(long_name="longitude",units="degrees")
    extra_vars.append(Bunch(name="lon",data=info.lon_data,dims=dims_2d,dtype="f",attributes=atts))


    qvatts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    
    # write to output file
    mygis.write(filename=filename,varname="qv",data=ccsm.qv,attributes=qvatts,dtype="f",
                  extravars=extra_vars)#,history=" Produced by ccsm2icar v."+info.version)
