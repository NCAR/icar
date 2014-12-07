import numpy as np
import mygis
from bunch import Bunch

def write_file(date,info,cmip):
    """writes cmip input data to a netcdf file"""
    
    filename=info.output_file+str(date).replace(" ","_")
    print("Outputting: "+filename)
    dims=("time","level","lat","lon")
    dims_3d=("time","lat","lon")
    dims_2d=("lat","lon")
    
    extra_vars=[]
    # 3D variables (+time)
    # cloud,ice,qv,u,v,t,p
    # z
    # 2D variables (constant in time)
    # hgt, latitude, longitude
    
    atts=Bunch(long_name="Cloud liquid water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="cloud",data=cmip["cloud"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Cloud ice water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="ice",data=cmip["ice"],dims=dims,dtype="f",attributes=atts))

    # used as primary variable in io.write
    # atts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    # extra_vars.append(Bunch(name="qv",data=cmip["qv"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="U (E/W) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="u",data=cmip["u"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="V (N/S) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="v",data=cmip["v"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Potential Temperature",units="kg kg**-1")
    extra_vars.append(Bunch(name="theta",data=cmip["t"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Pressure",units="Pa")
    extra_vars.append(Bunch(name="p",data=cmip["p"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Layer thicknesses",units="m")
    extra_vars.append(Bunch(name="dz",data=cmip["dz"].astype("f"),dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Layer height",units="m")
    extra_vars.append(Bunch(name="z",data=cmip["z"].astype("f"),dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Topographic Height",units="m")
    # print(cmip["hgt"].shape,dims_2d,cmip.qv.shape)
    extra_vars.append(Bunch(name="hgt",data=cmip["hgt"],dims=dims_2d,dtype="f",attributes=atts))

    atts=Bunch(long_name="latitude",units="degrees")
    extra_vars.append(Bunch(name="lat",data=info.lat_data,dims=dims_2d,dtype="f",attributes=atts))
    
    atts=Bunch(long_name="longitude",units="degrees")
    extra_vars.append(Bunch(name="lon",data=info.lon_data,dims=dims_2d,dtype="f",attributes=atts))

    atts=Bunch(long_name="xland",units="")
    extra_vars.append(Bunch(name="xland",data=cmip["land"],dims=dims_2d,dtype="f",attributes=atts))


    qvatts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    
    # write to output file
    mygis.write(filename=filename,varname="qv",data=cmip.qv,dims=dims, attributes=qvatts,dtype="f",
                  extravars=extra_vars)#,history=" Produced by cmip2icar v."+info.version)
