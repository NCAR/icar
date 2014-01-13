import numpy as np
import swim_io
from bunch import Bunch

def write_file(date,info,erai):
    """writes ERAi input data to a netcdf file"""
    
    filename=str(date).replace(" ","_")
    dims=("level","lat","lon")
    
    extra_vars=[]
    # 3D variables
    # cloud,ice,qv,u,v,t,p
    # 2D variables
    # hgt,latent_heat,PBL_height,sensible_heat,sfc_hgt (sfc_hgt not used currently should be ~the same as hgt)
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

    atts=Bunch(long_name="Potential Temperature",units="kg kg**-1")
    extra_vars.append(Bunch(name="theta",data=erai["t"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Pressure",units="Pa")
    extra_vars.append(Bunch(name="p",data=erai["p"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Atmospheric Elevation",units="m")
    extra_vars.append(Bunch(name="z",data=erai["z"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Topographic Height",units="m")
    extra_vars.append(Bunch(name="hgt",data=erai["hgt"],dims=dims[1:],dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Latent Heat flux (positive up)",units="W m**-2")
    extra_vars.append(Bunch(name="latent_heat",data=erai["latent_heat"],dims=dims[1:],dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Sensible Heat flux (positive up)",units="W m**-2")
    extra_vars.append(Bunch(name="sensible_heat",data=erai["sensible_heat"],dims=dims[1:],dtype="f",attributes=atts))

    atts=Bunch(long_name="Planetary Boundary Layer Height",units="m")
    extra_vars.append(Bunch(name="PBL_height",data=erai["PBL_height"],dims=dims[1:],dtype="f",attributes=atts))

    atts=Bunch(long_name="latitude",units="degrees")
    extra_vars.append(Bunch(name="lat",data=info.lat_data,dims=dims[1:],dtype="f",attributes=atts))
    
    atts=Bunch(long_name="longitude",units="degrees")
    extra_vars.append(Bunch(name="lon",data=info.lon_data,dims=dims[1:],dtype="f",attributes=atts))


    qvatts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    
    # write to output file
    swim_io.write(filename=filename,varname="qv",data=erai.qv,attributes=qvatts,dtype="f",
                  extravars=extra_vars,history=" Produced by erai2swm v."+info.version)
