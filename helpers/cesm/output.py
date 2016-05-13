import numpy as np
import mygis
from bunch import Bunch

def write_file(date,info,cesm):
    """writes cesm input data to a netcdf file"""
    
    filename=info.output_file+str(date).replace(" ","_")
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
    
    # used as primary variable in io.write
    # atts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    # extra_vars.append(Bunch(name="qv",data=cesm["qv"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Cloud liquid water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="cloud",data=cesm["cloud"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Cloud ice water content",units="kg kg**-1")
    extra_vars.append(Bunch(name="ice",data=cesm["ice"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="U (E/W) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="u",data=cesm["u"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="V (N/S) wind speed",units="m s**-1")
    extra_vars.append(Bunch(name="v",data=cesm["v"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Potential Temperature",units="kg kg**-1")
    extra_vars.append(Bunch(name="theta",data=cesm["t"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Pressure",units="Pa")
    extra_vars.append(Bunch(name="p",data=cesm["p"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Layer thicknesses",units="m")
    extra_vars.append(Bunch(name="dz",data=cesm["dz"].astype("f"),dims=(dims[1],),dtype="f",attributes=atts))

    atts=Bunch(long_name="Layer height",units="m")
    extra_vars.append(Bunch(name="z",data=cesm["z"].astype("f"),dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Topographic Height",units="m")
    print(cesm["hgt"].shape,dims_2d,cesm.qv.shape)
    extra_vars.append(Bunch(name="hgt",data=cesm["hgt"],dims=dims_2d,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Shortwave Radiation (positive down)",units="W m**-2")
    extra_vars.append(Bunch(name="swdown",data=cesm["sw"],dims=dims_3d,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Longwave Radiation (positive down)",units="W m**-2")
    extra_vars.append(Bunch(name="lwdown",data=cesm["lw"],dims=dims_3d,dtype="f",attributes=atts))

    atts=Bunch(long_name="Skin Temperature",units="K")
    extra_vars.append(Bunch(name="tskin",data=cesm["ts"],dims=dims_3d,dtype="f",attributes=atts))

    atts=Bunch(long_name="latitude",units="degrees")
    extra_vars.append(Bunch(name="lat",data=info.lat_data,dims=dims_2d,dtype="f",attributes=atts))
    
    atts=Bunch(long_name="longitude",units="degrees")
    extra_vars.append(Bunch(name="lon",data=info.lon_data,dims=dims_2d,dtype="f",attributes=atts))

    atts=Bunch(long_name="xland",units="")
    extra_vars.append(Bunch(name="xland",data=cesm["land"],dims=dims_2d,dtype="f",attributes=atts))


    for k in cesm.keys():
        print(k,cesm[k].shape)

    qvatts=Bunch(long_name="Specific Humidity",units="kg kg**-1")
    
    print(" ")
    print(" ")
    print("Writing:"+filename)
    print(" ")
    print(" ")
    # write to output file
    mygis.write(filename=filename,varname="qv",data=cesm.qv,dims=dims, attributes=qvatts,dtype="f",
                  extravars=extra_vars)#,history=" Produced by cesm2icar v."+info.version)
