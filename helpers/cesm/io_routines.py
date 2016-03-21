import subprocess,os,glob
import numpy as np
import netCDF4
from bunch import Bunch
import gc,sys

g=9.8
atmvarlist=["T","Q","U","V","Z3"]
icar_atm_var=["t","qv","u","v","z"]

# from mygis, modified to work with netCDF4
def read_nc(filename,var="data",proj=None,returnNCvar=False):
    '''read a netCDF file and return the specified variable

    output is a structure :
        data:raw data as an array
        proj:string representation of the projection information
        atts:data attribute dictionary (if any)
    if (returnNCvar==True) then the netCDF4 file is note closed and the netCDF4 
        representation of the variable is returned instead of being read into 
        memory immediately.  
    '''
    d=netCDF4.Dataset(filename, mode='r',format="nc")
    outputdata=None
    if var != None:
        data=d.variables[var]
        attributes=d.variables[var].__dict__
        if returnNCvar:
            outputdata=data
        else:
            outputdata=data[:]
            # ntimes=365*4
            # if len(data.shape)>2:
                # outputdata=data[:ntimes,...]
            # else:
                # outputdata=data[:]
    outputproj=None
    if proj!=None:
        projection=d.variables[proj]
        outputproj=str(projection)
    
    
    if returnNCvar:
        return Bunch(data=outputdata,proj=outputproj,ncfile=d,atts=attributes)
    d.close()
    return Bunch(data=outputdata,proj=outputproj,atts=attributes)

def find_atm_file(time,varname,info):
    file_base= info.atmdir+info.atmfile
    file_base= file_base.replace("_VAR_",varname)
    file_base= file_base.replace("_Y_",str(info.start_year))
    file_base= file_base.replace("_EXP_",info.experiment)
    atm_file = file_base.replace("_ENS_",info.ensemble)
    
    print(atm_file)
    return glob.glob(atm_file)[0]


def load_atm(time,info,starttime,endtime):
    """Load atmospheric variable from a netcdf file"""
    
    outputdata=Bunch()

    for s,v in zip(icar_atm_var,atmvarlist):
        atmfile=find_atm_file(time,v,info)
        print(atmfile)
        sys.stdout.flush()
        nc_data=read_nc(atmfile,v,returnNCvar=True)
        outputdata[s]=nc_data.data[starttime:endtime,:,info.ymin:info.ymax,info.xmin:info.xmax]
        nc_data.ncfile.close()

    atmfile=find_atm_file(time,"PS",info)
    nc_data=read_nc(atmfile,"PS",returnNCvar=True)
    outputdata.ps=nc_data.data[starttime:endtime,info.ymin:info.ymax,info.xmin:info.xmax]
    nc_data.ncfile.close()
    del nc_data
    print(gc.collect())
    sys.stdout.flush()
    
    a=read_nc(atmfile,"hyam").data
    b=read_nc(atmfile,"hybm").data
    p0=read_nc(atmfile,"P0").data
    #p_(i,j,k)= A_k * P_0 + B_k P_s(i,j)  from http://www.cesm.ucar.edu/models/atm-cam/docs/usersguide/node25.html
    outputdata.p = a[np.newaxis,:,np.newaxis,np.newaxis]*p0+b[np.newaxis,:,np.newaxis,np.newaxis]*outputdata.ps[:,np.newaxis,:,:]
    
    outputdata.ntimes=outputdata.p.shape[0]
    
    return outputdata

def load_sfc(time, info,starttime,endtime):
    """docstring for load_sfc"""
    outputdata=Bunch()
    basefile="/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/USGS-gtopo30_0.9x1.25_remap_c051027.nc"
    outputdata.hgt=read_nc(basefile,"PHIS").data[info.ymin:info.ymax,info.xmin:info.xmax]/g

    outputdata.land=np.zeros(outputdata.hgt.shape)
    landfrac=read_nc(basefile,"LANDFRAC").data[info.ymin:info.ymax,info.xmin:info.xmax]
    outputdata.land[landfrac>=0.5]=1
    
    tsfile=find_atm_file(time, "TS", info)
    outputdata.ts=read_nc(tsfile,"TS").data[starttime:endtime,info.ymin:info.ymax,info.xmin:info.xmax]

    swfile=find_atm_file(time, "FSDS", info)
    tmp=read_nc(swfile,"FSDS",returnNCvar=True)
    print(tmp.data.shape)
    tmp.ncfile.close()
    tmp=read_nc(swfile,"FSDS").data
    print(tmp.shape, starttime, endtime)
    outputdata.sw=tmp[starttime:endtime,info.ymin:info.ymax,info.xmin:info.xmax]
    print(swfile, starttime, endtime, info.xmin,info.xmax, info.ymin, info.ymax)
    print(outputdata.sw.shape)
    print(outputdata.sw[0].max(),outputdata.sw[0].min())
    print(outputdata.sw.max(),outputdata.sw.min())

    lwfile=find_atm_file(time, "FLDS", info)
    outputdata.lw=read_nc(lwfile,"FLDS").data[starttime:endtime,info.ymin:info.ymax,info.xmin:info.xmax]
    
    return outputdata

def load_data(time,info,starttime,endtime):
    """docstring for load_data"""
    print(time,starttime,endtime)
    atm=load_atm(time,info,starttime,endtime)
    sfc=load_sfc(time,info,starttime,endtime)
    return Bunch(sfc=sfc,atm=atm)
