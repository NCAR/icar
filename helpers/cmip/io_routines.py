import subprocess,os,glob
import gc,sys
import numpy as np
import netCDF4
from bunch import Bunch
import mygis

g=9.8
atmvarlist=["ta","hus","ua","va"]
icar_atm_var=["t","qv","u","v"]

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
            # outputdata=data[:]
            ntimes=365*4
            if len(data.shape)>2:
                outputdata=data[:ntimes,...]
            else:
                outputdata=data[:]
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
    file_base= file_base.replace("_GCM_",info.gcm_name)
    file_base= file_base.replace("_VAR_",varname)
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_EXP_",info.experiment)
    atm_file = file_base.replace("_ENS_",info.ensemble)

    print(atm_file)
    filelist = glob.glob(atm_file)
    filelist.sort()
    return filelist

def find_sst_file(time,info):
    file_base= info.atmdir+info.atmfile
    file_base= file_base.replace("_GCM_",info.gcm_name)
    file_base= file_base.replace("_VAR_","sst")
    file_base= file_base.replace("6hrLev","day")
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_EXP_",info.experiment)
    sst_file = file_base.replace("_ENS_",info.ensemble)

    print(sst_file)
    filelist=glob.glob(sst_file)
    filelist.sort()
    return filelist


def load_atm(time,info):
    """Load atmospheric variable from a netcdf file"""

    outputdata=Bunch()

    for s,v in zip(icar_atm_var,atmvarlist):
        atmfile_list=find_atm_file(time,v,info)
        for atmfile in atmfile_list:
            nc_data=read_nc(atmfile,v)#,returnNCvar=True)
            newdata=nc_data.data[:,:,info.ymin:info.ymax,info.xmin:info.xmax]
            if s in outputdata:
                outputdata[s]=np.concatenate([outputdata[s],newdata])
            else:
                outputdata[s]=newdata

    outputdata.ntimes=0
    for atmfile in atmfile_list:
        varname="ps"
        nc_data=read_nc(atmfile,varname)
        newdata=nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
        if varname in outputdata:
            outputdata[varname]=np.concatenate([outputdata[varname],newdata])
        else:
            outputdata[varname]=newdata

        varname="p"
        newdata=info.read_pressure(atmfile)[:,:,info.ymin:info.ymax,info.xmin:info.xmax]
        if varname in outputdata:
            outputdata[varname]=np.concatenate([outputdata[varname],newdata])
        else:
            outputdata[varname]=newdata

        outputdata.ntimes = outputdata.p.shape[0]

    # outputdata.times=info.read_time(atmfile)
    try:
        calendar = mygis.read_attr(atmfile_list[0], "calendar", varname="time")
    except (KeyError,IndexError):
        calendar = None

    outputdata.calendar = calendar

    return outputdata

def load_sfc(time,info):
    """docstring for load_sfc"""
    outputdata=Bunch()
    basefile=info.orog_file
    outputdata.hgt=read_nc(basefile,"orog").data[info.ymin:info.ymax,info.xmin:info.xmax]
    # sstfile=find_sst_file(time)
    # outputdata.sst=read_nc(sstfile,"sst").data[:,info.ymin:info.ymax,info.xmin:info.xmax]

    outputdata.land=np.zeros(outputdata.hgt.shape)
    landfrac=read_nc(basefile,"sftlf").data[info.ymin:info.ymax,info.xmin:info.xmax]
    outputdata.land[landfrac>=0.5]=1
    return outputdata

def load_data(time,info):
    """docstring for load_data"""
    print(time)
    atm=load_atm(time,info)
    sfc=load_sfc(time,info)
    return Bunch(sfc=sfc,atm=atm)
