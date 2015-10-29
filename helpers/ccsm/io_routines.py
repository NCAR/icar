import subprocess,os,glob
import numpy as np
import netCDF4
from bunch import Bunch
import gc,sys


sfcvarlist=["hfss","hfls"]
icar_sfc_var=["sensible_heat","latent_heat"]

atmvarlist=["ta","hus","ua","va"]
icar_atm_var=["t","qv","u","v"]

converted_sfc_files=[]
sfc_ncfiles=dict()

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
    outputproj=None
    if proj!=None:
        projection=d.variables[proj]
        outputproj=str(projection)
    
    
    if returnNCvar:
        return Bunch(data=outputdata,proj=outputproj,ncfile=d,atts=attributes)
    d.close()
    return Bunch(data=outputdata,proj=outputproj,atts=attributes)

def find_sfc_file(time,varname,info):
    file_base= info.sfcdir+info.sfcfile
    file_base= file_base.replace("_VAR_",varname)
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(1))
    file_base= file_base.replace("_D_","{0:02}".format(1))
    
    filestart=None
    for i in range(info.ntimes):
        if (info.times[i].year==time.year) and (filestart==None):
            filestart=i
        if (info.times[i]==time):
            offset=i-filestart
    
    print(file_base)
    return glob.glob(file_base)[0],offset*2

def find_atm_file(time,varname,info):
    file_base= info.atmdir+info.atmfile
    file_base= file_base.replace("_VAR_",varname)
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    atm_file = file_base.replace("_D_","{0:02}".format(time.day))
    
    print(atm_file)
    return glob.glob(atm_file)[0]

def load_sfc(time,info,ntimes):
    """load surface forcing from a grib file (or netcdf file if it has been converted previously)"""
    
    outputdata=Bunch()
    ny=info.ymax-info.ymin
    nx=info.xmax-info.xmin
    for s,v in zip(icar_sfc_var,sfcvarlist):
        inputfile,start=find_sfc_file(time,v,info)
        stop=start+ntimes*2
        print(inputfile, start,stop)
        sys.stdout.flush()
        nc_data=read_nc(inputfile,v)#,returnNCvar=True)
        curdata=nc_data.data[start:stop,info.ymin:info.ymax,info.xmin:info.xmax]
        # nc_data.ncfile.close()
        outputdata[s]=curdata.reshape(ntimes,2,ny,nx).mean(axis=1)
    
    return outputdata

def load_atm(time,info):
    """Load atmospheric variable from a netcdf file"""
    
    outputdata=Bunch()

    for s,v in zip(icar_atm_var,atmvarlist):
        atmfile=find_atm_file(time,v,info)
        print(atmfile)
        sys.stdout.flush()
        nc_data=read_nc(atmfile,v)#,returnNCvar=True)
        outputdata[s]=nc_data.data[:,:,info.ymin:info.ymax,info.xmin:info.xmax]
        # nc_data.ncfile.close()

    nc_data=read_nc(atmfile,"ps")#,returnNCvar=True)
    outputdata.ps=nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
    # nc_data.ncfile.close()
    del nc_data
    print(gc.collect())
    sys.stdout.flush()
    
    a=read_nc(atmfile,"a").data
    b=read_nc(atmfile,"b").data
    p0=read_nc(atmfile,"p0").data
    outputdata.p = a[np.newaxis,:,np.newaxis,np.newaxis]*p0+b[np.newaxis,:,np.newaxis,np.newaxis]*outputdata.ps[:,np.newaxis,:,:]
    
    outputdata.ntimes=outputdata.p.shape[0]
    
    return outputdata


def load_data(time,info):
    """docstring for load_data"""
    print(time)
    atm=load_atm(time,info)
    sfc=load_sfc(time,info,atm.ntimes)
    return Bunch(sfc=sfc,atm=atm)


