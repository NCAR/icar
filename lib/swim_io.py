# from netCDF4 import Dataset
import time
import os
import numpy as np
import Nio
from bunch import Bunch
import glob

def Dataset(filename,mode="r",format="nc"):
    return Nio.open_file(filename,mode=mode,format=format)

def read_files(pattern,var="data",returnNCvar=False,axis=None):
    if type(pattern)==list:
        files=pattern
    else:
        files=glob.glob(pattern)
    files.sort()
    d=[]
    for f in files:
        d.append(read_nc(f,var=var,returnNCvar=returnNCvar).data)
    if axis!=None:
        d=np.concatenate(d,axis=axis)
    return d
    

def read_nc(filename,var="data",proj=None,returnNCvar=False):
    '''read a netCDF file and return the specified variable

    output is a structure :
        data:raw data as an array
        proj:string representation of the projection information
        atts:data attribute dictionary (if any)
    if (returnNCvar==True) then the Nio file is note closed and the Nio 
        representation of the variable is returned instead of being read into 
        memory immediately.  
    '''
    d=Nio.open_file(filename, mode='r',format="nc")
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


def _write1d(NCfile,data,varname="data",units=None,dtype='f'):
    nx=data.size
    NCfile.create_dimension('x', nx)
    NCfile.create_variable(varname,dtype,('x',))
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units

def _write2d(NCfile,data,varname="data",units=None,dtype='f'):
    (ny,nx)=data.shape
    NCfile.create_dimension('x', nx)
    NCfile.create_dimension('y', ny)
    NCfile.create_variable(varname,dtype,('y','x'))
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units

def _write3d(NCfile,data,varname="data",units=None,dtype='f'):
    (nz,ny,nx)=data.shape
    NCfile.create_dimension('x', nx)
    NCfile.create_dimension('y', ny)
    NCfile.create_dimension('z', nz)
    NCfile.create_variable(varname,dtype,('z','y','x'))
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units

def _write4d(NCfile,data,varname="data",units=None,dtype='f'):
    (nt,nz,ny,nx)=data.shape
    NCfile.create_dimension('x', nx)
    NCfile.create_dimension('y', ny)
    NCfile.create_dimension('z', nz)
    NCfile.create_dimension('t', nt)
    NCfile.create_variable(varname,dtype,('t','z','y','x'))
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units


def addvar(NCfile,data,varname,dims,dtype='f',attributes=None):
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    if attributes:
        for k in attributes.keys():
            NCfile.variables[varname].__setattr__(k,attributes[k])

def write(filename,data,dtype='f',varname="data",units=None,lat=None,lon=None,extravars=None):
    """write(filename,data,dtype='f',varname="data",units=None,lat=None,lon=None,extravars=None)"""
    history = 'Created : ' + time.ctime() +'\nusing simple ncio.write by:'+os.environ['USER']
    NCfile=Nio.open_file(filename,mode="w",format="nc",history=history)
    if len(data.shape)==1:
        _write1d(NCfile,data,varname=varname,units=units,dtype=dtype)
    if len(data.shape)==2:
        _write2d(NCfile,data,varname=varname,units=units,dtype=dtype)
    if len(data.shape)==3:
        _write3d(NCfile,data,varname=varname,units=units,dtype=dtype)
    if len(data.shape)==4:
        _write4d(NCfile,data,varname=varname,units=units,dtype=dtype)
    
    if lat!=None:
        if len(lat.shape)>1:
            NCfile.create_variable("lat",'f',('y','x'))
        else:
            NCfile.create_variable("lat",'f',('y',))
        NCfile.variables["lat"][:]=lat.astype('f')
    if lon!=None:
        if len(lon.shape)>1:
            NCfile.create_variable("lon",'f',('y','x'))
        else:
            NCfile.create_variable("lon",'f',('x',))
        NCfile.variables["lon"][:]=lon.astype('f')
    
    if extravars:
        for e in extravars:
            addvar(NCfile,e.data,e.name,e.dims,e.dtype,e.attributes)
    
    NCfile.close()


class NC_writer(object):
    NCfile=None
    curVar=None
    nx=0
    ny=0
    nz=None
    def __init__(self, filename,nx,ny,nz=None,var=None,dtype='f'):
        history = 'Created : ' + time.ctime() + '\nby:'+os.environ['USER']+" using NC_writer Class"
        self.NCfile=Nio.open_file(filename,mode='w',format="nc",history=history)
        self.NCfile.create_dimension('time', 0)
        self.NCfile.create_dimension('lat', ny)
        self.ny=ny
        self.NCfile.create_dimension('lon', nx)
        self.nx=nx
        if nz:
            self.NCfile.create_dimension('level',nz)
            self.nz=nz
        self.NCfile.create_variable('time','l',('time',))
        if var: self.addVar(var,dtype=dtype)
            
    
    def addVar(self,varname,dtype='f'):
        if self.NCfile:
            if self.nz:
                self.NCfile.create_variable(varname,dtype,('time','level','lat','lon'))
            else:
                self.NCfile.create_variable(varname,dtype,('time','lat','lon'))
            self.curVar=varname
    
    def appendToVar(self,data,varname=None,date=None,pos=None,dtype='f'):
        if varname==None:varname=self.curVar
        var=self.NCfile.variables[varname]
        if pos:
            n=pos
        else:
            n=self.NCfile.dimensions['time']
            if n==None:
                n=0
        if self.nz:
            var[n,:,:,:]=data.astype(dtype)
        else:
            var[n,:,:]=data.astype(dtype)
        if date:self.NCfile.variables['time'][n]=long(date)
    
    def close(self):
        if self.NCfile:
            self.NCfile.close()
            self.NCfile=None
            
    def __del__(self):
        self.close()
        # if self.NCfile:
        #     self.NCfile.close()
        #     self.NCfile=None
        # super(NC_writer,self).__del__()
            
