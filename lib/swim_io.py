# from netCDF4 import Dataset
import time
import os
import numpy as np
import Nio
from bunch import Bunch
import glob

def Dataset(filename,mode="r",format="nc"):
    return Nio.open_file(filename,mode=mode,format=format)

def read_files(pattern,var=None,returnNCvar=False,axis=None):
    files=glob.glob(pattern)
    d=[]
    for f in files:
        d.append(read_nc(f,var=var,returnNCvar=returnNCvar).data)
    if axis!=None:
        d=np.concatenate(d,axis=axis)
    return d
    

def read_nc(filename,var=None,proj=None,returnNCvar=False):
    '''read a netCDF file and return the specified variable

    output is a structure :
        data:raw data as an array
        proj:string representation of the projection information
    '''
    d=Nio.open_file(filename, mode='r',format="nc")
    outputdata=None
    if var != None:
        data=d.variables[var]
        if returnNCvar:
            outputdata=data
        else:
            outputdata=data[:]
    outputproj=None
    if proj!=None:
        projection=d.variables[proj]
        outputproj=projection.__str__()
    
    
    if returnNCvar:
        return Bunch(data=outputdata,proj=outputproj,ncfile=d)
    d.close()
    return Bunch(data=outputdata,proj=outputproj)


def write1d(NCfile,data,varname="data",units=None):
    nx=data.size
    NCfile.create_dimension('x', nx)
    NCfile.create_variable(varname,'f',('x',))
    NCfile.variables[varname][:]=data.astype('f')
    if units!=None:
        NCfile.variables[varname].units=units

def write2d(NCfile,data,varname="data",units=None):
    (ny,nx)=data.shape
    NCfile.create_dimension('x', nx)
    NCfile.create_dimension('y', ny)
    NCfile.create_variable(varname,'f',('y','x'))
    NCfile.variables[varname][:]=data.astype('f')
    if units!=None:
        NCfile.variables[varname].units=units

def write3d(NCfile,data,varname="data",units=None):
    (nz,ny,nx)=data.shape
    NCfile.create_dimension('x', nx)
    NCfile.create_dimension('y', ny)
    NCfile.create_dimension('z', nz)
    NCfile.create_variable(varname,'f',('z','y','x'))
    NCfile.variables[varname][:]=data.astype('f')
    if units!=None:
        NCfile.variables[varname].units=units

def write(filename,data,varname="data",units=None,lat=None,lon=None):
    history = 'Created : ' + time.ctime() +'\nusing simple ncio.write by:'+os.environ['USER']
    NCfile=Nio.open_file(filename,mode="w",format="nc",history=history)
    if len(data.shape)==1:
        write1d(NCfile,data,varname=varname,units=units)
    if len(data.shape)==2:
        write2d(NCfile,data,varname=varname,units=units)
    if len(data.shape)==3:
        write3d(NCfile,data,varname=varname,units=units)
    
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
    
    NCfile.close()


class NC_writer(object):
    NCfile=None
    curVar=None
    nx=0
    ny=0
    nz=None
    def __init__(self, filename,nx,ny,nz=None,var=None):
        # self.NCfile=Dataset(filename,'w')
        history = 'Created : ' + time.ctime() + '\nby:'+os.environ['USER']+", using NC_writer Class"
        self.NCfile=Nio.open_file(filename,mode='w',format="nc",history=history)
        # self.NCfile.creator = 'NC_write Class by:'+os.environ['USER']
        self.NCfile.create_dimension('time', 0)
        self.NCfile.create_dimension('lat', ny)
        self.ny=ny
        self.NCfile.create_dimension('lon', nx)
        self.nx=nx
        if nz:
            self.NCfile.create_dimension('level',nz)
            self.nz=nz
        self.NCfile.create_variable('time','l',('time',))
        if var: self.addVar(var)
            
    
    def addVar(self,varname):
        if self.NCfile:
            if self.nz:
                self.NCfile.create_variable(varname,'f',('time','level','lat','lon'))
            else:
                self.NCfile.create_variable(varname,'f',('time','lat','lon'))
            self.curVar=varname
    
    def appendToVar(self,data,varname=None,date=None,pos=None):
        if varname==None:varname=self.curVar
        var=self.NCfile.variables[varname]
        if pos:
            n=pos
        else:
            n=self.NCfile.dimensions['time']
            if n==None:
                n=0
        if self.nz:
            var[n,:,:,:]=data.astype("f")
        else:
            var[n,:,:]=data.astype("f")
        if date:self.NCfile.variables['time'][n]=long(date)
    
    def close(self):
        if self.NCfile:
            self.NCfile.close()
            self.NCfile=None
            
    def __del__(self):
        if self.NCfile:
            self.NCfile.close()
            self.NCfile=None
        # super(NC_writer,self).__del__()
            
