#!/usr/bin/env python
# encoding: utf-8
"""
mygis.py

Series of routines for reading / writing various GIS formats 
    (added as I need them)

Requires GDAL and netCDF4 (both the libraries and the Python bindings)
    if GDAL is not available it will still import, but libraries that require gdal will raise an import error

Created by Ethan Gutmann on 2011-06-17.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""
import time
import os, glob

try:
    from osgeo import gdal,osr
    import osgeo.gdalconst as GDC
    gdalloaded=True
except:
    gdalloaded=False

NETCDF4="netCDF4"
NIO="Nio"
try:
    from netCDF4 import Dataset
    nclib=NETCDF4
except ImportError:
    import Nio
    nclib=NIO

from bunch import Bunch
import numpy as np


def ll2utm(lon,lat,zone=None,north=None):
    '''Convert a given longitude and latitude into a UTM coordinate
    
    (UTMx,UTMy,UTMz)=ll2utm(lon,lat,[zone=],[north=])
    If not specified the UTM zone is calculated from the longitude.
    If not specified north/south is calculated from the latitude.
    
    Either zone or hemisphere can be forced via zone=n or 
    north=1/0 (1=northern hemisphere, 0=southern hemisphere)
    
    Assumes a WGS84 datum/spheroid for both projections
    '''
    if not gdalloaded:
        raise ImportError("OSR not available")
    
    if zone==None:
        zone=np.int(np.ceil((lon+180)/6))
        if zone==0: zone=1
    if north==None:
        north=np.int(lat>=0)
    # create a UTM zone X projectiong reference
    utm_proj=osr.SpatialReference()
    # utm_proj.ImportFromEPSG(32613)
    utm_proj.SetUTM(zone,north)
    utm_proj.SetWellKnownGeogCS( 'WGS84' )
    
    # create a geographic reference in WGS84
    geog_ref=osr.SpatialReference()
    geog_ref.ImportFromEPSG(4326)
    
    # create a transfomration object between the two reference systems
    transform=osr.CoordinateTransformation(geog_ref,utm_proj)
    # transform the input coordinates to lat lon
    xy=transform.TransformPoint(lon,lat)
    return xy


def proj2ll(x=None,y=None,points=None,proj=None):
    '''Convert an arbitrary projection point to latitude and longitude. 
    
    (lon,lat)=proj2ll(x,y,proj)
    
    where proj is either: 
        <osgeo.osr.SpatialReference> object
        Well Known Text (WKT)
        Proj4
        EPSG
    
    Assumes a WGS84 datum/spheroid for both
    '''
    if not gdalloaded:
        raise ImportError("OSR not available")
    
    if proj==None:
        raise TypeError("Must specify a projection")
    
    if type(proj)==str:
        projWKT=proj
        proj=osr.SpatialReference()
        try:
            proj.ImportFromWkt(projWKT)
        except:
            proj.ImportFromProj4(projWKT)
            
    if type(proj)==int:
        projEPSG=proj
        proj=osr.SpatialReference()
        proj.ImportFromEPSG(projEPSG)
    
    # create a geographic reference in WGS84
    geog_ref=osr.SpatialReference()
    geog_ref.ImportFromEPSG(4326)
    
    # create a transfomration object between the two reference systems
    transform=osr.CoordinateTransformation(proj,geog_ref)
    # transform the input coordinates to lat lon
    if points!=None:
        data=np.array(transform.TransformPoints(points))
        lon=data[:,0]
        lat=data[:,1]
    else:
        lon,lat,z=transform.TransformPoint(x,y)
        
    return lat,lon


def ll2proj(lon=None,lat=None,points=None,proj=None):
    '''Convert latitude and longitude to an arbitrary projection point. 
    
    (x,y)=ll2proj(lon,lat,proj)
    
    where proj is either: 
        <osgeo.osr.SpatialReference> object
        Well Known Text (WKT)
        Proj4
        EPSG
    
    Assumes a WGS84 datum/spheroid for both
    '''
    if not gdalloaded:
        raise ImportError("OSR not available")
    
    if proj==None:
        raise TypeError("Must specify a projection")
    
    if type(proj)==str:
        projWKT=proj
        proj=osr.SpatialReference()
        try:
            proj.ImportFromWkt(projWKT)
        except:
            proj.ImportFromProj4(projWKT)
            
    if type(proj)==int:
        projEPSG=proj
        proj=osr.SpatialReference()
        proj.ImportFromEPSG(projEPSG)
    
    # create a geographic reference in WGS84
    geog_ref=osr.SpatialReference()
    geog_ref.ImportFromEPSG(4326)
    
    # create a transfomration object between the two reference systems
    transform=osr.CoordinateTransformation(geog_ref,proj)
    # transform the input coordinates to lat lon
    if points!=None:
        data=np.array(transform.TransformPoints(points))
        x=data[:,0]
        y=data[:,1]
    else:
        x,y,z=transform.TransformPoint(lon,lat)
        
    return x,y



def utm2ll(utmx,utmy,zone=13,north=1):
    '''Convert a given UTM point to latitude and longitude. 
    
    (lon,lat,elev)=utm2ll(utmx,utmy,[zone=],[north=])
    
    Defaults to Zone 13 Northern hemisphere, but either can be set.
    
    Assumes a WGS84 datum/spheroid for both
    '''
    if not gdalloaded:
        raise ImportError("OSR not available")
    # create a UTM zone X projectiong reference
    utm_proj=osr.SpatialReference()
    # utm_proj.ImportFromEPSG(32613)
    utm_proj.SetUTM(zone,north)
    utm_proj.SetWellKnownGeogCS( 'WGS84' )
    
    # create a geographic reference in WGS84
    geog_ref=osr.SpatialReference()
    geog_ref.ImportFromEPSG(4326)
    
    # create a transfomration object between the two reference systems
    transform=osr.CoordinateTransformation(utm_proj,geog_ref)
    # transform the input coordinates to lat lon
    latlon=transform.TransformPoint(utmx,utmy)
    return latlon


def read_img(filename):
    if gdalloaded:
        return gdal.Open(filename).ReadAsArray()
    else:
        raise ImportError("GDAL not available")

def read_tiff(filename,xmin=0,xmax=None,ymin=0,ymax=None,bounds=None):
    '''read a GeoTIFF and return the data and geographic information
    
    output is a structure :
        data:tiff data as an array
        geo: geographic transform data (top-left-X,dx,rot,top-left-y,rot,dy)
        proj:OSR spatial reference for the projection, use proj.ExportToPrettyWkt() to view
        topleft:topleft grid cell coordinate (x,y)
        dx:grid resolution in E-W direction
        dy:grid resolution in N-S direction
    '''
    if not gdalloaded:
        raise ImportError("GDAL not available")

    dataset = gdal.Open(filename, GDC.GA_ReadOnly)
    xsize=None
    ysize=None
    if xmax:
        xsize=xmax-xmin+1
    if ymax:
        ysize=ymax-ymin+1
    data=dataset.ReadAsArray(xoff=xmin,yoff=ymin,xsize=xsize,ysize=ysize)
    geo=dataset.GetGeoTransform()
    projWKT=dataset.GetProjection()
    proj=osr.SpatialReference()
    proj.ImportFromWkt(projWKT)

    return Bunch(proj=proj,geo=geo,data=data,topleft=(geo[0],geo[3]),dx=geo[1],dy=geo[5])
    #     GeoTransform[0] /* top left x */
    #     GeoTransform[1] /* w-e pixel resolution */
    #     GeoTransform[2] /* rotation, 0 if image is "north up" */
    #     GeoTransform[3] /* top left y */
    #     GeoTransform[4] /* rotation, 0 if image is "north up" */
    #     GeoTransform[5] /* n-s pixel resolution */


def write_tiff(filename,data,res=None,origin=None,zone=None):
    '''write a geotiff:WARNING not well tested and likely to break!
    
    write_tiff(filename,data,res=,origin=,zone=)
        for now all inputs are required
        filename= name of outputfile (string)
        data=output data (2D array of numbers)
        res=resolution (number)
        origin=UTM coordinate of top left [NE] gridcell (number,number)
        zone=UTM zone (number)
    '''
    if res==None:
        print("Error: must set output resolution")
    if origin==None:
        print("Error: must set origin (in UTM projection)")
    if zone==None:
        print("Error: For now, UTM Only, and you must set the zone")
    
    format = "GTiff"
    if not gdalloaded:
        raise ImportError("GDAL not available")
        
    driver = gdal.GetDriverByName( format )
    metadata = driver.GetMetadata()
    
    if data.dtype==np.int8:
        datatype=gdal.GDT_Byte
    elif data.dtype==np.int16:
        datatype=gdal.GDT_Int16
    elif data.dtype==np.int32:
        datatype=gdal.GDT_Int32
    elif data.dtype==np.uint16:
        datatype=gdal.GDT_UInt16
    elif data.dtype==np.uint32:
        datatype=gdal.GDT_UInt32
    elif data.dtype==np.float32:
        datatype=gdal.GDT_Float32
    elif data.dtype==np.float64:
        datatype=gdal.GDT_Float64
    else:
        print("UNKNOWN Datatype:"+str(data.dtype)+" Please add it to the code")
    
    xs=data.shape[1]
    ys=data.shape[0]
    dst_ds = driver.Create( filename, xs, ys, 1, datatype )
    dst_ds.SetGeoTransform( [ origin[0], res, 0, origin[1], 0, -1*res ] )
    
    srs = osr.SpatialReference()
    srs.SetUTM( zone, 1 )
    srs.SetWellKnownGeogCS( 'WGS84' )
    dst_ds.SetProjection( srs.ExportToWkt() )
    
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds=None


def read_geo(filesearch,outputdim=2):
    latnames=["lat","latitude","XLAT","XLAT_M"]
    lonnames=["lon","longitude","XLONG","XLONG_M"]
    latdat=None
    londat=None
    
    filename=glob.glob(filesearch)[0]
    
    for lat,lon in zip(latnames,lonnames):
        try:
            latdat=read_nc(filename,lat).data
            londat=read_nc(filename,lon).data
            
            try:
                latatts = read_atts(filename,lat)
                lonatts = read_atts(filename,lon)
            except:
                latatts=None
                lonatts=None
        except Exception as e:
            pass
    
    if latdat==None:
        # we probably weren't looking at a netcdf file...
        raise IOError("Unable to read latitude or longitude data from file:"+filename)
    
    if londat.max()>180:
        londat[londat>180]=londat[londat>180]-360
    if (len(londat.shape)==1) and (outputdim>1):
        londat,latdat=np.meshgrid(londat,latdat)
    if outputdim==3:
        londat=londat[np.newaxis,...]
        latdat=latdat[np.newaxis,...]
    
    return Bunch(lat=latdat,lon=londat, latatts=latatts, lonatts=lonatts)


def read_files(pattern,var="data",returnNCvar=False,axis=None,catch_exceptions=False,adddim=False,verbose=False):
    if type(pattern)==list:
        files=pattern
    else:
        files=glob.glob(pattern)
    files.sort()
    d=[]
    for f in files:
        if verbose:
            print(f)
        try:
            d.append(read_nc(f,var=var,returnNCvar=returnNCvar).data)
            if adddim:
                d[-1]=d[-1][np.newaxis,...]
        except Exception as e:
            if not catch_exceptions:
                raise e
            else:
                pass
    if axis!=None:
        d=np.concatenate(d,axis=axis)
    return d
    
def read_dims(filename,varname=None):
    """docstring for read_attr"""
    if nclib==NETCDF4:
        ncfile=Dataset(filename)
        if varname==None:
            dims=[str(d) for d in ncfile.dimensions]
        else:
            dims=[str(d) for d in ncfile.variables[varname].dimensions]
        
    # elif nclib==NIO:
    #     ncfile=Nio.open_file(filename,format="nc")
        # if varname==None:
        #     dims=[str(d) for d in ncfile.dimensions]
        # else:
        #     dims=[str(d) for d in ncfile.variables[varname].dimensions]
            
    ncfile.close()
    return dims


def read_attr(filename,attribute_name,varname=None):
    """docstring for read_attr"""
    if nclib==NETCDF4:
        ncfile=Dataset(filename)
        if varname==None:
            att_val=ncfile.getncattr(attribute_name)
        else:
            att_val=ncfile.variables[varname].getncattr(attribute_name)
        
    elif nclib==NIO:
        ncfile=Nio.open_file(filename,format="nc")
        if varname==None:
            att_val=ncfile.__dict__[attribute_name]
        else:
            att_val=ncfile.variables[varname].__dict__[attribute_name]
            
    ncfile.close()
    return att_val

def read_atts(filename,varname=None,global_atts=False):
    if varname==None and not global_atts:
        print("WARNING: no variable name or global attributes requested")
        return Bunch()
    
    atts=Bunch()
    if nclib==NETCDF4:
        ncfile=Dataset(filename)
    elif nclib==NIO:
        ncfile=Nio.open_file(filename,format="nc")
        
    if varname!=None:
        data=ncfile.variables[varname]
    if global_atts:
        data=ncfile
    
    
    for k in data.__dict__.keys():
        atts[k]=data.__dict__[k]
    
    ncfile.close()
    return atts
    

def read_nc(filesearch,var="data",proj=None,returnNCvar=False, format="NETCDF4"):
    '''read a netCDF file and return the specified variable

    output is a structure :
        data:raw data as an array
        proj:string representation of the projection information
        atts:data attribute dictionary (if any)
    if (returnNCvar==True) then the netCDF file is not closed and the netCDF4 
        representation of the variable is returned instead of being read into 
        memory immediately.  
    '''
    filename=glob.glob(filesearch)[0]
    if not filename: 
        print(filesearch)
    if nclib==NETCDF4:
        d=Dataset(filename, mode='r',format=format)
    else:
        d=Nio.open_file(filename,format="nc")
    outputdata=None
    if var != None:
        data=d.variables[var]
        attributes=d.variables[var].__dict__
        if returnNCvar:
            outputdata=data
        else:
            if len(data.shape)==0:
                if nclib==NIO:
                    outputdata=data.get_value()
                else:
                    outputdata=data[:]
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


def _write_one_var(NCfile,data,varname="data",units=None,dtype='f',dims=('t','z','y','x'),attributes=None, record_dim=None):
    dimlen=list(data.shape)
    if record_dim==None:
        for i in range(len(dimlen)):
            if dims[i]=="time":dimlen[i]=None
    else:
        dimlen[record_dim]=None
    dimlen=tuple(dimlen)
    for n,dim in zip(dimlen,dims):
        NCfile.createDimension(dim, n)
    try:
        fill_value=attributes["_FillValue"]
    except:
        fill_value=None
    NCfile.createVariable(varname,dtype,dims,fill_value=fill_value)
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units
    if attributes:
        for k in attributes.keys():
            if k!="_FillValue":
                NCfile.variables[varname].__setattr__(k,attributes[k])


def addvar(NCfile,data,varname,dims,dtype='f',attributes=None, record_dim=None):
    for i,d in enumerate(dims):
        if not(d in NCfile.dimensions):
            if (record_dim!=None) and (i==record_dim):
                NCfile.createDimension(d,None)
            else:
                NCfile.createDimension(d,data.shape[i])
    
    try:
        fill_value=attributes["_FillValue"]
    except:
        fill_value=None
    NCfile.createVariable(varname,dtype,dims,fill_value=fill_value)
    NCfile.variables[varname][:]=data.astype(dtype)
    if attributes:
        for k in attributes.keys():
            if k!="_FillValue":
                NCfile.variables[varname].__setattr__(k,attributes[k])

def write(filename,data,dtype='f',varname="data",dims=None,units=None,attributes=None, record_dim=None,
          lat=None,lon=None,extravars=None,history="",global_attributes=None, format="NETCDF4"):
    """write a netcdf file 
    
    filename = name of output netcdf file (.nc will be appended automatically)
    data = data to write to file
    dtype = data type of data (see below)
    varname = name of variable to create in the netcdf file
    units = units for the data
    lat = a latitude variable to add (here for legacy reasons, use extravars)
    lon = ditto    
    attribues: bunch/dictionary with key/values pairs to be added as attributes
    extravars is a list of variables to add
        each variable is a bunch or other class with attribues:
            data = data to write
            name = name of variable
            dims = dimensions to use: i.e. [[[[t],z],y],x] matching dimensions in primary data
            dtype= data type of output variable
                'd': 64 bit float
                'f': 32 bit float
                'l': long
                'i': 32 bit integer
                'h': 16 bit integer
                'b': 8 bit integer
                'S1': character
            attribues: bunch/dictionary with key/values pairs to be added as attributes
    """
    history = 'Created : ' + time.ctime() +'\nusing simple io.write by:'+os.environ['USER']+"  "+history
    
    
    if dims==None:
        if len(data.shape)==1:
            dims=('x',)
        if len(data.shape)==2:
            dims=('y','x')
        if len(data.shape)==3:
            dims=('z','y','x')
        if len(data.shape)==4:
            dims=('t','z','y','x')
        if len(data.shape)>4:
            dims=tuple([chr(dimname) for dimname in range(97,97+len(data.shape))])
    if nclib==NIO:
        return _write_nio(filename,data,dtype,varname,dims,units,attributes,lat,lon,extravars,history,global_attributes)
    else:
        NCfile=Dataset(filename,mode="w",format=format)
        _write_one_var(NCfile,data,varname=varname,units=units,dtype=dtype,dims=dims,attributes=attributes,record_dim=record_dim)
    
    if lat!=None:
        print("WARNING: explicitly passing lat is deprecated, use extravars")
        if len(lat.shape)>1:
            NCfile.createVariable("lat",'f',(dims[-2],dims[-1]))
        else:
            NCfile.createVariable("lat",'f',(dims[-2],))
        NCfile.variables["lat"][:]=lat.astype('f')
    if lon!=None:
        print("WARNING: explicitly passing lon is deprecated, use extravars")
        if len(lon.shape)>1:
            NCfile.createVariable("lon",'f',(dims[-2],dims[-1]))
        else:
            NCfile.createVariable("lon",'f',(dims[-1],))
        NCfile.variables["lon"][:]=lon.astype('f')
    
    if extravars:
        for e in extravars:
            if 'record_dim' in e.keys():
                this_record_dim=e.record_dim
            else:
                this_record_dim=None
            addvar(NCfile,e.data,e.name,e.dims,e.dtype,e.attributes, this_record_dim)

    if global_attributes!=None:
        for k in global_attributes.keys():
            if k=="history":
                NCfile.__setattr__(k,history+global_attributes[k])
            else:
                NCfile.__setattr__(k,global_attributes[k])
    else:
        NCfile.history=history
    NCfile.close()

    
def _write_nio(filename,data,dtype,varname,dims,units,attributes,lat,lon,extravars,history,global_attributes):
    NCfile=Nio.open_file(filename,mode="w",format="nc",history=history)
    
    for n,dim in zip(data.shape,dims):
        NCfile.create_dimension(dim, n)
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    
    if extravars:
        for e in extravars:
            _nio_addvar(NCfile,e.data,e.name,e.dims,e.dtype,e.attributes)

    if global_attributes!=None:
        for k in global_attributes.keys():
            NCfile.__setattr__(k,global_attributes[k])
    
    
    NCfile.close()


class NC_writer(object):
    NCfile=None
    curVar=None
    nx=0
    ny=0
    nz=None
    def __init__(self, filename,nx,ny,nz=None,var=None,dtype='f'):
        history = 'Created : ' + time.ctime() + '\nby:'+os.environ['USER']+" using NC_writer Class"
        self.NCfile=Dataset(filename,mode='w',format="NETCDF4",history=history)
        self.NCfile.createDimension('time', 0)
        self.NCfile.createDimension('lat', ny)
        self.ny=ny
        self.NCfile.createDimension('lon', nx)
        self.nx=nx
        if nz:
            self.NCfile.createDimension('level',nz)
            self.nz=nz
        self.NCfile.createVariable('time','l',('time',))
        if var: self.addVar(var,dtype=dtype)
            
    
    def addVar(self,varname,dtype='f'):
        if self.NCfile:
            if self.nz:
                self.NCfile.createVariable(varname,dtype,('time','level','lat','lon'))
            else:
                self.NCfile.createVariable(varname,dtype,('time','lat','lon'))
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
            
