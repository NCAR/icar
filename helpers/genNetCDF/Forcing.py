from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import numpy as np
import datetime
import math
from genNetCDF import fixType
from sys import exit

# Create NetCDF file containing the forcing data
class Forcing:
    # from ideal test
    attributes = {"history": "Dec 01 00:00:00 2020"}


    def __init__(self, nz=10, nx=2, ny=2, sealevel_pressure=100000.0,
                 rh=0.9, u_val=0.5, v_val=0.5, w_val=0.0,
                 water_vapor_val=0.001, theta_val=300.0, nt=2,
                 height_value=500, dx=10, dy=10, dz_value=500.0):

        self.setup_class_variables(nz, nx, ny, nt, sealevel_pressure)

        # --------------------------------------------------------------
        # Create and define variables for datafile
        # --------------------------------------------------------------
        # create time
        nt = 2
        t0 = datetime.datetime(2020,12,1)
        time = xr.DataArray([t0+datetime.timedelta(dt*100) for dt in range(nt)], name="time",
                            dims=["time"])

        # create longitude and latitude
        lat_tmp = np.zeros([nt, nx, ny])
        lon_tmp = np.zeros([nt, nx, ny])
        lat_flat = np.arange(39,39+nx*dx, dx)
        lon_flat = np.arange(-107,-107+ny*dy, dy)

        self.define_data_variables(nt, nz, nx, ny, height_value, lat_flat,
                                   lon_flat, dz_value)

        # --------------------------------------------------------------
        # Combine variables, create dataset and write to file
        # --------------------------------------------------------------
        data_vars = dict(
            u = self.u,
            v = self.v,
            theta = self.theta,
            qv = self.qv,
            height = self.height,
            z = self.z,
            pressure = self.pressure,
            lat_m = self.lat,
            lon_m = self.lon,
            time = time)

        ds = xr.Dataset(
            data_vars = data_vars,
            attrs = self.attributes
        )

        ds.to_netcdf("forcing.nc", "w", "NETCDF4", unlimited_dims='time')


    def set_water_vapor(self, water_vapor, temperature, pressure):
        for k in range(1,self.nz):
            for i in range(0,self.nx):
                for j in range(0,self.ny):
                    water_vapor[k,i,j] = sat_mr(temperature[k,i,j],
                                                pressure[k,i,j])
        return water_vapor


    def setup_class_variables(self, nz, nx, ny, nt, sealevel_pressure):
        self.nt = nt
        self.nz = nz
        self.nx = nx
        self.ny = ny
        self.sealevel_pressure = sealevel_pressure
        dimensions4d = {
            "time": nt,
            "level": nz,
            "lat": nx,
            "lon": ny
        }
        dimensions3d = {
            "level": nz,
            "lat": nx,
            "lon": ny
        }
        dimensions3d_t = {
            "time": nt,
            "lat": nx,
            "lon": ny
        }
        dimensions2d = {
            "lat": nx,
            "lon": ny
        }
        dimensions1d = {
            "time": 1

        }
        self.dimensions4d = dimensions4d
        self.dimensions3d = dimensions3d
        self.dimensions3d_t = dimensions3d_t
        self.dimensions2d = dimensions2d
        self.dimensions1d = dimensions1d


    def define_data_variables(self, nt, nz, nx, ny, height_value,lat_flat,
                              lon_flat, dz_value):
        dims2d = ["lat", "lon"]
        dims4d = ["time","level","lat", "lon"]

        # --- u variable
        self.u = xr.Variable(dims4d,
                             np.full([nt, nz, nx, ny], 0.5),
                             {'long_name':'U (E/W) wind speed', 'units':"m s**-1"})
        # --- v variable
        self.v = xr.Variable(dims4d,
                             np.full([nt, nz, nx, ny], 0.25),
                             {'long_name':'V (N/S) wind speed', 'units':"m s**-1"})

        # --- potential temperature variable
        self.theta = xr.Variable(dims4d,
                                 np.full([nt, nz, nx, ny], 270.),
                                 {'long_name':'Potential Temperature', 'units':"K"})
        # --- qv variable
        self.qv = xr.Variable(dims4d,
                              np.full([nt, nz, nx, ny], 0.1),
                              {'long_name':'Relative Humidity', 'units':"kg kg**-1"})
        # --- height
        self.height = xr.Variable(dims2d,
                                  np.full([nx, ny], height_value),
                                  {'long_name':'Topographic Height',
                                   'units':'m'})

        # --- Atmospheric Elevation
        dz = np.full([nz,nx,ny], dz_value)
        z_data = np.full([nt,nz,nx,ny], height_value)
        for k in range(1,nz):
            for i in range(0,nx):
                for j in range(0,ny):
                    z_data[:,k,i,j] = z_data[:,k-1,i,j] + dz[k,i,j]
                    self.z = xr.Variable(dims4d,
                                         z_data,
                                         {'long_name':'Atmospheric Elevation',
                                          'units':'m',
                                          'positive':'up'})

        # --- Pressure
        pressure_data = np.zeros([nt,nz,nx,ny])
        for k in range(0,nz):
            for i in range(0,nx):
                for j in range(0,ny):
                    pressure_data[:,k,i,j] = self.sealevel_pressure * \
                        (1 - 2.25577E-5 * z_data[0,k,i,j])**5.25588
                    self.pressure = xr.Variable(dims4d,
                                                pressure_data,
                                                {'long_name':'Pressure',
                                                 'units':'Pa'})
        del(pressure_data)

        # --- Latitude
        self.lat = xr.Variable(["lat"],
                               lat_flat,
                               {'long_name':'latitude',
                                'units':'degree_north'}
                               )
        # --- Longitude
        self.lon = xr.Variable(["lon"],
                               lon_flat,
                               {'long_name':'longitude',
                                'units':'degree_east'}
                               )


# Taken from atm_utilities.f90
def sat_mr(temperature,pressure):
    if (temperature < 273.15):
        a = 21.8745584
        b = 7.66
    else:
        a = 17.2693882
        b = 35.86
    e_s = 610.78 * math.exp(a * (temperature - 273.16) / (temperature - b))
    if ((pressure - e_s) <= 0):
        e_s = pressure * 0.99999
    sat_mr_val = 0.6219907 * e_s / (pressure - e_s)
    return sat_mr_val

def calc_exner(pressure):
    po=100000; Rd=287.058; cp=1003.5
    return (pressure / po) ** (Rd/cp)
