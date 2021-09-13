from netCDF4 import Dataset
import xarray as xr
import numpy as np
import math
from genNetCDF import fixType
from sys import exit

class Forcing:
    nx = 2
    ny = 2
    nz = 10
    # from ideal test
    sealevel_pressure = 100000.0 # pressure at sea level [Pa]
    hill_height       = 1000.0   # height of the ideal hill(s) [m]
    n_hills           = 1.0      # number of hills across the domain
    dz_value          = 500.0    # thickness of each model gridcell   [m]
    surface_z         = 500.0    # Wont work with dry bv test
    z_interface       = None
    z                 = None
    dz_mass           = None
    pressure          = None

    attributes = {
        "history": "Mon Dec 01 00:00:00 2020: ncrcat 2000-01-01_00:00:00 2000-01-31_0:00:00 "
    }

    variableNames = ["u", "v", "theta", "p", "lat", "lon"] # others


    def __init__(self, nz=10, nx=2, ny=2, sealevel_pressure=100000.0,
                 rh=0.9, u_val=0.5, v_val=0.5, w_val=0.0,
                 water_vapor_val=0.001, theta_val=300.0):

        # self.forcing_f = Dataset("forcing.nc", "w", format="NETCDF4")
        # self.createDimensions()
        # self.createAttributes()

        self.nz = nz
        self.nx = nx
        self.ny = ny

        dimensions3d = {
            "level": nz,
            "lat": nx,
            "lon": ny
        }
        dimensions2d = {
            "lat": nx,
            "lon": ny
        }

        self.dimensions3d = dimensions3d
        self.dimensions2d = dimensions2d

        lon = xr.DataArray(np.arange(dimensions3d['lon']), dims=['lon'])
        lat = xr.DataArray(np.arange(dimensions3d['lat']), dims=['lat'])

        print(lon)
        print("---")
        print(lat)
        print("---")
        print("nznxny = ", nz,nx,ny)
        # self.dz_mass = self.gen3dDataArray(self.dz_value)
        # self.dz_mass[:,0,:] = self.dz_mass[:,0,:]/2

        dz = self.gen3dDataArray(self.dz_value)
        z = self.gen3dDataArray(0)
        height = self.gen2dDataArray(0)
        for i in range(1,nz):
            z[i,:,:] = z[i-1,:,:] + dz[i,:,:]

        # calc pressure
        pressure = self.gen3dDataArray(0)
        pressure = sealevel_pressure * (1 - 2.25577E-5 * z)**5.25588

        # calc potential temperature
        potential_temperature = self.gen3dDataArray(theta_val)

        # copied from exner_function
        exner = self.calc_exner(pressure)
        temperature = exner * potential_temperature
        water_vapor = self.gen3dDataArray(water_vapor_val)
        water_vapor = self.set_water_vapor(water_vapor, temperature, pressure)

        data_vars = dict(
            height = height,
            dz = dz,
            z = z,
            pressure = pressure,
            temperature = temperature,
            potential_temperature = potential_temperature,
            water_vapor = water_vapor,
            u = self.gen3dDataArray(u_val),
            v = self.gen3dDataArray(v_val),
            w = self.gen3dDataArray(w_val),
            lat = lat,
            lon = lon
            )

        ds = xr.Dataset(
            data_vars = data_vars,
            attrs = self.attributes
        )
        ds.to_netcdf("forcing.nc", "w", "NETCDF4")


    def gen3dDataArray(self, val):
        data = xr.DataArray(np.full((self.nz,self.nx,self.ny), val),
                            dims=self.dimensions3d.keys())
        return data

    def gen2dDataArray(self, val):
        data = xr.DataArray(np.full((self.nx,self.ny), val),
                            dims=self.dimensions2d.keys())
        return data

    def getDimension(self):
        return self.dimensions

    def getAttributes(self):
        return self.attributes

    def getVariableNames(self):
        return self.variableNames

    def createDimensions(self):
        for name, val in self.dimensions.items():
            self.forcing_f.createDimension(name,val)

    def createAttributes(self):
        for name, val in self.attributes.items():
            fixed_val = fixType(val)
            self.forcing_f.setncattr(name, fixed_val)

    def createVariables(self):
        for var in self.variableNames:
            print(var)

    def calc_exner(self,pressure):
        po=100000; Rd=287.058; cp=1003.5
        return (pressure / po) ** (Rd/cp)


    def set_water_vapor(self, water_vapor, temperature, pressure):
        for k in range(1,self.nz):
            for i in range(0,self.nx):
                for j in range(0,self.ny):
                    water_vapor[k,i,j] = sat_mr(temperature[k,i,j],
                                                pressure[k,i,j])
        return water_vapor


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
