from netCDF4 import Dataset
from genNetCDF import fixType
import pandas as pd
import xarray as xr
import numpy as np
import math
from sys import exit


# Create NetCDF file containing the initial conditions topography
# Currently creates flat domain
class Topography:
    hill_height = 1000.0   # height of the ideal hill(s) [m]
    def __init__(self, nz=10, nx=100, ny=100, mult_factor=2.5, dx=0.1, dy=0.1,
                 n_hills=0.0, height_value=500, f_name="init.nc"):
        print("todo: nx and ny because of C")

        # initialize program variables
        nt = 1
        self.n_hills = n_hills
        nx,ny = self.setup_class_variables(nz, nx, ny, nt, mult_factor)
        self.setup_attributes(nx,ny)


        # --------------------------------------------------------------
        # Create and define variables for datafile
        # --------------------------------------------------------------
        # create time
        time_start = 20010101
        time_end   = time_start
        time_series = pd.date_range(start=str(time_start),
                                    end=str(time_end),
                                    freq='D').strftime('%Y-%m-%d_%H:%M:%S')
        time = time_series.astype(np.unicode_)

        # create longitude and latitude
        lon_flat = np.arange(-105,-105+(nx*dx),dx)[:nx]
        lat_flat = np.arange(40,40+(ny*dy),dy)[:ny]
        lat_tmp = np.zeros([nt, nx, ny])
        lon_tmp = np.zeros([nt, nx, ny])
        lat_tmp[:,:,:] = lat_flat.reshape([1,1,ny])
        lon_tmp[:,:,:] = lon_flat.reshape([1,nx,1])

        self.define_data_variables(lat_tmp, lon_tmp, height_value)

        # --------------------------------------------------------------
        # Combine variables, create dataset and write to file
        # --------------------------------------------------------------
        data_vars = dict(
            lat_hi = self.lat_m,
            lon_hi = self.lon_m,
            lat_u = self.lat_u,
            lon_u = self.lon_u,
            lat_v = self.lat_v,
            lon_v = self.lon_v,
            hgt_hi = self.hgt_m,
            Times = time)

        ds = xr.Dataset(
            data_vars = data_vars,
            attrs = self.attributes
        )

        ds.to_netcdf(f_name, "w", "NETCDF4", unlimited_dims='time')



    # Define individual variables for datafile
    def define_data_variables(self, lat_tmp, lon_tmp, height_value):
        # dimensions of variables
        dims2d = ["lat", "lon"]
        dims3d = ["time","lat", "lon"]

        # --- xlat_m
        self.lat_m = xr.Variable(dims3d,
                                 lat_tmp,
                                 {'FieldType':'104',
                                  'units':'degrees latitude',
                                  'description':'Latitude on mass grid',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })

        # --- xlon_m
        self.lon_m = xr.Variable(dims3d,
                                 lon_tmp,
                                 {'FieldType':'104',
                                  'units':'degrees longitude',
                                  'description':'Longitude on mass grid',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })

        # --- hgt_m
        hgt = np.full([self.nt,self.nx,self.ny], height_value)
        self.hgt_m = xr.Variable(dims3d,
                                 hgt,
                                 {'FieldType':'104',
                                  'MemoryOrder':'XY',
                                  'units':'meters MSL',
                                  'description':'topography height',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })


        # --- xlat_u
        self.lat_u = xr.Variable(dims3d,
                                 lat_tmp,
                                 {'FieldType':'104',
                                  'units':'degrees latitude',
                                  'description':'Latitude on U grid',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })

        # --- xlon_u
        self.lon_u = xr.Variable(dims3d,
                                 lon_tmp,
                                 {'FieldType':'104',
                                  'units':'degrees longitude',
                                  'description':'Longitude on U grid',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })

        # --- xlat_v
        self.lat_v = xr.Variable(dims3d,
                                 lat_tmp,
                                 {'FieldType':'104',
                                  'units':'degrees latitude',
                                  'description':'Latitude on V grid',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })

        # --- xlon_v
        self.lon_v = xr.Variable(dims3d,
                                 lon_tmp,
                                 {'FieldType':'104',
                                  'units':'degrees longitude',
                                  'description':'Longitude on V grid',
                                  'stagger':'M',
                                  'sr_x':'1',
                                  'sr_y':'1'
                                  })



    def close(self):
        self.topography_f.close()


    def setup_class_variables(self, nz, nx, ny, nt, mult_factor):
        self.nt = nt
        self.nz = nz
        self.nx = round(nx * mult_factor)
        self.ny = round(ny * mult_factor)

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
        return self.nx, self.ny


    # hasn't been integrated and tested
    def genHill(self):
        k = 0
        ids = 1
        ide = self.nx + 1
        jds = 1
        jde = self.ny + 1
        for i in range(ids,ide):
            for j in range(jds,jde):
                sine_curve = (math.sin((i-ids)/ \
                              np.single((ide-ids)/self.n_hills) \
                              * 2*3.14159 - 3.14159/2) + 1) / 2 \
                              * (math.sin((j-jds)/np.single((jde-jds) \
                                / self.n_hills) * 2*3.14159 - 3.14159/2)+1)/2
                self.z_interface[k,i-1,j-1] = self.surface_z + sine_curve * self.hill_height
        self.z[k,:,:] = self.z_interface[k,:,:] + self.dz_value/2


    def setup_attributes(self,nx,ny):
        self.attributes = {
            "TITLE": "OUTPUT FROM CONTINUOUS INTEGRATION TEST",
            "SIMULATION_START_DATE": "0000-00-00_00:00:00",
            "WEST-EAST_GRID_DIMENSION": nx,
            "SOUTH-NORTH_GRID_DIMENSION": ny,
            "BOTTOM-TOP_GRID_DIMENSION": 0,
            "GRIDTYPE": "C",
            "DX": 100.0,
            "DY": 100.0
        }
