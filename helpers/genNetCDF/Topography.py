import pandas as pd
import xarray as xr
import numpy as np


# Create NetCDF file containing the initial conditions topography
# Currently creates flat domain
class Topography:
    def __init__(self,
                 nz=32,
                 nx=100,
                 ny=100,
                 f_name="init.nc",
                 mult_factor=1,
                 dx=0.01,
                 dy=0.01,
                 n_hills=0.0,
                 height_value=500,
                 hill_height = 2000.0,
                 time_start = 20010101,
                 lat0 = 39.5,
                 lon0 = -105):
        print("todo: nx and ny because of C")

        # initialize program variables
        nt = 1
        self.n_hills = n_hills
        nx,ny = self.setup_class_variables(nx, ny, nt, mult_factor)
        self.setup_attributes(nx,ny)


        # --------------------------------------------------------------
        # Create and define variables for datafile
        # --------------------------------------------------------------
        # create time
        time_end   = time_start
        time_series = pd.date_range(start=str(time_start),
                                    end=str(time_end),
                                    freq='D').strftime('%Y-%m-%d_%H:%M:%S')
        time = time_series.astype(np.unicode_)

        # create longitude and latitude
        lon_tmp = np.arange(lon0,lon0+(nx*dx),dx)[:nx] #[np.newaxis,:nx].repeat(ny,axis=0)
        lat_tmp = np.arange(lat0,lat0+(ny*dy),dy)[:ny] #[:ny,np.newaxis].repeat(nx,axis=1)

        lon_tmp, lat_tmp = np.meshgrid(lon_tmp, lat_tmp)

        self.define_data_variables(lat_tmp, lon_tmp, height_value, hill_height)

        # --------------------------------------------------------------
        # Combine variables, create dataset and write to file
        # --------------------------------------------------------------
        data_vars = dict(
            lat_hi = self.lat_m,
            lon_hi = self.lon_m,
            hgt_hi = self.hgt_m,
            Times = time)

        ds = xr.Dataset(
            data_vars = data_vars,
            attrs = self.attributes
        )

        ds.to_netcdf(f_name, "w", "NETCDF4", unlimited_dims='time')



    # Define individual variables for datafile
    def define_data_variables(self, lat_tmp, lon_tmp, height_value,
                              hill_height):
        # dimensions of variables
        dims2d = ["lat", "lon"]
        # dims3d = ["time","lat", "lon"]

        # --- xlat_m
        print(lat_tmp.shape)
        print(dims2d)
        self.lat_m = xr.Variable(dims2d,
                                 lat_tmp,
                                 {'units':'degrees latitude',
                                  'description':'Latitude on mass grid',
                                  })

        # --- xlon_m
        self.lon_m = xr.Variable(dims2d,
                                 lon_tmp,
                                 {'units':'degrees longitude',
                                  'description':'Longitude on mass grid',
                                  })

        # --- hgt_m
        # hgt = np.full([self.nt,self.nx,self.ny], height_value)
        hgt = self.genHill(hill_height)
        self.hgt_m = xr.Variable(dims2d,
                                 hgt,
                                 {'units':'meters MSL',
                                  'description':'topography height',
                                  })


    def close(self):
        self.topography_f.close()


    def setup_class_variables(self, nx, ny, nt, mult_factor):
        self.nt = nt
        self.nx = round(nx * mult_factor)
        self.ny = round(ny * mult_factor)

        return self.nx, self.ny


    # hasn't been integrated and tested
    def genHill(self, hill_height):
        i = (np.arange(self.nx) - self.nx/2) / self.nx * np.pi * 2
        j = (np.arange(self.nx) - self.ny/2) / self.ny * np.pi * 2

        ig, jg = np.meshgrid(i,j)

        hgt = ((np.cos(ig)+1) * (np.cos(jg)+1))/4 * hill_height

        return hgt



    def setup_attributes(self,nx,ny):
        self.attributes = {
            "TITLE": "OUTPUT FROM CONTINUOUS INTEGRATION TEST",
            "SIMULATION_START_DATE": "0000-00-00_00:00:00",
            "WEST-EAST_GRID_DIMENSION": nx,
            "SOUTH-NORTH_GRID_DIMENSION": ny,
            "GRIDTYPE": "C",
            "DX": 1000.0,
            "DY": 1000.0
        }
