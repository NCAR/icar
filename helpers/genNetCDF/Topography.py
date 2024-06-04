import datetime as dt
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
                 lon0 = -105,
                 Schaer_test=False):
        # print("todo: nx and ny because of C") #

        # initialize program variables
        nt = 1
        self.n_hills = n_hills
        nx,ny = self.setup_class_variables(nx, ny, nt, mult_factor)
        self.setup_attributes(nx,ny)


        # --------------------------------------------------------------
        # Create and define variables for datafile
        # --------------------------------------------------------------
        # create time
        time_s = dt.datetime.strptime(str(time_start), '%Y%m%d')
        time = [time_s.strftime('%Y-%m-%d_%H:%M:%S')]

        # create longitude and latitude
        # lon_tmp = np.arange(lon0,lon0+(nx*dx),dx)[:nx] #[np.newaxis,:nx].repeat(ny,axis=0)
        # lat_tmp = np.arange(lat0,lat0+(ny*dy),dy)[:ny] #[:ny,np.newaxis].repeat(nx,axis=1)

        ## If your displacements aren't too great (less than a few kilometers) and you're not right at the poles,
        # use the quick and dirty estimate that 111,111 meters (111.111 km) in the y direction is 1 degree
        # (of latitude) and 111,111 * cos(latitude) meters in the x direction is 1 degree (of longitude).
        lon_tmp = np.arange(lon0-(nx/2*dx/111111/np.cos(np.radians(lat0))),
                    lon0+(nx/2*dx/111111/np.cos(np.radians(lat0))),
                    dx/111111/np.cos(np.radians(lat0))
                  )[:nx]
        lat_tmp = np.arange(lat0-(ny/2*dy/111111),
                            lat0+(ny/2*dy/111111),
                            dy/111111
                        )[:ny]
        # print(" lat_hi min/max: ", np.amin(lat_tmp),  np.amax(lat_tmp))
        lon_tmp, lat_tmp = np.meshgrid(lon_tmp, lat_tmp)


        i = (np.arange(self.nx) - self.nx/2) * dx
        j = (np.arange(self.ny) - self.ny/2) * dy    # dx=dy

        X, Y = np.meshgrid(i,j)

        self.define_data_variables(lat_tmp, lon_tmp, X, Y, height_value, hill_height, n_hills, Schaer_test, dx)

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
    def define_data_variables(self, lat_tmp, lon_tmp, X, Y, height_value,
                              hill_height, n_hills, Schaer_test, dx):
        # dimensions of variables
        dims2d = ["lat", "lon"]
        # dims3d = ["time","lat", "lon"]

        # --- xlat_m
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

        print("   hires lon/lat min/max:  ", np.min(lon_tmp), np.max(lon_tmp), np.min(lat_tmp), np.max(lat_tmp) )

        # --- hgt_m
        # hgt = np.full([self.nt,self.nx,self.ny], height_value)
        # hgt = self.genHill(hill_height)
        if Schaer_test==True:
            hgt=self.gen_adv_test_topo(hill_height, dx)
        elif n_hills == 1:
            hgt = self.genHill(hill_height)
        elif n_hills >1:
            hgt = self.gen_n_Hills(hill_height, n_hills)
        elif n_hills ==0:
            hgt = self.genHill(hill_height=0)


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


    # # generate a simple mountain range
    def gen_n_Hills(self, hill_height, n_hills):
        i = (np.arange(self.nx) - self.nx/2) / self.nx * np.pi * 2
        j = (np.arange(self.nx) - self.ny/2) / self.ny * np.pi * 2

        ig, jg = np.meshgrid(i,j)

        # should become arguments, but for now:
        c = 0.15                  # fraction of domain taken up by the hill(s)
        sigma = n_hills**2       # amount of cosines, (very) roughly

        hgt= (
            ( np.cos(ig/c) )**2 * np.exp(-(ig/c)**2/sigma) *
            ( np.cos(jg/c) )**2 * np.exp(-(jg/c)**2/sigma)
        ) * hill_height
        print("   generated ", n_hills," hills w max hgt: ", np.amax(hgt), " (hh=", hill_height, ")")
        return hgt


    # Topo for Schär's advection test
    def gen_adv_test_topo(self, hill_height, dx):
        i = (np.arange(self.nx) - self.nx/2) * dx  # / self.nx * np.pi * 2
        j = (np.arange(self.ny) - self.ny/2) * dx  # / self.ny * np.pi * 2  # dx=dy
        ig, jg = np.meshgrid(i,j)

        lmbda = 8000
        a     = 25000
        adv_3D = False  # should become an argument (maybe)

        if adv_3D==True:
            hgt = (
                hill_height *
                self.h_x(ig,  lmbda) * self.h_x_star( ig,  a)
                * self.h_x(jg, lmbda) * self.h_x_star( jg,  a)
            )
            j_a =  np.where(abs(j)>a)[0]  # satisfy the hgt=0 for |x|>a  condition (eqn 26b in Schaer 2002)
            hgt[j_a] = 0
        elif adv_3D==False:  # generate 2D topo:
            hgt = (
                hill_height *
                self.h_x(ig,  lmbda) * self.h_x_star( ig,  a)
            )

        # satisfy the hgt=0 for |x|>a  condition  (eqn 26b): (could be done more elegantly)
        i_a =  np.where(abs(i)>a)[0]
        hgt[:,i_a] = 0


        print("   generated Schaer Topo w max hgt: ", np.amax(hgt), " (hh=", hill_height, ")")
        return hgt

    def h_x(self, x,lmbda):
        h_x = (np.cos(np.pi*x /lmbda))**2
        return h_x

    def h_x_star(self, x, a):

        h_x = (np.cos(np.pi*x/2/a))**2
        return h_x


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
