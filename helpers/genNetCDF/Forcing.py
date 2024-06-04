import datetime
import numpy as np
import xarray as xr
import math

# Create NetCDF file containing the forcing data
class Forcing:
    # from ideal test
    attributes = {"history": "Dec 01 00:00:00 2020"}


    def __init__(self,nt=10, nz=10, nx=2, ny=2, sealevel_pressure=100000.0,
                 u_val=0.5, v_val=0.5, w_val=0.0,
                 water_vapor_val=0.001, theta_val=300.0,
                 height_value=500, dx=10, dy=10, dz_value=500.0,
                 qv_val=0.1, weather_model='basic',
                 pressure_func='calc_pressure_from_sea',
                 hill_height=0, lat0 = 39.5,
                 lon0 = -105,
                 Schaer_test=False):
        print(weather_model.capitalize(), "weather model in use")
        self.setup_class_variables(nz, nx, ny, nt, sealevel_pressure)

        # --- Create and define variables for datafile
        # lat_flat = np.arange(39,39+nx*dx, dx)
        # lon_flat = np.arange(-107,-107+ny*dy, dy)

        # center around lat0, lon0 (i.s.o. corner), so we can align with hi-res grid.
        lon_flat = np.arange(lon0-(nx/2*dx/111111/np.cos(np.radians(lat0))),
                    lon0+(nx/2*dx/111111/np.cos(np.radians(lat0))),
                    dx/111111/np.cos(np.radians(lat0))
                   )[:nx]
        lat_flat = np.arange(lat0-(ny/2*dy/111111),
                            lat0+(ny/2*dy/111111),
                            dy/111111
                        )[:ny]

        x_m = np.arange(-nx*dx/2,nx*dx/2, dx)


        print( "   forcing lon/lat min/max:  ", np.min(lon_flat), np.max(lon_flat), np.min(lat_flat), np.max(lat_flat) )
        # print(" forcing lat min/max: ", np.amin(lat_flat), np.amax(lat_flat))
        # lon_flat, lat_flat = np.meshgrid(lon_tmp, lat_tmp)

        self.define_data_variables(nt, nz, nx, ny, height_value, lat_flat,
                                   lon_flat, dz_value, theta_val, u_val,
                                   v_val, qv_val, weather_model, pressure_func,
                                   hill_height, Schaer_test, dx, x_m)

        # define time
        t0 = datetime.datetime(2020,12,1)
        time = xr.DataArray([t0+datetime.timedelta(hours=dt) for dt in range(nt)], name="time",
                            dims=["time"])

        # --- Write all variable to netcdf file
        self.write_netcdf_file(time)


    # Combine variables, create dataset and write to file
    def write_netcdf_file(self, time):
        data_vars = dict(
            u = self.u,
            v = self.v,
            theta = self.theta,
            qv = self.qv,
            height = self.height,
            z = self.z,
            pressure = self.pressure,
            temperature = self.temperature,
            lat_m = self.lat,
            lon_m = self.lon,
            x_m = self.x_m,
            time = time)

        ds = xr.Dataset(
            data_vars = data_vars,
            attrs = self.attributes
        )

        ds.to_netcdf("forcing.nc", "w", "NETCDF4", unlimited_dims='time', encoding={'time': {'dtype': 'i4'}})


    def set_water_vapor(self, water_vapor, temperature, pressure):
        water_vapor = sat_mr(temperature, pressure)
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
            "lat": ny,
            "lon": nx
        }
        dimensions3d = {
            "level": nz,
            "lat": ny,
            "lon": nx
        }
        dimensions3d_t = {
            "time": nt,
            "lat": ny,
            "lon": nx
        }
        dimensions2d = {
            "lat": ny,
            "lon": nx
        }
        dimensions1d = {
            "time": 1

        }
        self.dimensions4d = dimensions4d
        self.dimensions3d = dimensions3d
        self.dimensions3d_t = dimensions3d_t
        self.dimensions2d = dimensions2d
        self.dimensions1d = dimensions1d
        self.dims4d = list(dimensions4d.keys())
        self.dims3d = list(dimensions3d.keys())
        self.dims2d = list(dimensions2d.keys())
        self.dims1d = list(dimensions1d.keys())



    def define_data_variables(self, nt, nz, nx, ny, height_value,lat_flat,
                              lon_flat, dz_value, theta_val, u_val, v_val,
                              qv_val, weather_model, pressure_func, hill_height,
                              Schaer_test, dx, x_m):

        # --- u variable
        # if advection test is selected, set the appropriate windfield:
        if Schaer_test==True:
            z1 = 4000. ; z2 = 5000. ; hill_height = 3000.0  ; u0=10
            u_val = np.array( [0]* int(z1/dz_value)
                + [u0* (np.sin(np.pi/2*(z1/dz_value+1 - z1/dz_value) / ((z2-z1)/dz_value) ))**2 ]
                + [u0* (np.sin(np.pi/2*(z1/dz_value+2 - z1/dz_value) / ((z2-z1)/dz_value) ))**2 ]
                + [u0] * nz #int(nz-z2/dz_value)
            )
            u_array=np.tile(u_val[:nz], (nt,nx,ny,1) )
            u_array = np.transpose(u_array,(0,3,2,1) )  # order?

        # if uval is given as a single float, make a uniform windfield:
        elif isinstance(u_val, float):
            u_array= np.full([nt, nz, ny, nx], u_val[:nz])
        # if u_val is given as a vector, interpret this a vector in the z direction (bottom-top):
        elif isinstance(u_val, np.ndarray):
            u_array=np.tile(u_val[:nz], (nt,nx,ny,1) )
            u_array = np.transpose(u_array,(0,3,2,1) )
            # print(U.shape)
            print("   Treating u_test_val as a u field in z-direction")

        self.u = xr.Variable(self.dims4d,
                             u_array,
                             {'long_name':'U (E/W) wind speed', 'units':"m s**-1"})

        # --- v variable
        if Schaer_test==True: v_val=0.
        self.v = xr.Variable(self.dims4d,
                             np.full([nt, nz, ny, nx], v_val),
                             {'long_name':'V (N/S) wind speed', 'units':"m s**-1"})

        # --- height
        self.height = xr.Variable(self.dims2d,
                                  np.full([ny, nx], height_value),
                                  {'long_name':'Topographic Height',
                                   'units':'m'})

        # --- Atmospheric Elevation
        dz = np.full([nz,ny,nx], dz_value)
        z_data = np.full([nt,nz,ny,nx], height_value)
        # dz[0,:,:] = [50.]
        # dz[1,:,:] = [75.]
        # dz[2,:,:] = [125.]
        # dz[3,:,:] = [200.]
        # dz[4,:,:] = [300.]
        # dz[5,:,:] = [400.]

        for k in range(1,nz):
            z_data[:,k,:,:] = z_data[:,k-1,:,:] + dz[k,:,:]
        self.z_data = z_data
        self.z = xr.Variable(self.dims4d,
                             z_data,
                             {'long_name':'Atmospheric Elevation',
                              'units':'m',
                              'positive':'up'})
        del(z_data)

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

        # --- x_m
        self.x_m = xr.Variable(["x_m"],
                                x_m,
                                {'long_name':'x distance from domain center',
                                'units':'meters'}
                                )


        # --- potential temperature variable
        self.set_theta(theta_val, weather_model)

        # --- Pressure
        self.set_pressure(weather_model, pressure_func)

        # --- Temperature
        self.set_temperature(weather_model)

        # --- qv variable
        if Schaer_test==True:
            # create a small blob of moisture, in an otherwise dry environment. Values from Schaer et al 2002
            qv_arr = np.zeros([nt,nz,ny,nx])
            z0 = int(9000/dz_value)
            x0 = int(-50000/dx + nx/2)  # -50km
            Ax = int(25000/dx)
            Az = int(3000/dz_value)
            print("   setting up advection test with a cloud of qv with half-width ",Ax,"km")
            if x0-Ax<0:
                print("   QV blob outside forcing domain; increase nx_lo and or dx_lo (currently",nx, " and ", dx)
                print("   x0-Ax=", x0-Ax)
            for r in np.arange(1,0,-0.05):
                for t in np.arange(0,np.pi*2,0.1):
                    qv_arr[0,
                            z0-int(r*Az*np.sin(t)):z0+int(r*Az*np.sin(t)),
                            :,
                            x0-int(np.cos(t)*r*Ax):x0+int(np.cos(t)*r*Ax),
                        ] = (np.cos(np.pi*r/2))**2 *qv_val
            print("   qv_arr min: ",np.amin(qv_arr), "  max:", np.amax(qv_arr))

        else:  # homogenous qv throughout domain
            qv_arr = np.full([nt,nz,ny,nx], qv_val)

        self.qv = xr.Variable(self.dims4d,
                            #   np.full([nt, nz, nx, ny], ),
                            qv_arr,
                              {'long_name':'Relative Humidity',
                               'units':"kg kg**-1"})



    def set_theta(self, theta_val, model='basic'):
        if model in ['basic']:
            theta = np.full([self.nt, self.nz, self.ny, self.nx], theta_val)
        elif model in ['WeismanKlemp']:
            print('Note: theta value of', theta_val,
                  'has been replaced with a profile generated for', model,
                  'model')
            theta = np.zeros([self.nt, self.nz, self.ny, self.nx])
            theta[0,:,:,:] = np.vectorize(calc_wk_theta)(self.z_data[0,:,:,:])
            theta[:,:,:,:] = theta[0,:,:,:]
        self.theta = xr.Variable(self.dims4d, theta,
                                 {'long_name':'Potential Temperature',
                                  'units':"K"})

    def set_pressure(self, model='basic',
                     pressure_func='calc_pressure_from_sea'):
        print('Pressure function used is', pressure_func)
        # basic is defined in ICAR's atm_utilities
        pressure_data = np.zeros([self.nt,self.nz,self.ny,self.nx])
        # print(self.z_data[0,:,0,0])
        if model in ['basic', 'WeismanKlemp']:
            if pressure_func == 'calc_pressure_from_sea':
                pressure_data[:,:,:,:] = np.vectorize(calc_pressure_from_sea)(
                    self.sealevel_pressure,
                    self.z_data[:,:,:,:])
            elif pressure_func in ['calc_pressure_dz_iter',
                                   'calc_pressure_1m_iter']:
                pressure_data[:,0,:,:] = \
                    np.vectorize(pressure_func_d[pressure_func])(
                        self.sealevel_pressure,
                        0,
                        self.z_data[0,0,:,:])
                for z in range(1,self.nz):
                    pressure_data[:,z,:,:] = \
                        np.vectorize(pressure_func_d[pressure_func])(
                            pressure_data[0,z-1,:,:],
                            self.z_data[0,z-1,:,:],
                            self.z_data[0,z,:,:])
            self.pressure = xr.Variable(self.dims4d,
                                        pressure_data,
                                        {'long_name':'Pressure',
                                         'units':'Pa'})
        else:
             print("Error: ", pressure_model, " is not defined")
             exit()
        # print("NX NY", self.nx, self.ny)
        for t in range(self.nt):
            for i in range(self.nx):
                for j in range(self.ny):
                    if not np.array_equal(pressure_data[t,:,0,0],
                                          pressure_data[t,:,j,i]):
                        print("ERROR: PRESSURE DATA NOT EQUAL THROUGHOUT")
                        sys.exit()
        # print("--- PRESSURE DATA EQUAL THROUGHOUT ---")
        # print(pressure_data[0,:,0,0])
        # exit()
        del(pressure_data)


    def set_temperature(self, model='basic'):
        if (model in ['basic', 'WeismanKlemp']):
            # --TODO--
            # get better equation with temp and humidity
            # look at ICAR
            temp = np.zeros([self.nt, self.nz, self.ny, self.nx])
            temp = np.vectorize(calc_temp)(self.pressure.values,
                                           self.theta.values)
        else:
             print("Error: ", weather_model, " temperature is not defined")
             exit()
        self.temperature = xr.Variable(self.dims4d, temp,
                                 {'long_name':'Temperature',
                                  'units':"K"})


#---
# Lambda like functions used for np.vectorize
#---
# Weisman Klemp Theta equation
# z is elevation in meters
def calc_wk_theta(z):
    z_tr =  12000. # m
    theta_0 = 300. #
    theta_tr = 343. # K
    T_tr = 213. # K
    WK_C_p = 1000.0
    # WK_C_p = 1003.5
    # q_v0 = 11 # g kg^-1
    # q_v0 = 16 # g kg^-1
    # q_v0 = 14 # g kg^-1
    if z <= z_tr:
        theta = theta_0 + (theta_tr - theta_0) * (z / z_tr) ** (5./4)
    else:
        theta = theta_tr * math.exp((gravity / (WK_C_p * T_tr)) * (z - z_tr))
    return theta

# ---
# Functions taken from or based on functions from atm_utilities.f90
# ---
# Constancts from icar_constants.f90
gravity = 9.81
R_d = 287.058
C_p = 1003.5
Rd_over_Cp = R_d / C_p
P_0 = 100000

# p input pressure dz below, t temperature in layer between
# qv water vapor in layer between
def compute_p_offset(p, dz, t, qv):
    return p * exp( -dz / (Rd / gravity * ( t * ( 1 + 0.608 * qv ) )))

def calc_pressure_from_sea(sealevel_pressure, z):
    return sealevel_pressure * (1 - 2.25577E-5 * z)**5.25588

def calc_pressure_dz_iter(base_pressure, from_z, to_z):
    return base_pressure * (1 - 2.25577E-5 * (to_z - from_z))**5.25588

def calc_pressure_1m_iter(pressure, from_z, to_z):
    dz = 1
    for i in range(from_z,to_z,dz):
        pressure = pressure * (1 - 2.25577E-5 * dz)**5.25588
    return pressure

pressure_func_d = {
    'calc_pressure_dz_iter':calc_pressure_dz_iter,
    'calc_pressure_1m_iter':calc_pressure_1m_iter}

# theta * exner
def calc_temp(pressure, theta):
    return theta * (pressure / P_0)**Rd_over_Cp

# Modified from atm_utilities.f90
def sat_mr(temperature,pressure):

    e_s = np.zeros(temperature.shape)

    freezing = (temperature < 273.16)
    a = 21.8745584
    b = 7.66
    e_s[freezing] = 610.78 * np.exp(a * (temperature[freezing] - 273.16) / (temperature[freezing] - b))

    a = 17.2693882
    b = 35.86
    freezing = not freezing
    e_s[freezing] = 610.78 * np.exp(a * (temperature[freezing] - 273.16) / (temperature[freezing] - b))

    # not quite sure what this is needed for, maybe very low pressure rounding errors?
    high_es = e_s > pressure
    e_s[high_es] = pressure[high_es] * 0.9999

    sat_mr_val = 0.6219907 * e_s / (pressure - e_s)

    return sat_mr_val

def calc_exner(pressure):
    return (pressure / p_0) ** Rd_over_Cp
