from netCDF4 import Dataset
import numpy as np
import math
from sys import exit

class IdealTest:
    # class attributes
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


    boundaryDimensions = {
        "Time": None,
        "DataStrLen": 19,
        "south_north": nx,
        "west_east": ny,
        # "south_north_stag": 350,
        # "west_east_stag": 200,
        # "land_cat": 21,
        # "soil_cat": 16,
        "month": 1,
        "num_urb_params": 132,
        "dust_erosion_dimension": 3
    }

    boundaryAttributes = {
        "TITLE": "OUTPUT FROM CONTINUOUS INTEGRATION TEST",
        "SIMULATION_START_DATE": "0000-00-00_00:00:00",
        "WEST-EAST_GRID_DIMENSION": nx,
        "SOUTH-NORTH_GRID_DIMENSION": ny,
        "BOTTOM-TOP_GRID_DIMENSION": 0,
        # "WEST-EAST_PATCH_START_UNSTAG": 1,
        # "WEST-EAST_PATCH_END_UNSTAG": 199,
        # "WEST-EAST_PATCH_START_STAG": 1,
        # "WEST-EAST_PATCH_END_STAG": 200,
        # "SOUTH-NORTH_PATCH_START_UNSTAG": 1,
        # "SOUTH-NORTH_PATCH_END_UNSTAG": 349,
        # "SOUTH-NORTH_PATCH_START_STAG": 1,
        # "SOUTH-NORTH_PATCH_END_STAG": 350,
        "GRIDTYPE": "C",
        "DX": 100.0,
        "DY": 100.0
    }


    boundaryVariableNames = ["XLAT_M", "XLONG_M"]


    forcingDimensions = {
        "time": None,
        "level": nz,
        "lat": nx,
        "lon": ny
    }

    forcingAttributes = {
        "history": "Mon Dec 01 00:00:00 2020: ncrcat 2000-01-01_00:00:00 2000-01-31_0:00:00 "
    }

    forcingVariableNames = ["u", "v", "theta", "p", "lat", "lon"] # others

    def __init__(self, nx=2, ny=2, nz=10, n_hills=1.0):
        self.boundary_f = Dataset("boundary.nc", "w", format="NETCDF4")
        self.forcing_f = Dataset("forcing.nc", "w", format="NETCDF4")
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.n_hills = 1.0
        self.z_interface  = np.zeros((nx,nz,ny))
        self.z = np.full((nx,nz,ny), self.dz_value)
        self.genHill()
        self.dz_mass = np.full((nx,nz,ny), self.dz_value)
        self.dz_mass[:,0,:] = self.dz_mass[:,0,:]/2
        self.pressure = np.zeros((nx,nz,ny))


    def close(self):
        self.boundary_f.close()
        self.forcing_f.close()


    def getFilename(self, filetype):
        if filetype == "boundary":
            return self.boundary_f
        elif filetype == "forcing":
            return self.forcing_f
        else:
            print("ERROR")
            sys.exit()

    def getDimension(self, filetype):
        if filetype == "boundary":
            return self.boundaryDimensions
        elif filetype == "forcing":
            return self.forcingDimensions

    def getAttributes(self, filetype):
        if filetype == "boundary":
            return self.boundaryAttributes
        elif filetype == "forcing":
            return self.forcingAttributes

    def getVariableNames(self, filetype):
        if filetype == "boundary":
            return self.boundaryVariableNames
        elif filetype == "forcing":
            return self.forcingVariableNames

    def fixType(self, val):
        val_t = type(val)
        if val_t == int:
            return np.int32(val)
        elif val_t == float:
            return np.single(val)
        else:
            return val

# this%pressure(:,kms,:)    = pressure_at_elevation(sealevel_pressure, this%z(:,kms,:))

# do i=kms+1,kme
#   this%z(:,i,:)           = this%z(:,i-1,:)           + this%dz_mass(:,i,:)
#   this%z_interface(:,i,:) = this%z_interface(:,i-1,:) + this%dz_interface(:,i,:)
#   this%pressure(:,i,:)    = pressure_at_elevation(sealevel_pressure, this%z(:,i,:))
# end do




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
                self.z_interface[i-1,k,j-1] = self.surface_z + sine_curve * self.hill_height
        self.z[:,k,:] = self.z_interface[:,k,:] + self.dz_value/2

    def forcingVariables(name):
        print(name)


    def createDimensions(self, filetype):
        f = self.getFilename(filetype)
        dimensions = self.getDimension(filetype)

        for name, val in dimensions.items():
            f.createDimension(name,val)


    def createAttributes(self, filetype):
        f = self.getFilename(filetype)
        attributes = self.getAttributes(filetype)
        for name, val in attributes.items():
            fixed_val = self.fixType(val)
            f.setncattr(name, fixed_val)

    def createVariables(self, filetype):
        f = self.getFilename(filetype)



        # double time(time) ;
        #         time:long_name = "time" ;
        #         time:units = "days since 1900-01-01" ;
        #         time:calendar = "gregorian" ;
 # time = 36524, 36524.25, 36524.5, 36524.75, 36525, 36525.25, 36525.5,
 #    36525.75, 36526, 36526.25, 36526.5, 36526.75, 36527, 36527.25, 36527.5,
 #    ...  36554.75





    # def XLAT_M():
    #     print("XLAT_M")
        # float XLAT_M(Time, south_north, west_east) ;
        #         XLAT_M:FieldType = 104 ;
        #         XLAT_M:MemoryOrder = "XY " ;
        #         XLAT_M:units = "degrees latitude" ;
        #         XLAT_M:description = "Latitude on mass grid" ;
        #         XLAT_M:stagger = "M" ;
        #         XLAT_M:sr_x = 1 ;
        #         XLAT_M:sr_y = 1 ;
        # float XLONG_M(Time, south_north, west_east) ;
        #         XLONG_M:FieldType = 104 ;
        #         XLONG_M:MemoryOrder = "XY " ;
        #         XLONG_M:units = "degrees longitude" ;
        #         XLONG_M:description = "Longitude on mass grid" ;
        #         XLONG_M:stagger = "M" ;
        #         XLONG_M:sr_x = 1 ;
        #         XLONG_M:sr_y = 1 ;
