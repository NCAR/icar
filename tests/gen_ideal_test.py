from netCDF4 import Dataset
import numpy as np
import math
import Topography as tg
import Forcing as fc
import ICARoptions as opt
from sys import exit

class IdealTest:
    # class attributes
    nz = 10
    nx = 2
    ny = 2
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


    def __init__(self, nz=10, nx=2, ny=2, n_hills=1.0):
        self.nz = nz
        self.nx = nx
        self.ny = ny
        self.n_hills = 1.0
        self.z_interface  = np.zeros((nz,nx,ny))
        self.z = np.full((nz,nx,ny), self.dz_value)
        self.genHill()
        self.dz_mass = np.full((nz,nx,ny), self.dz_value)
        self.dz_mass[0,:,:] = self.dz_mass[0,:,:]/2
        self.pressure = np.zeros((nz,nx,ny))

        self.topography = tg.Topography(nz,nx,ny,n_hills)

        sealevel_pressure = 100000.0
        rh = 0.9
        u_test_val = 0.5
        v_test_val = 0.5
        w_test_val = 0.0
        water_vapor_test_val = 0.001
        theta_test_val = 300.0
 # - forcing can be completely flat, potentialtempt,  pressure, mixing ratio,
 #   vertical coaringses, potential tempt constant, RH etc,

        self.forcing = fc.Forcing(nz, nx, ny, sealevel_pressure,
                                  rh, u_test_val, v_test_val, w_test_val,
                                  water_vapor_test_val, theta_test_val)



    def close(self):
        self.topography.close()
        # self.forcing.close()


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


def main():
    options = opt.ICARoptions()

    test = IdealTest(nz=40, nx=4, ny=4, n_hills=1.0)
    # test.topography.createDimensions()
    # test.topography.createVariables()
    # test.topography.createAttributes()

    test.close()

if __name__ == "__main__":
    main()



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
