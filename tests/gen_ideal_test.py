from netCDF4 import Dataset
import numpy as np
import math
import Boundary as bd
import Forcing as fc
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


    def __init__(self, nx=2, ny=2, nz=10, n_hills=1.0):
        self.boundary = bd.Boundary(nx,ny,nz,n_hills)
        self.forcing = fc.Forcing(nx,ny,nz,n_hills)
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
        self.boundary.close()
        self.forcing.close()


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



def main():

    test = IdealTest(nx=10, ny=10)
    test.boundary.createDimensions()
    test.boundary.createVariables()
    test.boundary.createAttributes()

    test.forcing.createDimensions()
    test.forcing.createVariables()
    test.forcing.createAttributes()

    test.close()

if __name__ == "__main__":
    main()


# this%pressure(:,kms,:)    = pressure_at_elevation(sealevel_pressure, this%z(:,kms,:))

# do i=kms+1,kme
#   this%z(:,i,:)           = this%z(:,i-1,:)           + this%dz_mass(:,i,:)
#   this%z_interface(:,i,:) = this%z_interface(:,i-1,:) + this%dz_interface(:,i,:)
#   this%pressure(:,i,:)    = pressure_at_elevation(sealevel_pressure, this%z(:,i,:))
# end do




    # def forcingVariables(name):
    #     print(name)


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
