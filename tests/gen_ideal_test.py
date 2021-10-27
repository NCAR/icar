from netCDF4 import Dataset
import numpy as np
import math
from sys import exit, path
from os import getcwd
path.insert(0, getcwd()+'/../helpers/genNetCDF')
import Topography as tg
import Forcing as fc
import ICARoptions as opt

# Python program generates an ideal case
class IdealTest:
    # from ideal test
    sealevel_pressure = 100000.0 # pressure at sea level [Pa]
    dz_value          = 500.0    # thickness of each model gridcell   [m]
    # hill values currently do nothing
    hill_height       = 1000.0   # height of the ideal hill(s) [m]
    n_hills           = 1.0      # number of hills across the domain

    def __init__(self, nz=10, nx=2, ny=2, n_hills=1.0):
        rh = 0.9
        u_test_val = 0.5
        v_test_val = 0.5
        w_test_val = 0.0
        water_vapor_test_val = 0.001
        theta_test_val = 300.0

        self.forcing = fc.Forcing(nz, nx, ny, self.sealevel_pressure,
                                  rh, u_test_val, v_test_val, w_test_val,
                                  water_vapor_test_val, theta_test_val,
                                  dz_value=self.dz_value)

        self.init = tg.Topography(nz, nx, ny)


def main():
    # ICAR Options generate the ICAR namelist
    options = opt.ICARoptions()
    test = IdealTest(nz=40, nx=40, ny=40, n_hills=1.0)

if __name__ == "__main__":
    main()
