#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
from sys import exit, path
from os import getcwd
path.insert(0, getcwd()+'/../helpers/genNetCDF')
import Topography as tg
import Forcing as fc
import ICARoptions as opt

# ---------------------------------------
# ----- Settings For Generating Files -----
# ---------------------------------------
# choose dimensions
nz = 40
nx = ny = 100
dx=dy=1000
# from ideal test
dz_value          = 200.0    # thickness of each (Forcing?) model gridcell   [m]
# hill values currently do amazing things :)
hill_height       = 2000.0   # height of the ideal hill(s) [m]
n_hills           = 5.0      # number of hills across the domain

# relative_humidity = 0.01
u_test_val = v_test_val = w_test_val = 5.0 #0.0
water_vapor_test_val = 0.000
mixing_ratio = 0.001 # water vapor # not if constant
qv_val = mixing_ratio

# --- choose function for creating pressure ---
pressure_func = 'calc_pressure_from_sea'
# pressure_func = 'calc_pressure_dz_iter'
# pressure_func = 'calc_pressure_1m_iter'
# --- choose weather model ---
# weather_model = 'basic'
weather_model = 'WeismanKlemp'


def main():
    # ICAR Options generate the ICAR namelist
    opt.ICARoptions(nz=nz,
                    output_vars=['pressure','temperature', 'lon', 'lat', 'z', 'dz_i', 'u', 'v', 'w', 'w_grid', ],
                    sleve= ".True.",
                    space_varying = ".True.")
    print("Generated icar_options.nml")

    tg.Topography(nz, nx, ny, n_hills=n_hills, hill_height=hill_height, dx=dx, dy=dy)
    print("Generated init.nc")

    # double check all passed variable get used
    forcing = fc.Forcing(nz, nx, ny,
                         u_val=u_test_val,
                         v_val=v_test_val,
                         w_val=w_test_val,
                         water_vapor_val=water_vapor_test_val,
                         dz_value=dz_value,
                         qv_val=qv_val,
                         weather_model=weather_model,
                         pressure_func=pressure_func)
    print("Generated forcing.nc")

if __name__ == "__main__":
    main()
