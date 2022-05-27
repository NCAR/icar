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
# choose dimensions of hi-res grid:
nz = 40
nx = ny = 100
dx=dy=1000
dz_levels= [50., 75., 125., 200., 300., 400.] + [500.] * 50 

# hill values currently do amazing things :)
hill_height       = 3000.0   # height of the ideal hill(s) [m]
n_hills           = 5.0      # number of hills across the domain

# perform an advection test, like the one in Schär 2002. Should have mp=pbl=0, adv=1 or 2. This will overwrite any u/v specified below. 
Schaer_test=True  

# ---- Forcing specs  -----
nt=12
# from ideal test
dz_value          = 500.0    # thickness of each (Forcing?) model gridcell   [m]  
# u field can be a constant, or an array of size nz
u_test_val = np.array([0., 0., 0., 0., 0., 0. , 0. ,2.] + [5.] *35)   # u field in z-direction
v_test_val = 0.0

# relative_humidity = 0.01
water_vapor_test_val = 0.000
mixing_ratio = 0.1 # water vapor # not if constant
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
                    output_vars=['pressure','temperature', 'lon', 'lat', 'z', 'dz_i', 'u', 'v', 'w', 'w_grid', 'qv' ],
                    dz_levels=dz_levels,
                    model_comment = 'flat_z_height=0',
                    flat_z_height=0, 
                    sleve= ".True.",
                    terrain_smooth_windowsize = 4,   
                    terrain_smooth_cycles = 5,
                    decay_rate_L_topo = 2.0,
                    decay_rate_S_topo = 5.0,
                    sleve_n = 1.35,
                    space_varying = ".True.",
                    phys_opt_mp = 0,
                    phys_opt_wind = 3,
                    use_agl_height = True,  # !   Use height above ground level to interpolate the wind field instead of height above sea level.
                    agl_cap = 400,        # !   Height at which we switch from AGL-interpolation to using ASL-interpolation
)
    print("Generated icar_options.nml")

    tg.Topography(nz, nx, ny, 
                  n_hills=n_hills, hill_height=hill_height, 
                  dx=dx, dy=dy,
                  lat0 = 39.5,lon0 = -105
                  )
    print("Generated init.nc")

    # double check all passed variable get used
    forcing = fc.Forcing(nt=nt, nz=nz, nx=nx+10, ny=ny+10,    # should/could become separate from hi-res nx/y/z 
                         u_val=u_test_val,
                         v_val=v_test_val,
                        #  w_val=w_test_val,
                         water_vapor_val=water_vapor_test_val,
                         dz_value=dz_value,
                         dx=dx, dy=dy,
                         qv_val=qv_val,
                         weather_model=weather_model,
                         pressure_func=pressure_func,
                         hill_height=hill_height,
                         lat0 = 39.5,lon0 = -105,
                         Schaer_test=Schaer_test  
                        )
    print("Generated forcing.nc")

if __name__ == "__main__":
    main()
