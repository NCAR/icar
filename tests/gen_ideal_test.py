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
nz = 80
nx = ny = 100
dx=dy=1000
dz_levels= [50., 75., 125., 200., 300., 400.] + [500.] * 50 
# dz_levels= [20., 25., 30., 40., 50., 75., 100., 150., 200., 250.,] + [300.] * 60 
# dz_levels=[20.0, 20.07, 20.32, 20.739999999999995, 21.38000000000001, 22.22, 23.269999999999996, 24.539999999999992, 26.03, 27.76000000000002, 29.7098, 31.89009999999996, 34.310100000000034, 36.96999999999997, 39.879800000000046, 43.02009999999996, 46.40989999999999, 50.05000000000001, 53.950199999999995, 58.08990000000006, 62.47989999999993, 67.14020000000005, 72.04999999999995, 77.21000000000004, 82.6499, 88.3401, 94.28999999999996, 100.50999999999999, 106.99000000000001, 113.75, 120.75999999999999, 128.0598, 135.61020000000008, 143.44989999999984, 151.55009999999993, 159.9299000000001, 168.58010000000013, 177.51999999999998, 186.7199999999998, 196.19990000000007, 205.9602, 216.01000000000022, 226.32989999999972, 236.9399000000003, 247.8171000000002, 259.0, 270.4531999999999, 282.1795999999995, 294.21880000000056, 306.5194999999994, 319.1289000000006, 332.0117999999993, 345.1796000000004, 358.64850000000024, 372.41010000000006, 386.4531999999999, 400.78899999999976, 415.4180000000006, 430.33979999999974, 445.5625, 461.0702999999994, 476.8672000000006, 492.98049999999967, 509.3711000000003, 526.0702999999994, 543.0586000000003, 560.3397999999997, 577.9415000000008, 595.8202999999994, 614.0077999999994, 632.5, 651.2812000000013, 670.3710999999985, 689.7695999999996, 709.4609000000019, 729.4491999999991, 749.7616999999991, 770.3594000000012, 791.2694999999985, 812.4805000000015]
decay_rate_L_topo = 2.0
decay_rate_S_topo = 5.0

# hill values currently do amazing things :)
hill_height       = 3000.0   # height of the ideal hill(s) [m]
n_hills           = 5.0      # number of hills across the domain

# perform an advection test, like the one in Schär 2002. Should have mp=pbl=0, adv=1 or 2. This will overwrite any u/v specified below. 
Schaer_test=True  
if Schaer_test==True:
    dx=dy=1000
    nz=50; nx =ny=300
    dz_levels = [500]*nz
    decay_rate_L_topo = 1.6667  # s1 = 15 km
    decay_rate_S_topo = 10.0    # s2 = 2.5 km



# ---- Forcing specs  -----
nt=3 # nr of timesteps (hours)
# from ideal test
dz_value          = 500.0    # thickness of each (Forcing?) model gridcell   [m]  
# u field can be a constant, or an array of size nz. When Schaer_test is chosen, this get overwritten. 
u_test_val = np.array([0., 0., 0., 0., 0., 0. , 0. ,2.] + [5.] *35)   # u field in z-direction
# u_test_val = 5.0
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
                    model_comment = 'flat_z_height=-10',
                    flat_z_height=-10, 
                    sleve= ".True.",
                    terrain_smooth_windowsize = 4,   
                    terrain_smooth_cycles = 2,
                    decay_rate_L_topo = decay_rate_L_topo,
                    decay_rate_S_topo = decay_rate_S_topo,
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
                  lat0 = 39.5,lon0 = -105,
                  Schaer_test=Schaer_test
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
