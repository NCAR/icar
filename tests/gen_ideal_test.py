#!/usr/bin/env python3
from ast import Or
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
nx = 300
ny = 20
dx=dy=1000

dz_levels= [50., 75., 125., 200., 300., 400.] + [500.] * 50
decay_rate_L_topo = 2.0
decay_rate_S_topo = 5.0
terrain_smooth_windowsize = 4
terrain_smooth_cycles = 5

# Paramters to generate hi-res topography:
hill_height       = 3000.0   # height of the ideal hill(s) [m]
n_hills           = 5.0      # number of hills across the domain (ignored for Schaer test)

# perform an advection test, like the one in Schär 2002. Should have mp=pbl=0, adv=1 or 2. This will overwrite any Forcing u/v specified below.
# More details in 'A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models' Christoph Schär et al, 2002, Monthly Weather Review vol 130.
Schaer_test=True
if Schaer_test==True:
    print(" - - - Setting up an idealized advection test - - - \n")
    # Values below are as specified in Schaer et al 2002.
    dx = dy = 1000                # 1000 m
    nx = 300; ny = 20             # 300 m, ..m ?
    nz = 50                       # 50 m
    dz_levels = [500]*nz          # 500 m , model top at 25 km
    decay_rate_L_topo = 1.6667    # s1 = 15 km
    decay_rate_S_topo = 13.0      # s2 = 2.5 km
    hill_height       = 3000.0    # height of the ideal hill(s) [m]



# ---- Forcing specs  -----
nt_lo         = 4 # nr of timesteps (hours) - this is also how long the ideal sim will run.
nz_lo         = 51
nx_lo = 300; ny_lo = 20
dx_lo = dy_lo = 1000     # make sure dx_lo*nx_lo => dx*nx & dy_lo*ny_lo => dy*ny
dz_lo         = 500.0    # thickness of each (Forcing?) model gridcell   [m]

if dx_lo*nx_lo < dx*nx or dy_lo*ny_lo < dy*ny or dz_lo*nz_lo < np.sum(dz_levels) :
    print("\n   ERROR: Forcing domain smaller than hi-res domain. Incease forcing domain size \n")

# u field can be a constant, or an array of size nz. When Schaer_test is chosen, this get overwritten.
u_test_val = np.array([0., 0., 0., 0., 0., 0. , 0. ,2.] + [5.] *35)   # u field in z-direction
# u_test_val = 5.0
v_test_val = 0.0

# relative_humidity = 0.01
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
                    output_vars=['pressure','temperature', 'lon', 'lat', 'z', 'dz', 'dz_i', 'u', 'v', 'w', 'w_grid', 'qv', 'terrain' ],
                    dz_levels=dz_levels,
                    model_comment = 'flat_z_height=-10',
                    flat_z_height=-10,
                    sleve= ".True.",
                    terrain_smooth_windowsize = terrain_smooth_windowsize,
                    terrain_smooth_cycles = terrain_smooth_cycles ,
                    decay_rate_L_topo = decay_rate_L_topo,
                    decay_rate_S_topo = decay_rate_S_topo,
                    sleve_n = 1.35,
                    space_varying = ".True.",
                    dx = dx,                       # <-   affects advection speed!
                    phys_opt_mp = 0,
                    phys_opt_adv = 1,
                    phys_opt_wind = 3,
                    smooth_wind_distance = dx_lo,  # Very important - has effect on vertical speeds!
                    use_agl_height = True,         # !   Use height above ground level to interpolate the wind field instead of height above sea level.
                    agl_cap = 400,                 # !   Height at which we switch from AGL-interpolation to using ASL-interpolation
                    output_file = 'icar_out_',
                    qv_is_relative_humidity ='false',
                    output_interval = 1200,
                    end_date = '2020-12-01 02:00:00',
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
    forcing = fc.Forcing(nt=nt_lo, nz=nz_lo, nx=nx_lo+10, ny=ny_lo+10,
                         u_val=u_test_val,
                         v_val=v_test_val,
                         water_vapor_val=water_vapor_test_val,
                         dz_value=dz_lo,
                         dx=dx_lo, dy=dy_lo,
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
