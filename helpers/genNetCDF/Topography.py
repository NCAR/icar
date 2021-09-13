from netCDF4 import Dataset
from genNetCDF import fixType
import numpy as np
import math
from sys import exit

class Topography:
    nx = 2
    ny = 2
    nz = 10

    dimensions = {
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

    attributes = {
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


    variableNames = ["XLAT_M", "XLONG_M"]



    def __init__(self, nx=2, ny=2, nz=10, n_hills=1.0):
        self.nx = nx
        self.nz = nz
        self.ny = ny
        self.n_hills = n_hills
        self.topography_f = Dataset("topography.nc", "w", format="NETCDF4")

    def close(self):
        self.topography_f.close()


    def getDimension(self):
        return self.dimensions

    def getAttributes(self):
        return self.attributes

    def getVariableNames(self):
        return self.variableNames

    def createDimensions(self):
        for name, val in self.dimensions.items():
            self.topography_f.createDimension(name,val)

    def createAttributes(self):
        for name, val in self.attributes.items():
            fixed_val = fixType(val)
            self.topography_f.setncattr(name, fixed_val)

    def createVariables(self):
        # f = self.getFilename(filetype)
        pass
