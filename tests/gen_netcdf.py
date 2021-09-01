from netCDF4 import Dataset
import numpy as np
import IdealTest as it


def createVariables(f):
    for var in it.ds_variable_names:
        print(var)

    #   This is the high-resolution input filename
    #   primary inputs from this file are lat, lon, and terrain, optionally
    #   soil and veg types

def main():
    boundary = "boundary"
    forcing = "forcing"

    test = it.IdealTest(nx=10, ny=10)

    test.createDimensions(forcing)
    test.createVariables(forcing)
    test.createAttributes(forcing)

    test.close()

if __name__ == "__main__":
    main()
