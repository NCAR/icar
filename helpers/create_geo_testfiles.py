#!/usr/bin/env python
import numpy as np
import Nio

nx=100
ny=80
nz=5

lofile=Nio.open_file("lo.nc","w")
lofile.create_dimension("x",nx)
lofile.create_dimension("y",ny)
lofile.create_dimension("z",nz)

lat=np.arange(ny)[:,np.newaxis].repeat(nx,axis=1)
lon=np.arange(nx)[np.newaxis,:].repeat(ny,axis=0)
lofile.create_variable("XLAT","f",("y","x"))
lofile.create_variable("XLONG","f",("y","x"))

lofile.variables["XLAT"][:]=lat.astype("f")
lofile.variables["XLONG"][:]=lon.astype("f")

lofile.close()

hifile=Nio.open_file("hi.nc","w")
hifile.create_dimension("x",nx)
hifile.create_dimension("y",ny)
hifile.create_dimension("z",nz)

lat=np.arange(ny)[:,np.newaxis].repeat(nx,axis=1)/30.0+3.01
lon=np.arange(nx)[np.newaxis,:].repeat(ny,axis=0)/30.0+3.25
hifile.create_variable("XLAT","f",("y","x"))
hifile.create_variable("XLONG","f",("y","x"))

hifile.variables["XLAT"][:]=lat.astype("f")
hifile.variables["XLONG"][:]=lon.astype("f")

hifile.close()
