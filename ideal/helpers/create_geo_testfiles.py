#!/usr/bin/env python
import numpy as np
import Nio

lofile=Nio.open_file("lo.nc","w")
lofile.create_dimension("x",10)
lofile.create_dimension("y",10)
lofile.create_dimension("z",10)

lat=np.arange(10)[:,np.newaxis].repeat(10,axis=1)
lon=np.arange(10)[np.newaxis,:].repeat(10,axis=0)
lofile.create_variable("XLAT","f",("y","x"))
lofile.create_variable("XLONG","f",("y","x"))

lofile.variables["XLAT"][:]=lat.astype("f")
lofile.variables["XLONG"][:]=lon.astype("f")

lofile.close()

hifile=Nio.open_file("hi.nc","w")
hifile.create_dimension("x",10)
hifile.create_dimension("y",10)
hifile.create_dimension("z",10)

lat=np.arange(10)[:,np.newaxis].repeat(10,axis=1)/3.0+3.01
lon=np.arange(10)[np.newaxis,:].repeat(10,axis=0)/3.0+3.25
hifile.create_variable("XLAT","f",("y","x"))
hifile.create_variable("XLONG","f",("y","x"))

hifile.variables["XLAT"][:]=lat.astype("f")
hifile.variables["XLONG"][:]=lon.astype("f")

hifile.close()
