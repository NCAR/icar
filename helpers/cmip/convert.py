import numpy as np
import units
from bunch import Bunch

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

def convert_atm(data,sfc):
    output_data=Bunch()
    #                  [time,z,ns,ew]
    output_data.u  = data.u                   # m/s
    output_data.v  = data.v                   # m/s
    output_data.p  = data.p                   # Pa
    output_data.qv = data.qv                  # kg/kg
    
    pii = (100000.0 / output_data.p)**(R / cp)
    output_data.t = data.t * pii              # K (converted to potential temperature)
    
    
    if "z" in data.keys():
        output_data.z=data.z                        # m
    else:
        print(data.ps.shape,
            sfc.hgt.shape,
            data.t.shape)
        print(data.p[0,0,...].mean(),        data.p[0,-1,...].mean())

        data.slp=units.calc_slp(data.ps,sfc.hgt[np.newaxis,...],ts=data.t[:,0,...],mr=output_data.qv[:,0,...],method=2)
        output_data.z=np.zeros(data.t.shape)
        for z_time in range(data.t.shape[0]):
            output_data.z[z_time,...]=units.calc_z(data.slp[z_time],
                                                  output_data.p[z_time],
                                                  t=output_data.t[z_time],
                                                  mr=output_data.qv[z_time])
    
    # now calculate layer thicknesses
    output_data.dz=np.zeros(output_data.z.shape)
    
    output_data.dz[:,0,:,:]= 2 * (output_data.z[:,0,:,:]-sfc.hgt[np.newaxis,:,:])
    for i in range(1,output_data.z.shape[1]):
        output_data.dz[:,i,:,:]= 2 * np.mean(output_data.z[:,i,:,:]-output_data.z[:,i-1,:,:])-output_data.dz[:,i-1]
    output_data.dz[0]=output_data.dz[1]
    
    output_data.cloud= np.zeros(data.qv.shape)
    output_data.ice  = output_data.cloud
    
    
    
    
    return output_data

def cmip2icar(data):
    output_data=Bunch()
    atm=convert_atm(data.atm,data.sfc)
    
    for k in atm.keys():
        output_data[k]=atm[k]
    for k in data.sfc.keys():
        output_data[k]=data.sfc[k]
    
    return output_data