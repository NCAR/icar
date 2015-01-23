import numpy as np
from bunch import Bunch

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

def convert_atm(data):
    output_data=Bunch()
    output_data.u=data.u[:,::-1,:,:]                        # m/s
    output_data.v=data.v[:,::-1,:,:]                        # m/s
    output_data.p=data.p[:,::-1,:,:]                        # Pa
    output_data.z=data.z[:,::-1,:,:]                        # m

    output_data.dz=np.zeros(output_data.z.shape[1])
    
    output_data.dz[0]= 2 * np.mean(output_data.z[0,0,:,:])
    for i in range(1,output_data.z.shape[1]):
        output_data.dz[i]= 2 * np.mean(output_data.z[:,i,:,:]-output_data.z[:,i-1,:,:])-output_data.dz[i-1]
    output_data.dz[0]=output_data.dz[1]
    
    pii=(100000.0/output_data.p)**(R/cp)
    output_data.t=data.t[:,::-1,:,:]*pii                    # K (converted to potential temperature)
    
    output_data.qv=data.qv[:,::-1,:,:]
    output_data.cloud=np.zeros(data.qv.shape)
    output_data.ice=np.zeros(data.qv.shape)
    
    return output_data

def cesm2icar(data):
    output_data=Bunch()
    atm=convert_atm(data.atm)
    
    for k in atm.keys():
        output_data[k]=atm[k]
    for k in data.sfc.keys():
        output_data[k]=data.sfc[k]
    
    return output_data