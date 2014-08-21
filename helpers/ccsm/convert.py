import numpy as np
from bunch import Bunch

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

def convert_atm(data):
    output_data=Bunch()
    output_data.u=data.u                        # m/s
    output_data.v=data.v                        # m/s
    output_data.p=data.p                                     # Pa

    output_data.hgt=(1-(data.ps.mean(axis=0)/101600)**0.190263)/2.25577e-5
    output_data.hgt[output_data.hgt<0]=0
    z=(1-(data.p.mean(axis=0)/101600)**0.190263)/2.25577e-5 - output_data.hgt[np.newaxis,:,:]
    output_data.dz=np.zeros(z.shape[0])
    
    output_data.dz[0]= 2 * np.mean(z[0,:,:])
    for i in range(1,z.shape[0]):
        output_data.dz[i]= 2 * np.mean(z[i,:,:]-z[i-1,:,:])-output_data.dz[i-1]
    
    pii=(100000.0/output_data.p)**(R/cp)
    output_data.t=data.t*pii                    # K (converted to potential temperature)
    
    output_data.qv=data.qv
    output_data.cloud=np.zeros(data.qv.shape)
    output_data.ice=np.zeros(data.qv.shape)
    
    return output_data

def ccsm2icar(data):
    output_data=Bunch()
    atm=convert_atm(data.atm)
    
    for k in atm.keys():
        output_data[k]=atm[k]
    for k in data.sfc.keys():
        output_data[k]=data.sfc[k]
    
    return output_data