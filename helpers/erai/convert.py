import numpy as np
from bunch import Bunch

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

# icar_atm_var=["u","v","gph","t","qv","ln_p_sfc","cloud","ice","sigma"]
def convert_atm(data):
    output_data=Bunch()
    output_data.u=data.u[::-1,::-1,:]                        # m/s
    output_data.v=data.v[::-1,::-1,:]                        # m/s
    output_data.hgt=data.gph[::-1,:]/g                       # (m^2/s^2) / (m/s^2) = m

    # calculate pressure in Pa from ln(sfc_press) and hybrid sigma coordinates
    output_data.p=np.zeros((data.u.shape))
    ps=np.exp(data.ln_p_sfc[::-1,:])                         # Pa
    for i in range(len(data.sigma_a)):
        # see http://rda.ucar.edu/datasets/ds627.0/docs/Eta_coordinate/
        # notes on http://aaron.boone.free.fr/aspdoc/node7.html might help...
        output_data.p[len(data.sigma_a)-i-1,:,:]=(data.sigma_a[i]*data.P0+data.sigma_b[i]*ps)

    psl=ps/((1 - 2.25577E-5*output_data.hgt)**5.25588)
    output_data.z=np.zeros(output_data.p.shape)
    for i in range(output_data.p.shape[0]):
        output_data.z[i,...]=((output_data.p[i,...]/psl)**(1/5.25588)-1) / (-2.25577E-5)

    pii=(100000.0/output_data.p)**(R/cp)
    output_data.t=data.t[::-1,::-1,:]*pii                    #K (convertred to potential temperature)
    
    output_data.qv=data.qv[::-1,::-1,:]                      # kg/kg
    output_data.cloud=data.cloud[::-1,::-1,:]                # kg/kg
    output_data.ice=data.ice[::-1,::-1,:]                    # kg/kg
    
    return output_data

# icar_sfc_var=["sensible_heat","latent_heat","hgt_98","PBL_height"]
def convert_sfc(data):
    dt= -3.*60.*60.
    output_data=Bunch()
    output_data.sensible_heat=data.sensible_heat[::-1,:]/dt  # W/m^2
    output_data.latent_heat=data.latent_heat[::-1,:]/dt      # W/m^2
    output_data.sfc_hgt=data.hgt_98[::-1,:]/g                # (m^2/s^2) / (m/s^2) = m
    output_data.PBL_height=data.PBL_height[::-1,:]           # m
    return output_data

def era2icar(data):
    output_data=Bunch()
    atm=convert_atm(data.atm)
    sfc=convert_sfc(data.sfc)
    
    for k in atm.keys():
        output_data[k]=atm[k]
    for k in sfc.keys():
        output_data[k]=sfc[k]
    
    return output_data