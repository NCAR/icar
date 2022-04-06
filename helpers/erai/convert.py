import numpy as np
from bunch import Bunch

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2
global last_longwave
last_longwave = None
global last_shortwave
last_shortwave = None

# icar_atm_var=["u","v","gph","t","qv","ln_p_sfc","cloud","ice","sigma"]
def convert_atm(data):
    output_data=Bunch()
    output_data.u   = data.u[np.newaxis, ::-1,::-1,:]           # m/s
    output_data.v   = data.v[np.newaxis, ::-1,::-1,:]           # m/s
    output_data.hgt = data.gph[::-1,:]/g                        # (m^2/s^2) / (m/s^2) = m

    # calculate pressure in Pa from ln(sfc_press) and hybrid sigma coordinates
    output_data.p = np.zeros((output_data.u.shape))
    ps = np.exp(data.ln_p_sfc[::-1,:])                          # Pa
    for i in range(len(data.sigma_a)):
        # see http://rda.ucar.edu/datasets/ds627.0/docs/Eta_coordinate/
        # notes on http://aaron.boone.free.fr/aspdoc/node7.html might help...
        output_data.p[0,len(data.sigma_a)-i-1,:,:]=(data.sigma_a[i]*data.P0+data.sigma_b[i]*ps)

    psl=ps/((1 - 2.25577E-5*output_data.hgt)**5.25588)
    output_data.z=np.zeros(output_data.p.shape)
    for i in range(output_data.p.shape[1]):
        output_data.z[0,i,...]=((output_data.p[0,i,...]/psl)**(1/5.25588)-1) / (-2.25577E-5)

    pii=(100000.0/output_data.p)**(R/cp)
    output_data.t=data.t[np.newaxis,::-1,::-1,:]*pii                    # K (converted to potential temperature)

    output_data.qv    = data.qv[np.newaxis,::-1,::-1,:]                 # kg/kg
    output_data.cloud = data.cloud[np.newaxis,::-1,::-1,:]              # kg/kg
    output_data.ice   = data.ice[np.newaxis,::-1,::-1,:]                # kg/kg

    return output_data

def bfill(arr):
    ''' from https://stackoverflow.com/questions/41190852/most-efficient-way-to-forward-fill-nan-values-in-numpy-array
    '''
    mask = np.isnan(arr)
    idx = np.where(~mask, np.arange(mask.shape[1]), mask.shape[1] - 1)
    idx = np.minimum.accumulate(idx[:, ::-1], axis=1)[:, ::-1]
    out = arr[np.arange(idx.shape[0])[:,None], idx]
    return out


def numpy_fill(arr):
    '''modified from Solution provided by Divakar.
    from https://stackoverflow.com/questions/41190852/most-efficient-way-to-forward-fill-nan-values-in-numpy-array
    '''
    for i in range(arr.shape[0]):
        mask = np.isnan(arr[i])
        idx = np.where(~mask,np.arange(mask.shape[1]),0)
        np.maximum.accumulate(idx,axis=1, out=idx)
        out = arr[i,np.arange(idx.shape[0])[:,None], idx]
        arr[i] = bfill(out) # in case there are still missing values on the left side

    return arr

# icar_sfc_var=["sensible_heat","latent_heat","hgt_98","PBL_height"]
def convert_sfc(data):
    global last_longwave
    global last_shortwave
    dt = 3.0 * 60.0 * 60.0
    output_data=Bunch()
    output_data.sensible_heat   = data.sensible_heat[np.newaxis,::-1,:]/dt  # W/m^2
    output_data.latent_heat     = data.latent_heat[np.newaxis,::-1,:]/dt    # W/m^2
    output_data.sfc_hgt         = data.hgt_98[::-1,:]/g                     # (m^2/s^2) / (m/s^2) = m
    output_data.PBL_height      = data.PBL_height[np.newaxis,::-1,:]        # m
    output_data.tskin           = data.tskin[np.newaxis,::-1,:]             # K
    output_data.sw              = data.sw[np.newaxis,::-1,:] / dt   # convert from Joules to W /m^2
    output_data.lw              = data.lw[np.newaxis,::-1,:] / dt   # convert from Joules to W /m^2
    output_data.cp              = data.cp[np.newaxis,::-1,:] * 1000 # convert m to mm

    output_data.landmask = data.landmask[np.newaxis,::-1,:]
    # landval = data.tskin[np.argmax(data.landmask)] # ~273.15, alternatively, tskin[landmask>0.99].mean()
    #  above seems to always create an array, and sometimes with very different values in it ... e.g. >300...
    landval = 273.16
    output_data["sst"]           = (data.tskin[np.newaxis,::-1,:] - (output_data.landmask * landval)) / (1 - output_data.landmask)
    output_data["sst"][output_data.landmask>0.25] =np.nan
    output_data["sst"] = numpy_fill(output_data["sst"])
    # this is now handled in io so it can just use the last value in the file, much simple
    #  ... though in some ways what is below is better as it integrates over a longer time period
    # if last_longwave==None:
    #     last_shortwave = np.zeros(output_data.sw.shape) + output_data.sw
    #     last_longwave = np.zeros(output_data.lw.shape) + output_data.lw
    #     # initial values for the day were only integrated over a 3hr period
    #     output_data.sw*=2
    #     output_data.lw*=2
    # else:
    #     if output_data.lw.max() > last_longwave.max():
    #         temp = np.zeros(output_data.lw.shape) + output_data.lw
    #         output_data.lw -= last_longwave
    #         last_longwave = temp
    #         temp = np.zeros(output_data.sw.shape) + output_data.sw
    #         output_data.sw -= last_shortwave
    #         last_shortwave = temp
    #     else:
    #         last_longwave[:] = output_data.lw
    #         last_shortwave[:] = output_data.sw
    #         # initial values for the day were only integrated over a 3hr period
    #         output_data.sw*=2
    #         output_data.lw*=2
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
