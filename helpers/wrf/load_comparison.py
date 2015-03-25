#!/usr/bin/env python
import glob
import numpy as np
import mygis

from bunch import Bunch

# icar[:,4:,1:])*10-wrf[:,:-6,:-4]

def load_wrf(filename, preciponly=False):
    """docstring for load_wrf"""
    precip=mygis.read_nc(filename,"RAINNC").data
    try:
        precip+=mygis.read_nc(filename,"I_RAINNC").data*100 # try to add bucket data
    except KeyError:
        pass
    if preciponly:
        return Bunch(precip=precip,hgt=None)
    
    t=mygis.read_nc(filename,"T").data+300
    qv=mygis.read_nc(filename,"QVAPOR").data
    p=mygis.read_nc(filename,"P").data
    p+=mygis.read_nc(filename,"PB").data
    
    u=mygis.read_nc(filename,"U").data
    v=mygis.read_nc(filename,"V").data
    w=mygis.read_nc(filename,"W").data
    
    z=mygis.read_nc(filename,"PH").data
    z+=mygis.read_nc(filename,"PHB").data
    z/=9.8
    hgt=mygis.read_nc(filename,"HGT").data
    
    hgt=hgt[0,:-6,:-4]
    p=p[:,:,:-6,:-4]
    t=t[:,:,:-6,:-4]
    qv=qv[:,:,:-6,:-4]
    u=u[:,:,:-6,:-4]
    v=v[:,:,:-6,:-4]
    w=w[:,:,:-6,:-4]
    precip=precip[:,:-6,:-4]
    return Bunch(w=w,z=z,hgt=hgt,u=u,v=v,precip=precip,p=p,qv=qv,t=t)

def load_icar(filename,preciponly=False):
    """docstring for load_wrf"""
    
    precip=mygis.read_nc(filename,"rain").data
    if preciponly:
        return Bunch(precip=precip)
        
    u=mygis.read_nc(filename,"u").data
    v=mygis.read_nc(filename,"v").data
    w=mygis.read_nc(filename,"w").data
    
    z=mygis.read_nc(filename,"z").data
    dz=np.zeros(z.shape)
    
    dx=4000.0
    dz[:,0,:]=z[:,0,:]*2
    for i in range(1,z.shape[1]):
        dz[:,i,:] = 2*(z[:,i,:]-z[:,i-1,:]) - dz[:,i-1,:]
        w[:,:,i,:] *= dz[np.newaxis,:,i,:] / dx
    
    dzdx=np.diff(z,axis=2)
    dzdy=np.diff(z,axis=0)
    
    w_u=u[:,:,:,:-2]*dzdx[np.newaxis,:,:,:]/dx
    w_u=(w_u[:,:,:,1:]+w_u[:,:,:,:-1])/2.0

    w_v=v[:,:-2,:,:]*dzdy[np.newaxis,:,:,:]/dx
    w_v=(w_v[:,1:,:,:]+w_v[:,:-1,:,:])/2.0
    
    w[:,1:-1,:,1:-1]+=w_u[:,1:-1,:,:]+w_v[:,:,:,1:-1]
    

    qv=mygis.read_nc(filename,"qv").data
    t=mygis.read_nc(filename,"th").data
    p=mygis.read_nc(filename,"p").data

    w=w[:,4:,:,1:]
    v=v[:,4:,:,1:]
    u=u[:,4:,:,1:]
    z=z[4:,:,1:]
    precip=precip[:,4:,1:]
    t=t[:,4:,:,1:]
    qv=qv[:,4:,:,1:]
    p=p[:,4:,:,1:]
    
    return Bunch(z=z,w=w,u=u,v=v,qv=qv,t=t,p=p,precip=precip)
    
def load_multi(fnames,load_func=load_icar):
    """docstring for load_multi"""
    if type(fnames)!=list:
        fnames=glob.glob(fnames)
    fnames.sort()
    
    data=[]
    for f in fnames:
        print(f)
        data.append(load_func(f))
    
    master_data=Bunch()
    print(data[0].keys())
    for k in data[0].keys():
        print(k)
        if (k=="z") or (k=="hgt"):
            master_data[k]=data[0][k]
        else:
            if k=="precip":
                nt,ny,nx=data[0][k].shape
                master_data[k]=np.zeros((nt*len(data),ny,nx))
            else:
                nt,ny,nz,nx=data[0][k].shape
                master_data[k]=np.zeros((nt*len(data),ny,nz,nx))

            for i in range(len(data)):
                master_data[k][i*nt:(i+1)*nt]=data[i][k]
    
    return master_data
            
