from __future__ import print_function
import numpy as np
import time

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

def get_real_t(atm,potentialT):
    return potentialT/((100000.0/atm.p)**(R/cp))
def get_potential_t(atm):
    return atm.ta*(100000.0/atm.p)**(R/cp)

def swim2d(domain,weather,swim,options):
    # modified to work with new calling sequence.  
    # Eventually this whole routine should get cleaned up too
    # currently just recasting arrays to order="F" is taking ~1sec, 
    # with simplified microphysics this is becoming important
    t0=time.time()
    physics=options.physics
    oldatm=weather.old
    atm=weather.new
    z=domain.topo
    dx=options.wrfres
    timestep=options.timestep
    
    N=z.shape
    if timestep==None:
        timestep=3*60.0*60.0 # seconds
    # dTdt=np.array(dTdt.transpose([1,0,2]),order="F")
    th=np.array(oldatm.th.transpose([1,0,2]), order='F',dtype=np.float32) # theta = potential temperature
    # all of these should be mixing ratios kg/kg-air (not kg/kg-dry-air)
    qv=np.array(oldatm.qv.transpose([1,0,2]),order='F',dtype=np.float32) # vapor pressure
    qc=np.array(oldatm.qc.transpose([1,0,2]),order='F',dtype=np.float32) # cloud water
    qr=np.array(oldatm.qr.transpose([1,0,2]),order='F',dtype=np.float32) # rain in air
    nr=np.array(oldatm.nr.transpose([1,0,2]),order='F',dtype=np.float32) # number concentration of rain drops
    qs=np.array(oldatm.qs.transpose([1,0,2]),order='F',dtype=np.float32) # snow in air
    qi=np.array(oldatm.qi.transpose([1,0,2]),order='F',dtype=np.float32) # cloud ice
    ni=np.array(oldatm.ni.transpose([1,0,2]),order='F',dtype=np.float32) # number concentration of ice crystals
    qg=np.array(oldatm.qg.transpose([1,0,2]),order='F',dtype=np.float32) # graupel in air
    pptrain=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptsnow=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptgraul=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptrainV=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptsnowV=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptgraulV=np.zeros(N[0:2],order='F',dtype=np.float32)
    SR=np.zeros(N[0:2],order='F',dtype=np.float32)
    
    dz=np.array(np.gradient(atm.hgt)[0].transpose([1,0,2]),order="F",dtype=np.float32)
    
    U=np.array(oldatm.u.transpose([1,0,2]),order='F',dtype=np.float32) # vapor pressure
    V=np.array(oldatm.v.transpose([1,0,2]),order='F',dtype=np.float32) # cloud water
    W=np.array(oldatm.w.transpose([1,0,2]),order='F',dtype=np.float32) # rain in air
    newU=np.array(atm.u.transpose([1,0,2]),order='F',dtype=np.float32)
    newV=np.array(atm.v.transpose([1,0,2]),order='F',dtype=np.float32)
    newW=np.array(atm.w.transpose([1,0,2]),order='F',dtype=np.float32)
    p=np.array(oldatm.p.transpose([1,0,2]), order='F',dtype=np.float32)
    p2=np.array(atm.p.transpose([1,0,2]), order='F',dtype=np.float32)
    
    t=0.0
    itimestep=0
    
    # calculate the maximum advection dt which still meets the Courant condition (here 0.5 for 2d 0.33 for 3d)
    # dt=0.49/np.max(np.abs(U[1:,:,:])/dx+np.abs(V[:,:,1:])/dx) # for 2d advection
    dt=0.3/max((np.max(np.abs(W[1:,:,1:])+np.abs(U[1:,:,:])+np.abs(V[:,:,1:]))/dx),
            (np.max(np.abs(newW[1:,:,1:])+np.abs(newU[1:,:,:])+np.abs(newV[:,:,1:]))/dx)) # for 3d advection?
    # physics could be speed up by converting to 1D advection (courant must be <1 so ~3x fewer time steps)
    #  but then physics should provide alternating direction advection calculations.
    # and since we try to keep dt<~60 for microphysics anyway, this may not matter too much
    # maybe look at different microphysics (e.g. cuda implementation of WSM6?) :)
    if dt>60:dt=60.0 # maximum dt = 120? 60? 30? seconds even if winds happen to be 0 so that the microphysics will be stable? 
    # ensure that dt goes into timestep evenly
    ntimes=np.ceil(timestep/dt)
    dt=timestep/ntimes
    n=U.shape
    
    dth=np.array(((atm.th-oldatm.th)/ntimes).transpose([1,0,2]),order="F",dtype=np.float32)
    dqv=np.array(((atm.qv-oldatm.qv)/ntimes).transpose([1,0,2]),order="F",dtype=np.float32)
    du=np.array((newU-U)*dt/dx/ntimes,order="F",dtype=np.float32)# *dt/timestep
    dv=np.array((newV-V)*dt/dx/ntimes,order="F",dtype=np.float32)# *dt/timestep
    dw=np.array((newW-W)*dt/dx/ntimes,order="F",dtype=np.float32)# *dt/timestep
    dp=np.array((p2-p)/ntimes,order="F",dtype=np.float32)
    # set up domain parameters and assume a two layer model 
    ids=1;ide=N[0];jds=1;jde=N[1];kds=1;kde=n[1];
    ims=ids;ime=ide;jms=jds;jme=jde;kms=kds;kme=kde;
    # for now set it up to run the entire domain (except the outer border that is used for the boundary conditions)
    its=ids+1;ite=ide-1;jts=jds+1;jte=jde-1;kts=kds;kte=kde;
    # its=50;ite=350;jts=200;jte=400;kts=kds;kte=kde;
    
    qv[qv<1E-15]=1E-15
    qc[qc<1E-15]=1E-15
    qr[qr<1E-15]=1E-15
    nr[nr<1E-15]=1E-15
    qs[qs<1E-15]=1E-15
    qi[qi<1E-15]=1E-15
    ni[ni<1E-15]=1E-15
    qg[qg<1E-15]=1E-15
    pii=np.array(1.0/((100000.0/p)**(R/cp)),order="F",dtype=np.float32)
    U*=dt/dx
    V*=dt/dx
    W*=dt/dx # note vertical coordinates are in pressure space (dz!=dx), 
             # but w was calculating to balance u,v assuming dz=dx
    # W*=-1
    # dw*=-1
    # print("W*=-1")
    print("dt="+str(dt))

    # SURFACE FLUXES
    sensible_heat=np.array(oldatm.sfc.sensible_heat,order="F",dtype="f")
    latent_heat=np.array(oldatm.sfc.latent_heat,order="F",dtype="f")
    pblh=np.array(oldatm.sfc.pblh,order="F",dtype="i")
    # SURFACE FLUXES
    physics=int(physics)
    t1=time.time()
    print("FORTRAN setup="+str(t1-t0))
    # print(U[100:104,13,99:104]/dt*dx)
    # print(V[100:104,13,99:104]/dt*dx)
    # print(W[100:104,13,99:104]/dt*dx)
    swim.swim_step.timestep(ntimes,U,V,W,
                sensible_heat, latent_heat,pblh,
                qv,qc,qr,qi,qs,qg,ni,nr, 
                th, pii, p, dz, dt, itimestep, 
                dth,dqv,du,dv,dw,dp,
                pptrain, pptrainV, pptsnow, pptsnowV, pptgraul, pptgraulV, SR, physics, 
                int(ids),int(ide), int(jds),int(jde), int(kds),int(kde), # domain dims
                int(ims),int(ime), int(jms),int(jme), int(kms),int(kme), # memory dims
                int(its),int(ite), int(jts),int(jte), int(kts),int(kte))
    t2=time.time()
    print("FORTRAN runtime="+str(t2-t1))
                
    atm.th=th.transpose([1,0,2])
    atm.qv=qv.transpose([1,0,2])
    atm.qc=qc.transpose([1,0,2])
    atm.qr=qr.transpose([1,0,2])
    atm.nr=nr.transpose([1,0,2])
    atm.qi=qi.transpose([1,0,2])
    atm.ni=ni.transpose([1,0,2])
    atm.qs=qs.transpose([1,0,2])
    atm.qg=qg.transpose([1,0,2])

    # print("FORTRAN takedown="+str(time.time()-t2))
    return pptrain+pptgraul+pptsnow
    
