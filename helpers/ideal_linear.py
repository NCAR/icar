#!/usr/bin/env python
from __future__ import print_function
import sys

import matplotlib.pyplot as plt #graphics library
import numpy as np              #numerical library (FFT etc)

from bunch import Bunch         # simple structure type
import Nio                  # NCAR python io library

# set up the physical domain of the experiment to be run
#    zs   = topography (Fzs = FFT(zs)) 
#    dx   = grid cell size (m)
#    U    = wind speed (m/s)
#    Hw   = water vapor scale height
#    Ndsq = Brunt Vaisalla freq. squared (s^-2)
# also sets up parameters for the linear model (e.g. tauf,c, cw)
#    tauf = average hydrometeor fall time (seconds) ~500s? for rain w/Hw=5000, 1500s? for snow w/Hw=3000
#    tauc = average hydrometeor formation time (seconds) ~500s for rain, longer for snow, faster for warmer
def setup_experiment(wind=2,experiment=1, verbose=False):
    U     = [5.0,10.0,15.0,25.0][wind]       #  wind speed
    # print("wind=",U,"  exp=",experiment)
    # Ndsq  = 3.6e-5;                          #  dry BV freq sq. #original
    # Ndsq  = 0.002**2    #4e-6                  #  dry BV freq sq.
    # Ndsq  = 0.00011     #1e-4                  #  dry BV freq sq.
    Ndsq  = 0.000011     #1e-4                  #  dry BV freq sq.


    # experiment                        D1      D2      D3
    # h         height of the hill    1800    1400    1040   [meters]
    # sigma     half-width              60      40       3.1 [grid cells]
    # z0        base of the hill      1700    2000    2200   [meters]
    # G         number of grids        420     250      52   [grid cells]

    Nx    = [420.,250.,52.][experiment]*2            #  length of domain  (grid cells)
    hm    = [1800.0,1400.0,1040.0][experiment]    #  mnt height (m)
    xm    = Nx/2.0                                #  mountain location in domain (grid cell)
    am    = [60.0,40.0,3.1][experiment]           #  mountain half-width (grid cells)

    dx    = 2000.0                                #  grid spacing (m)
    Lx    = Nx*dx                                 #  length of domain (m)
    x     = np.linspace(0,Lx,Nx)                  #  distance array (m)
    zo    = [1700.0,2000.0,2200.0][experiment]    # mountain base height (m)  NOT REALLY USED CORRECTLY YET
    zo    = 0.0
    p0    = 101325 * (1 - 2.25577e-5*zo)**5.25588   # compute base pressure
    # p0    = 101325.0
    T2m   = 268.0# 270.56 #needs to be selected for each experiment...?

    # hw    = 4000.0                                #  scale of water vapor (see formula in SB04,appendix)
    base_mr=[0.00325,0.0032,0.003][experiment]    # mixing ratio at the base of the mountain
    base_mr = 0.003255
    base_mr = 0.0025687
    # hw=hw-zo
    # hw=3000.0
    # zo=0.0
    # ----------------------------------------------------------------
    #  Make the mountain (theoretical)
    # 
    #   zs=hm*exp(-(x-xm).^2/am^2);   Gaussian
    zs=hm/(1.0+((x/dx-xm)/am)**2.)  # eqn from Trude
    zs=zs-zs[Nx/4] #sets the zero point to be 1/4 of the way in because we have doubled the size of the domain
    zs[zs<0]=0          #set anything below 0 to 0
    # zs+=zo
    # -----------------------------------------------------------------
    #  put zs in Fourier space
    Fzs  = np.fft.fftshift(np.fft.fft(zs))/Nx

    #  linear model paramters (see calculations in SB'04):
    #  -------------------------------------------------------------------------------
    t0 = 273.16
    # p0 = 100000
    # p0 =  82000.0 (now calculated above from z0)
    L   = 2.5e6
    ratio = 18.015/28.964
    R  = 287.0
    Rv = 461.0
    cp = 1004.0
    g  = 9.81
    es = 611.21*np.exp(17.502*(T2m-t0)/(T2m-32.19));
    qs0 = ratio *es/(p0-es);
    cap_gamma = -(g * (1.+(L*qs0)/(R*T2m)) / (cp + (L**2*qs0*ratio) / (R*T2m**2)));
    env_gamma = Ndsq*T2m/g + cap_gamma #Ndsq calculated from potential temperature profile cap_gamma converts to real temp?
    hw = np.abs((Rv*T2m**2)/(L*env_gamma))
    # if env_gamma pulled from model, enforce reasonable values with Ndsq=min(Ndsq,0.012)
    #     cw below calculated from Trude's profile, Ndsq=0.00011, T2m=271K, p0=820mb,  dth/dz=0.004K/m dt/dz=-0.0054
    # could calculate Ndsq as = (env_gamma-cap_gamma)*g/T2m
    Ndsq=(-0.0054-cap_gamma)*g/T2m
    #  -------------------------------------------------------------------------------
    # cw  = 1.9                # sensitivity (commonly set to 1.1, see paper SB04) = cap_gamma / env_gamma
    cw  = cap_gamma/env_gamma
    # using base_mr from profile, but maybe it should be qs0 above?
    cwqv= cw*base_mr            # =sensitivity times q_vs (set to: 0.01 kg/kg, see paper SB04)
    # print(cwqv,cw,base_mr,qs0,p0,hw)
    z0  = 0                 # at the height (AGL?) where you want the precip

    vterm = 2.0              # vertical terminal velocity for e.g. snow = 2m/s rain=10m/s
    tauf= hw/vterm           # BS'11: =zg/Vt ( =hw/Vt for moisture level: around 500s (->750s) is normally good)
    tauc= 2000.0             # cloud->hydrometero conversion time.  probably around 500s for rain, 
                             # shorter for warmer condition, longer for snow?
    if verbose:
        print("   Ndsq=",Ndsq)
        print("   Environmental lapse rate=0.004K/m")
        print("   \"Dry\" lapse rate=0.0017K/m")
        print("   Base MR=",base_mr)
        print("   Scale height=",hw)
        print("   tauc=",tauc)
        print("   tauf=",tauf)
        print("   cwqv=",cwqv)
    # ---------------------------------------------------------------------------------
    params=Bunch(cw=cw,cwqv=cwqv,z0=z0,tauf=tauf,tauc=tauc,hw=hw,Ndsq=Ndsq)
    return (x,zs,Fzs,U,dx,params)

def get_params(T2m,U,Ndsq,zs,env_gamma,verbose=False):
    """docstring for get_params"""
    Nx    = len(zs)            #  length of domain  (grid cells)
    # hm    = [1800.0,1400.0,1040.0][experiment]    #  mnt height (m)
    # xm    = Nx/2.0                                #  mountain location in domain (grid cell)
    # am    = [60.0,40.0,3.1][experiment]           #  mountain half-width (grid cells)

    # dx    = 2000.0                                #  grid spacing (m)
    # Lx    = Nx*dx                                 #  length of domain (m)
    # x     = np.linspace(0,Lx,Nx)                  #  distance array (m)
    # zo    = [1700.0,2000.0,2200.0][experiment]    # mountain base height (m)  NOT REALLY USED CORRECTLY YET
    zo    = 0.0
    p0    = 101325 * (1 - 2.25577e-5*zo)**5.25588   # compute base pressure
    # p0    = 101325.0
    # T2m   = 268.0# 270.56 #needs to be selected for each experiment...?

    # hw    = 4000.0                                #  scale of water vapor (see formula in SB04,appendix)
    # base_mr=[0.00325,0.0032,0.003][experiment]    # mixing ratio at the base of the mountain
    # base_mr = 0.003255
    # base_mr = 0.0025687
    # -----------------------------------------------------------------
    #  put zs in Fourier space
    Fzs  = np.fft.fftshift(np.fft.fft(zs))/Nx

    #  linear model paramters (see calculations in SB'04):
    #  -------------------------------------------------------------------------------
    t0 = 273.16
    # p0 = 100000
    # p0 =  82000.0 (now calculated above from z0)
    L   = 2.5e6
    ratio = 18.015/28.964
    R  = 287.0
    Rv = 461.0
    cp = 1004.0
    g  = 9.81
    es = 611.21*np.exp(17.502*(T2m-t0)/(T2m-32.19))
    qs0 = ratio *es/(p0-es)
    base_mr=qs0
    
    cap_gamma = -(g * (1.+(L*qs0)/(R*T2m)) / (cp + (L**2*qs0*ratio) / (R*T2m**2)))
    # env_gamma = Ndsq*T2m/g + cap_gamma #Ndsq calculated from potential temperature profile cap_gamma converts to real temp?
    hw = np.abs((Rv*T2m**2)/(L*env_gamma))
    # if env_gamma pulled from model, enforce reasonable values with Ndsq=min(Ndsq,0.012)
    #     cw below calculated from Trude's profile, Ndsq=0.00011, T2m=271K, p0=820mb,  dth/dz=0.004K/m dt/dz=-0.0054
    # could calculate Ndsq as = (env_gamma-cap_gamma)*g/T2m
    # Ndsq=(-0.0054-cap_gamma)*g/T2m
    #  -------------------------------------------------------------------------------
    # cw  = 1.9                # sensitivity (commonly set to 1.1, see paper SB04) = cap_gamma / env_gamma
    cw  = cap_gamma/env_gamma
    if verbose:
        print(cap_gamma, env_gamma, cw, qs0)
    # using base_mr from profile, but maybe it should be qs0 above?
    cwqv= cw*base_mr            # =sensitivity times q_vs (set to: 0.01 kg/kg, see paper SB04)
    # print(cwqv,cw,base_mr,qs0,p0,hw)
    z0  = 0                 # at the height (AGL?) where you want the precip

    vterm = 2.0              # vertical terminal velocity for e.g. snow = 2m/s rain=10m/s
    tauf= hw/vterm           # BS'11: =zg/Vt ( =hw/Vt for moisture level: around 500s (->750s) is normally good)
    tauc= 2000.0             # cloud->hydrometero conversion time.  probably around 500s for rain, 
                             # shorter for warmer condition, longer for snow?
    if verbose:
        print("   Ndsq=",Ndsq)
        print("   Environmental lapse rate=0.004K/m")
        print("   \"Dry\" lapse rate=0.0017K/m")
        print("   Base MR=",base_mr)
        print("   Scale height=",hw)
        print("   tauc=",tauc)
        print("   tauf=",tauf)
        print("   cwqv=",cwqv)
    # ---------------------------------------------------------------------------------
    params=Bunch(cw=cw,cwqv=cwqv,z0=z0,tauf=tauf,tauc=tauc,hw=hw,Ndsq=Ndsq)
    return (Fzs,params)
    

# 
#     This is a 2-D precipitation solution for the linear model. 
# 
def solve(Fzs,U,dx,params, zlevels=None):
    cw  = params.cw   #1.1    #efficiency factor
    cwqv= params.cwqv
    z0  = params.z0   #0      #meters
    tauf= params.tauf #1500   #seconds
    tauc= params.tauc #1000   #seconds
    hw  = params.hw   #3000   #m water vapor scale height, maybe this shouldn't be in params...
    Ndsq= params.Ndsq #0.0011 # 1/s^2 Brunt-Vaisalla frequency squared
    
    Nx=Fzs.shape[0]
    m=np.ones(Nx)
    FSterm=np.zeros(Nx).astype('complex')
    FSraw=np.zeros(Nx).astype('complex')
    FPterm=np.zeros(Nx).astype('complex')
    i=np.sqrt(np.complex(-1))

    k    = np.linspace(-np.pi/dx,np.pi/dx,Nx)
    sig  = U*k
    denom = sig**2 # -cori**2 # don't bother with coriollis force in 2D case
    # discrim=Ndsq/U**2-k**2
    # for ik in np.arange(Nx):
    #  nonhydrostatic version, needed for smaller-scale mountains(?)
    #   if (discrim(ik) > 0)
    #     m(ik) = sign(sig(ik))*sqrt(discrim(ik));
    #   else
    #     m(ik) = sqrt(discrim(ik));                       Check before use
    #   end
    # nonhydrostatic version needed for smaller-scale mountains(?). 
    msq = (Ndsq/denom * (k**2)).astype('complex')          #  vertical wave number, hydrostatic
    mimag=np.zeros(Nx).astype('complex')
    mimag.imag=(np.sqrt(-msq)).real
    m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)
    # m = np.sign(sig)*np.abs(np.sqrt(Ndsq)/U)          #  vertical wave number, hydrostatic (eqn 14 in SB'04 (?))

    # FSraw  = i*cwqv*Fzs*sig*np.exp(-z0/hw)                            #  the raw upslope (~Smith'79)
    # FSterm = (i*cwqv*Fzs*sig/(1-i*m*hw))                                #  simple, hydrostatic solution (eqn 16 in SB'04)
    FSterm = (i*cwqv*Fzs*np.exp(i*m*z0)*sig*np.exp(-z0*(1-i*m*hw)/hw)/  #  simple, hydrostatic solution (eqn 16 in SB'04)
                 (1-i*m*hw))                                            # see BS'11 for z0 component
    FPterm = FSterm/((1+i*sig*tauc)*(1+i*sig*tauf))                     #  this also include time delays and dynamics (eqn 5 in SB'04)

    # calculate "3D" wind field for Nz levels ranging from 0m to Hw*5 m AGL
    if zlevels!=None:
        Nz=len(zlevels)
    else:
        Nz=20
        Zrange=hw*5.0
    
    Fw=np.zeros((Nz,Nx)).astype("complex")
    Fu=np.zeros((Nz,Nx)).astype("complex")
    z3d=np.zeros((Nz,Nx))
    k2=k**2
    k2[k2==0]=1E-15
    sig[sig==0]=1E-15
    for curz in range(Nz):
        if zlevels!=None:
            zlevel=zlevels[curz]
        else:
            zlevel=curz*Zrange/Nz
        z3d[curz,:]=zlevel
        neta=Fzs*np.exp(i*m*zlevel)
        # Fw[curz,:] = i*sig*Fzs*np.exp(i*m*zlevel)
        Fw[curz,:]=i*sig*neta
        # Fu[curz,:]=(-m*(sig*k-i*l*f)*i*neta)/kl
        Fu[curz,:]=(-m*(sig*k)*i*neta)/k2

    #  bring all back from Fourier space:
    # S      = Nx*np.fft.ifft(np.fft.ifftshift(FSterm))
    # Sraw   = Nx*np.fft.ifft(np.fft.ifftshift(FSraw))
    Pt     = np.real(Nx*np.fft.ifft(np.fft.ifftshift(FPterm)))
    w=np.zeros((Nz,Nx))
    u_hat=np.zeros((Nz,Nx))
    for curz in range(Nz):
        w[curz,:]=Nx*np.real(np.fft.ifft(np.fft.ifftshift(Fw[curz,:])))
        u_hat[curz,:]=Nx*np.real(np.fft.ifft(np.fft.ifftshift(Fu[curz,:])))

    #  trucate negative values as unphysical:
    Pt[np.where(Pt<=0)]=0
    return (Pt,w,U+u_hat,z3d)
    
def get_winds(x,z,xi,zi,U,W):
    Nz,Nx=z.shape
    if len(x.shape)==1:
        x=x[np.newaxis,:].repeat(Nz,axis=0)
    i,j=np.unravel_index(np.argmin((x-xi)**2+(z-zi)**2),(Nz,Nx))
    try:
        if len(U)>1:
            curU=U[i,j]
    except TypeError:
        curU=U
    curW=W[i,j]
    return(curU,curW)

def bilin_weights(yi,y,xi,x):
    x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1.0-x0 # equivalent to: np.abs((xi-x[1])/(x[1]-x[0]))
    x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1.0-x2 # equivalent to: np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    if (y6-y5)==0:
        f1=0.5
    else:
        f1=(yi-y5)/(y6-y5)
    f2=1.0-f1 # equivalent to: (y6-yi)/(y6-y5)
    return np.array([x1*f2,x0*f2,x3*f1,x2*f1])

def get_bilin_winds(x,z,xi,zi,U,W):
    Nz,Nx=z.shape
    if len(x.shape)==1:
        x=x[np.newaxis,:].repeat(Nz,axis=0)
    i,j=np.unravel_index(np.argmin((x-xi)**2+(z-zi)**2),(Nz,Nx))
    if x[i,j]>xi:
        if z[i,j]>zi:
            xpos=np.array([j,j-1,j,j-1])
            zpos=np.array([i,i,i-1,i-1])
        else:
            xpos=np.array([j,j-1,j,j-1])
            zpos=np.array([i,i,i+1,i+1])
    else:
        if z[i,j]>zi:
            xpos=np.array([j,j+1,j,j+1])
            zpos=np.array([i,i,i-1,i-1])
        else:
            xpos=np.array([j,j+1,j,j+1])
            zpos=np.array([i,i,i+1,i+1])
    xpos[xpos<0]=0
    zpos[zpos<0]=0
    xpos[xpos>=Nx]=Nx-1
    zpos[zpos>=Nz]=Nz-1
    weights=bilin_weights(zi,z[zpos,xpos],xi,x[zpos,xpos])

    try:
        if len(U)>1:
            curU=(U[zpos,xpos]*weights).sum()
    except TypeError:
        curU=U
    curW=(W[zpos,xpos]*weights).sum()
    # curW=W[z0,x0]*w0+W[z1,x1]*w1+W[z2,x2]*w2+W[z2,x2]*w2
    return(curU,curW)

def plot_stream(x,z,start_level,allU,allw):
    dx=x[1]-x[0]
    dz=z[1,0]-z[0,0]
    Nz=z.shape[0]
    x=x[np.newaxis,:].repeat(Nz,axis=0)
    # set dt to at most move half a grid cell between time steps. 
    # dt=(dx/allU+np.abs(dz/allw)).min() /20.0
    dt=10.0
    i,j=start_level,0
    curx=0.0
    curz=float(z[i,j])
    allx,allz=[],[]
    # generate the stream by tracking a particle
    # through the domain in a lagrangian model
    while (curx<x[-1,-1]):
        # print(curx)
        allx.append(curx)
        allz.append(curz)
        # u,w=get_winds(x,z,curx,curz,allU,allw)
        u,w=get_bilin_winds(x,z,curx,curz,allU,allw)
        curx+=u*dt
        curz+=w*dt
        # print(curx,u,dt)
        # print(curz,w,dt)
    # print(np.min(allx),np.min(allz),np.max(allx),np.max(allz))
    plt.plot(np.array(allx)/1000.0,allz,color='grey')

# note this is painfully slow given the huge number of lines/"arrows" generated
# could be trimmed down by decreasing the number of x cells looked at. could also use quiver 
def plot_arrows(x,z,U,w):
    Nz,Nx=z.shape
    maxw,maxu=w.max(),U.max()
    maxz,maxx=z.max(),x.max()
    dx=maxx/50.0
    dz=maxz/50.0
    allx=np.arange(0,maxx,dx)
    allz=np.arange(0,maxz,dz)
    x3d=x[np.newaxis,:].repeat(Nz,axis=0)
    for xi in allx:
        for zi in allz:
            ui,wi=get_bilin_winds(x3d,z,xi,zi,U,w)
            curpos=np.argmin(np.abs(x-xi))
            if zi<=z[0,curpos]:
                ui=0
                wi=0
            plt.plot([xi/1000.0,(xi+ui*dx/maxu*0.9)/1000.0],[zi,zi+wi*dz/maxw*0.9],color="black")
    
def plot_output(x,zs,Pt,z3d,w,U,title,fname):
    #  P-term:
    plt.clf()
    f=plt.gcf()
    # plot precipitation
    plt.plot(x/1000,3600*Pt.real*500,'b',label='Precip.',linewidth=3.0)

    # make a polygon to fill in for the topography
    xpoly=list()
    xpoly.append(0)
    xpoly.extend(list(x/1000))
    xpoly.append(10000)
    zpoly=list()
    zpoly.append(-100)
    zpoly.extend(list(zs))
    zpoly.append(-100)
    # plot the topography
    plt.fill(xpoly,zpoly,'#99ff99',label='Elevation')
    # set axis labels and limits
    plt.xlabel('Distance (km)')
    plt.ylabel('Elevation (m)  --   Precip. [mm/hr]*500')
    plt.xlim(x[0]/1000.0,x[-1]/1000.0)
    # plt.ylim(0,z3d.max())
    plt.ylim(0,10000)
    # set graph title
    # plt.title("Wind Speed="+str(U)+"m/s")
    plt.title(title)
    # add a legend
    plt.legend()
    # plot 3D wind field
    Nz=z3d.shape[0]
    z3d2=z3d+zs[np.newaxis,:]
    plot_arrows(x,z3d2,U,w)
    # alternatively, trace streamlines, but this takes a long time
    #  and streamlines sometimes diverge
    # for i in range(Nz):
    #     print(float(i)/Nz*100.0,end="% ")
    #     plot_stream(x,z3d2,i,U,w)
    #     sys.stdout.flush()
    # print("100%")
    
    # save figure as an image
    print(fname)
    f.savefig(fname)
    
def write_output(fname,x,z,p,w,u=None):
    sz=z.shape
    xpad=sz[1]/4
    Nz,Nx=z.shape
    Nx-=xpad*2
    ncf=Nio.open_file(fname,mode="w",format="nc")
    ncf.create_dimension("x",Nx)
    ncf.create_dimension("z",Nz)

    varname="x"
    x=x[:sz[1]/2]
    ncf.create_variable(varname,'f',('x',))
    ncf.variables[varname][:]=x.astype('f')
    ncf.variables[varname].units="m"
    ncf.variables[varname].description="x coordinate"

    varname="z"
    z=z[:,xpad:-xpad]
    ncf.create_variable(varname,'f',('z','x'))
    ncf.variables[varname][:]=z.astype('f')
    ncf.variables[varname].units="m"
    ncf.variables[varname].description="z coordinate"

    varname="precipitation"
    p=p[xpad:-xpad]
    ncf.create_variable(varname,'f',('x',))
    ncf.variables[varname][:]=p.astype('f')
    ncf.variables[varname].units="kg/m^2/s"
    ncf.variables[varname].description="precipitation rate"
    
    varname="w"
    w=w[:,xpad:-xpad]
    ncf.create_variable(varname,'f',('z','x'))
    ncf.variables[varname][:]=w.astype('f')
    ncf.variables[varname].units="m/s"
    ncf.variables[varname].description="vertical velocity"

    if u!=None:
        varname="u"
        u=u[:,xpad:-xpad]
        ncf.create_variable(varname,'f',('z','x'))
        ncf.variables[varname][:]=u.astype('f')
        ncf.variables[varname].units="m/s"
        ncf.variables[varname].description="horizontal velocity"

    ncf.close()

def write_2d_output(fname,lat,lon,p):
    
    Ny,Nx=p.shape
    ncf=Nio.open_file(fname,mode="w",format="nc")
    ncf.create_dimension("lat",Ny)
    ncf.create_dimension("lon",Nx)

    varname="XLAT"
    ncf.create_variable(varname,'f',('lat','lon'))
    ncf.variables[varname][:]=lat.astype('f')
    ncf.variables[varname].units="degrees"
    ncf.variables[varname].description="Latitude"

    varname="XLONG"
    ncf.create_variable(varname,'f',('lat','lon'))
    ncf.variables[varname][:]=lon.astype('f')
    ncf.variables[varname].units="degrees"
    ncf.variables[varname].description="Longitude"

    varname="precipitation"
    ncf.create_variable(varname,'f',('lat','lon'))
    ncf.variables[varname][:]=p.astype('f')
    ncf.variables[varname].units="kg/m^2/s"
    ncf.variables[varname].description="precipitation rate"
    
    ncf.close()

  
def main():
    wind=[3,2,1,0]
    experiments=[2,1,0]
    wind=[1]
    experiments=[1]
    # curwind=2
    # curexp=2
    for curwind in wind:
        for curexp in experiments:
            # print(curwind,curexp)
            (x,zs,Fzs,U,dx,params)=setup_experiment(curwind,curexp)
            # for nd in np.arange(1e-6,0.01,5e-5):
            # nd=0.002**2
            # params.Ndsq=nd
            (Pt,w,U3d,z3d)=solve(Fzs,U,dx,params)
            U3d[U3d<=1E-5]=1E-5
            fname="U_"+str(U)+"_exp_"+str(curexp)
            title="wind = "+str(U)+"m/s"+"  Experiment:"+str(curexp)
            print(title)
            plot_output(x,zs,Pt,z3d,w,U3d*0+U,title,fname+".png")
            write_output(fname,x,z3d+zs[np.newaxis,:],Pt,w,u=U3d)

if __name__ == '__main__':
    main()

# Some Linear Model papers :
#         SB'04  :: Smith and Barstad (2004); JAS
#         BGS'07 :: Barstad et al. 2007; J Hydrol.
#         BS'11  :: Barstad and Schuller (2011) Nov, JAS
