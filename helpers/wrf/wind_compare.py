#!/usr/bin/env python
import sys
import glob
import copy
import numpy as np
import matplotlib.pyplot as plt

# from lt_winds import linear_winds
import ideal_linear

import mygis
from bunch import Bunch

dz_levels=np.array([
		536.017,531.413,527.025,522.853,518.894,515.153,511.639,508.358,505.32,502.537
		,500.027,497.807,495.897,494.33,493.125,492.323,491.962,492.091,492.762,494.038
		,495.997,498.733,502.353,506.99,512.811,520.02,528.88,539.72,552.979,569.233
		,589.277,614.229,645.693,686.135,739.443,812.309,917.201,1080.76,1372.97,2066.42])
dx=2000.0

def load_wrf(filename, preciponly=False):
    """docstring for load_wrf"""
    time=15
    precip=mygis.read_nc(filename,"RAINNC").data
    precip=precip[time+25,1,:] - precip[time-5,1,:]
    precip/=15.0 # 30 outputsteps = 15 hours
    if preciponly:
        return Bunch(precip=precip,hgt=None)

    w=mygis.read_nc(filename,"W").data[time,:,1,:]
    u=mygis.read_nc(filename,"U").data[time,:,1,:]
    z=mygis.read_nc(filename,"PH").data[time,:,1,:]
    z+=mygis.read_nc(filename,"PHB").data[time,:,1,:]
    z/=9.81

    hgt=mygis.read_nc(filename,"HGT").data[time,1,:]

    return Bunch(w=w,z=z,hgt=hgt,u=u,precip=precip)

def load_icar(filename,h,preciponly=False):
    """docstring for load_wrf"""
    # run_num=filename.split("out")[-1]
    # last_file=filename.replace(run_num,("{0:0"+str(len(run_num))+"}").format(int(run_num)-1))
    if type(filename)==list:
        output=[]
        for f in filename:
            output.append(load_icar(f,h,preciponly=preciponly))
        return output

    precip=mygis.read_nc(filename,"rain").data
    precip=precip[11,1]-precip[10,1]
    if preciponly:
        return Bunch(precip=precip)

    u=mygis.read_nc(filename,"u").data[0,1,...]
    # dudx=-np.diff(u,axis=1)
    dudx=u[:,:-1]-u[:,1:]

    # w=mygis.read_nc(filename,"w").data[1,...]
    z=dudx*0
    z[0,:]=h+dz_levels[0]/2
    nz=z.shape[0]
    for i in range(1,nz):
        z[i,:]=z[i-1,:]+(dz_levels[i-1]+dz_levels[i])/2

    # print(z[:,0])
    dz=np.diff(h)
    # w1=u[:,1:-1]*np.sin(np.arctan(dz[np.newaxis,:]/dx))
    # w1=u[:,1:-1]*dz[np.newaxis,:]/np.sqrt(dz[np.newaxis,:]**2 + dx**2)
    w1=u[:,1:-1]*dz[np.newaxis,:]/dx
    w1=(w1[:,1:]+w1[:,:-1])/2.0 # * dz_levels[:,np.newaxis]/dx # *10.0

    w2=dudx*dz_levels[:,np.newaxis]/dx
    for i in range(1,nz):
        w2[i,:]+=w2[i-1,:]
    w2=w2[:,1:-1]

    w=z*0
    w[:,1:-1]=w1+w2
    # w[:,1:-1]=w2
    return Bunch(w=w,z=z,hgt=h,u=u,precip=precip)

def load_linear(zs,T2m=260,u=10.0,v=0,levels=np.array([250,750]),Ndsq=6e-5,dthdz=3.0):
    """docstring for load_linear"""

    mean_wind=np.round(np.mean(u))
    U = mean_wind # dont ask... need to figure this out though
    # print(U)
    # print(T2m)
    V = 0
    z=0.0
    dx,dy=2000.0,2000.0
    env_gamma = -dthdz/1000.0

    if (len(zs)%2) == 1:
        zs=zs[:-1]

    (Fzs,params) = ideal_linear.get_params(T2m,U,Ndsq,zs,env_gamma)
    # params.tauf*=2
    # print(params)
    (Pt,w,U3d,z3d)=ideal_linear.solve(Fzs,U,dx,params,zlevels=levels)

    return Bunch(w=w,z=levels,z3d=z3d,hgt=zs,u=U3d,precip=Pt*3600)


    # old code for use with lt_winds module, which doesn't seem to work
    # hgt=hgtin.repeat(2,axis=0)
    # Ny,Nx=hgt.shape
    # Fzs   = np.fft.fftshift(np.fft.fft2(hgt))/(Nx*Ny)
    # dz=np.diff(hgt,axis=1)
    # u_out=np.zeros(levels.shape)
    # w_out=np.zeros(levels.shape)
    # # p,u_hat,v_hat,w_hat = linear_winds(Fzs, U,V,z,dx=dx,dy=dy,Ndsq=Ndsq,outputP=True)
    # p*=3600.0
    # zlevels=levels[:,0]
    # i=0
    # for z in zlevels:
    #     # u_hat,v_hat,w_hat = linear_winds(Fzs, U,V,z,dx=dx,dy=dy,Ndsq=Ndsq)
    #     u_hat2=(u_hat[:,1:]+u_hat[:,:-1])/2 # put u_hat on the points between gridcells
    #     # w1=u[:,1:-1]*np.sin(np.arctan(dz[np.newaxis,:]/dx))
    #     # w1=u[:,1:-1]*dz[np.newaxis,:]/np.sqrt(dz[np.newaxis,:]**2 + dx**2)
    #     w1=(mean_wind+u_hat2)*dz/dx
    #     w1=(w1[:,1:]+w1[:,:-1])/2.0 # put w1 back on to grid centers
    #     w = np.zeros(hgt.shape) + w_hat
    #     w_out[i,...]=w[0,:]
    #     u_out[i,...]=mean_wind+u_hat[0,:]
    #     i+=1
    # return Bunch(w=w_out,z=levels,hgt=hgt,u=u_out,precip=p)



def interp2point(w,z,point,verbose=False):
    """docstring for interp2point"""
    if point>z[0]:
        for i in range(1,len(z)):
            if point<=z[i]:
                weight=(z[i]-point)/(z[i]-z[i-1])
                if verbose:print(weight,point,z[i-1])
                return w[i]*(1-weight)+w[i-1]*weight
    else:
        return w[0]

def vinterp(d1,d2):
    """docstring for vinterp"""
    newdata=copy.deepcopy(d2)
    nz,nx=newdata.z.shape
    newdata.precip=d1.precip
    newdata.u=d1.u
    for x in range(nx):
        for z in range(nz):
            newdata.w[z,x]=interp2point(w=d1.w[:,x],z=d1.z[:,x],point=d2.z[z,x])#,verbose=(x==nx/2))
            # newdata.u[z,x]=interp2point(w=d1.u[:,x],z=d1.z[:,x],point=d2.z[z,x])#,verbose=(x==nx/2))

    return newdata

def main(wrf_file,icar_file,output_file,makeplot=True,plot_legend=False,
         precip_max=2.0,getPrecip=True,Ndsq=None,t=None,add_file=None,master=False,
         text_font_size=8, legend_font_size=8, label_font_size=8,
         draw_xlabels=True, draw_ylabels=True):
    """docstring for main"""
    if Ndsq==None:
        Ndsq=icar_file.split("_ns")[1][:3]
        Ndsq=float(Ndsq.replace("e","e-"))
    if t==None:
        case=wrf_file.split("/")[0]
        t=float(case.split("_")[3])
    # print(Ndsq,t)
    dx=2.0 # km
    wrfdata=load_wrf(wrf_file,preciponly=False)
    icardata=load_icar(icar_file,wrfdata.hgt,preciponly=False)
    if add_file:
        icar2=load_icar(add_file,wrfdata.hgt,preciponly=False)
    linear_data=load_linear(wrfdata.hgt,T2m=t,u=icardata.u,v=0,levels=icardata.z[:,0],Ndsq=Ndsq,dthdz=3.0)

    if not getPrecip:
        iwrf=vinterp(wrfdata,icardata)
    else:
        iwrf=wrfdata


    # lim=10
    if makeplot:
        lim=2
        if not getPrecip:
            fig=plt.figure(figsize=(13,7))
            plt.subplot(221)
            # plt.imshow(icardata.u-icardata.u.mean(),cmap=plt.cm.seismic)
            plt.imshow(icardata.w,cmap=plt.cm.seismic)
            # plt.colorbar()
            plt.title("ICAR")
            plt.clim(-lim,lim)
            plt.xlim(150,250)

            plt.subplot(222)
            # plt.imshow(wrfdata.u-wrfdata.u.mean(),cmap=plt.cm.seismic)
            plt.imshow(iwrf.w,cmap=plt.cm.seismic)
            # plt.colorbar()
            plt.title("WRF")
            plt.clim(-lim,lim)
            plt.xlim(150,250)

            plt.subplot(223)
            # plt.imshow(wrfdata.u-wrfdata.u.mean(),cmap=plt.cm.seismic)
            # plt.imshow(iwrf.w,cmap=plt.cm.seismic)
            plt.imshow(linear_data.w,cmap=plt.cm.seismic)
            # plt.colorbar()
            plt.title("Linear Theory")
            plt.clim(-lim,lim)
            plt.xlim(150,250)

            plt.subplot(224)
            plt.plot(iwrf.w[0:5,150:250].mean(axis=0),label="WRF",linewidth=2,color="black")
            plt.plot(linear_data.w[0:5,150:250].mean(axis=0),label="Linear",linewidth=2,color="red")
            plt.plot(icardata.w[0:5,150:250].mean(axis=0),label="ICAR",linewidth=2,color="blue")
            plt.legend(loc=3,fontsize=10)
            plt.plot(iwrf.w[0,150:250],label="WRF0",color="black")
            plt.plot(linear_data.w[0,150:250],label="Linear0",color="red")
            plt.plot(icardata.w[0,150:250],label="ICAR0",color="blue")
            # plt.legend(loc=3,ncol=2,fontsize=10)
            plt.plot([0,100],[0,0],color="black",linestyle="--")

        if getPrecip:
            if not master:
                fig=plt.figure(figsize=(6,4.5))
            offset=0
        else:
            offset=-1.5
        precip_width=2
        # precip_max=3.5
        # precip_max=2.0
        x=np.arange(100)*dx

        plt.plot(x, linear_data.precip[150:250]+offset,color="red",label="Linear",linewidth=precip_width)
        plt.plot(x, wrfdata.precip[150:250]+offset,color="black",label="WRF",linewidth=precip_width)
        plt.plot(x, icardata.precip[150:250]+offset,color="green",label="ICAR$_t$",linewidth=precip_width)
        plt.plot([50*dx,50*dx],[offset,precip_max+offset],color="black",linestyle=":",linewidth=precip_width)
        print("case: "+" ".join(wrf_file.split("/")[0].split("_")[1:]))
        print("ICAR-t "+str(icardata.precip[150:250].mean()))
        print("WRF "+str(wrfdata.precip[150:250].mean()))
        print("linear "+str(linear_data.precip[150:250].mean()))
        if add_file:
            if type(add_file)==list:
                for i in range(len(add_file)):
                    # label=add_file[i].split("/")[-2].split("_")[0]
                    label=["ICAR$_l$","ICAR$_s$"][i]
                    color=["blue","orange"][i]
                    plt.plot(x, icar2[i].precip[150:250]+offset,color=color,label=label,linewidth=precip_width)
                    print(label+" "+str(icar2[i].precip[150:250].mean()))
            else:
                label=add_file.split("_")[0]
                label="ICAR$_l$"
                plt.plot(x, icar2.precip[150:250]+offset,color="blue",label=label,linewidth=precip_width)
                print("ICAR-l "+str(icar2.precip[150:250].mean()))

        if getPrecip:
            if master:
                plot=master
            else:
                plot = fig.add_subplot(111)

            case=wrf_file.split("/")[0].split("_")[1:]
            plt.text(10, (precip_max+offset)*0.75, " U   ={0[0]}m/s\n RH={0[1]}\n T   ={0[2]}K".format(case),fontsize=text_font_size)
            # if (case[1]=="0.75") and (case[2]=="260") and (case[0]=="5"):
            if plot_legend:
                plt.legend(fontsize=legend_font_size)

            # plot.tick_params(axis='both', which='major', labelsize=0)
            if draw_xlabels:
                # if case[0]=="20":
                plot.tick_params(axis='x', which='major', labelsize=label_font_size)
                plt.xlabel("Distance (km)", fontsize=label_font_size)
            else:
                plot.set_xticklabels("")
            if draw_ylabels:
                # if case[1]=="0.75" and case[2]=="260":
                plot.tick_params(axis='y', which='major', labelsize=label_font_size)
                plt.ylabel("Precipitation Rate (mm/hr)", fontsize=label_font_size)
            else:
                plot.set_yticklabels("")

        plt.ylim(offset,precip_max+offset)

        # plt.subplot(313)
        # plt.imshow(wrfdata.w)
        # plt.colorbar()
        # plt.title("WRF")
        # plt.xlim(150,250)

        case=wrf_file.split("/")[0]
        # print("OUTPUT = wfiles/"+case+"_"+output_file+"_w.png")
        if not master:
            plt.savefig("wfiles/"+case+"_"+output_file+"_w.png",dpi=300)
            plt.close()

    if add_file:
        return iwrf,icardata, icar2,linear_data
    else:
        return iwrf,icardata,linear_data


bad_cases=[
    "wind_15_0.99_280_3",
    "wind_10_0.99_280_3",
    "wind_5_0.99_280_3",
    "wind_10_0.9_280_3",
    "wind_5_0.9_280_3"]
# bad_cases=[]
def mean_precip(Nsquared="1e4"):
    """docstring for mean_precip"""
    cases=glob.glob("wind_*_*_*_3")
    output=[]
    full_data=[]
    for case in cases:
        if not (case in bad_cases):
            wrf_file=glob.glob(case+"/wrfout_*")[0]
            icar_file=glob.glob(case+"/thompson_output_2/icar_out2000_01_01*")[0]
            wrf,icar,linear=main(wrf_file,icar_file,Nsquared,makeplot=False,getPrecip=True,Ndsq=float(Nsquared.replace("e","e-")))

            output.append([wrf.precip[150:250].mean(),icar.precip[150:250].mean(),linear.precip[150:250].mean()])
            full_data.append([wrf.precip[150:250],icar.precip[150:250]])
        else:
            print("BAD:" + case)

    return np.array(output),full_data

if __name__ == '__main__':
    if len(sys.argv)==4:
        main(sys.argv[1],sys.argv[2],sys.argv[3])
    if len(sys.argv)==5:
        main(sys.argv[1],sys.argv[2],sys.argv[3],makeplot=True, getPrecip=False, Ndsq=6e-5,add_file=sys.argv[4])
        # main(sys.argv[1],sys.argv[2],sys.argv[3],makeplot=True, getPrecip=False, Ndsq=6e-5,add_file=sys.argv[4])
    #main(wrf_file,icar_file,output_file,makeplot=True,getPrecip=True,Ndsq=None,t=None,add_file=None)
