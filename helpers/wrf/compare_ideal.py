#!/usr/bin/env python
import sys
import glob
from bunch import Bunch
import numpy as np
import mygis
import matplotlib.pyplot as plt

def load_wrf(filename):
    """docstring for load_wrf"""
    data=[]
    datelist=[]
    qvdata=mygis.read_nc(filename,"QVAPOR").data
    qcdata=mygis.read_nc(filename,"QCLOUD").data \
            + mygis.read_nc(filename,"QICE").data \
            + mygis.read_nc(filename,"QSNOW").data \
            + mygis.read_nc(filename,"QRAIN").data
    tdata=mygis.read_nc(filename,"T").data + 300
    pdata=mygis.read_nc(filename,"P").data \
            + mygis.read_nc(filename,"PB").data
    rdata=mygis.read_nc(filename,"RAINNC").data
    udata=mygis.read_nc(filename,"U").data
    dates=mygis.read_nc(filename,"Times").data
    
    r0=np.zeros(rdata[0,1,:].shape)
    for i in range(0,qvdata.shape[0],2):
        outputdata=Bunch(qv=qvdata[i,:,1,:],
                         qc=qcdata[i,:,1,:],
                         t = tdata[i,:,1,:],
                         p = pdata[i,:,1,:],
                         r = rdata[i,1,:] - r0,
                         u = udata[i,:,1,:])
        r0=rdata[i,1,:]
        datelist.append("".join(dates[i]))
        data.append(outputdata)
    
    return data,datelist

def load_icar(filenames):
    filenames=glob.glob(filenames)
    filenames.sort()
    data=[]
    r0=mygis.read_nc(filenames[0],"rain").data[1,:]*0
    for f in filenames:
        r=mygis.read_nc(f,"rain").data[1,:]
        outputdata=Bunch(qv= mygis.read_nc(f,"qv").data[1,:,:],
                         qc= mygis.read_nc(f,"qc").data[1,:,:]
                            +mygis.read_nc(f,"qi").data[1,:,:]
                            +mygis.read_nc(f,"qs").data[1,:,:]
                            +mygis.read_nc(f,"qr").data[1,:,:],
                         t = mygis.read_nc(f,"th").data[1,:,:],
                         p = mygis.read_nc(f,"p").data[1,:,:],
                         r = r - r0,
                         u = mygis.read_nc(f,"u").data[1,:,:])
        r0=r
        data.append(outputdata)
    
    return data

def plot_panel(panel,d1,d0,title=""):
    """docstring for plot_panel"""
    nrows=len(d1)
    ncols=3
    currow=0
    for k in d1.keys():
        plt.subplot(nrows,ncols,currow*ncols+panel)
        if k!="r":
            delta=(d1[k]-d0[k])
            vrange=max(abs(delta.min()),delta.max()) * 0.95
            if vrange==0: 
                delta=d1[k]-d1[k].mean()
                vrange=max(abs(delta.min()),delta.max()) * 0.95
                
            plt.imshow(delta.repeat(5,axis=0),cmap=plt.cm.seismic)
            plt.clim(-vrange,vrange)
            plt.colorbar()
            if k==d1.keys()[0]:
                plt.title(title+k)
            else:
                plt.title(k)
        else:
            plt.plot(d1[k])
            plt.plot(d0[k])
            plt.title("Precip")
        currow+=1
    

def make_plots(wrf,icar,w0,i0,title=""):
    plot_panel(1,wrf,w0,title=title+" WRF : ")
    plot_panel(2,icar,i0,title="ICAR : ")
    plot_panel(3,wrf,icar,title="WRF - ICAR : ")
    plt.tight_layout()
        

def main(wrffile,icarfiles):
    """docstring for main"""
    if len(sys.argv)>1:
        icar_dir=sys.argv[1]+"/"
        if len(sys.argv)>2:
            title=sys.argv[2]
        else:
            title=""
    else:
        icar_dir="output/"
    
    print("Loading WRF data")
    wrf_data,dates=load_wrf(wrffile)
    print("Loading ICAR data")
    icar_data=load_icar(icar_dir+icarfiles)
    dt=11
    
    plt.figure(figsize=(16,13))
    for d,w,i,j in zip(dates[:24:dt],wrf_data[:24:dt],icar_data[:24:dt],range(0,24,dt)):
        print("Ploting : "+d)
        make_plots(w,i,w0=wrf_data[0],i0=icar_data[0],title="")
        plt.savefig(icar_dir+title+d+".png")
        plt.clf()
        


if __name__ == '__main__':
    main("wrfout_d01_0001-01-01_00:00:00","icar*")