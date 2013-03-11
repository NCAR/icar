#!/usr/bin/env python

"""
SYNOPSIS

    swim_wrf.py [-h] [--verbose] [-v, --version] [parameter_file]

DESCRIPTION
    Reads forcing data to drive the: 
        Simple Weather (Interpolation) Model  (SWIM)
        Pseudo-Dynamical Downscaling Model (PDDM)
        Physics based Orographic Precipitation model (POP)
        Hybrid Orographic Precipitation model (HOP)

EXAMPLES

    swim_wrf.py swim_parameters.txt

EXIT STATUS

    None

WARNING/KNOWN/LIKELY ISSUES
    This may fail on windows because windows does not have a proper "fork()"
    As a result, dictionaries in classes passed to sub-processes can effectively 
    get re-initialized (or somesuch). Not sure what that will do to the Bunch class which 
    relies on an internal __dict__.  IF this is an issue, it can be fixed by removing 
    the parallel wrapper around the IO code. 

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This program is in the public domain.

VERSION
    0.6.1 - provided commandline and inputfile option handling. 
            and (slightly) streamlined main()
    0.6   - adding WRF high res winds
    0.5.1 - SPEED (esp. with new microphysics [physics=0] ~50x faster)
    0.5   - added LT winds & fancy 3D grid(?)
    0.4   - adding PBL mixing and (stupid) radiative cooling 1.5K/day
    0.3   - WRF based (previously from swim_narr.py)
    0.2   - the hey it sort of works version!
    0.1   - skeletal, pieces work
"""
# Python standard modules
import sys
import os
import traceback
import argparse # requires python >=2.7
from multiprocessing import Pool,Queue,Process
import glob
import time
import copy
import re
import ast #used to convert strings in parameter file to python variables
# "Standard" 3rd party modules
import numpy.fft as fft
from scipy.signal import convolve
import numpy as np
# my common modules
# import units
from bunch import Bunch
# swim specific modules
import swim_io # this class combines the functionality of the three file IO classes below
import swim as swim_lib #this is the meet of the fortran code
import swim_fast as swim # a faster version of the hand off between python and fortran
import fast_mean # inline C to do a running spatial mean
import lt_winds #computes Linear Theory 3D winds
from wrf_reader import WRF_Reader

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2
physics=int(1) #1=thompson other=simple


def convert_p(p,h,dz):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals
    # h,dz in meters
    slp = p/(1 - 2.25577E-5*h)**5.25588
    pout=slp*(1 - 2.25577E-5*(h+dz))**5.25588
    return pout


def topo_adjust_weather(hitopo,lowtopo,weather):
    N=hitopo.shape
    # just adjust all pressures by dz (difference between coarse topo and high res topo)
    dz=hitopo-lowtopo
    if len(dz.shape)==2:
        weather.p=convert_p(weather.p,weather.hgt,dz[np.newaxis,:,:].repeat(weather.p.shape[0],axis=0))
    else:
        weather.p=convert_p(weather.p,weather.hgt,dz)

def fix_top_bottom(topo,dz,edgevalue=None):
    if edgevalue==None:
        # find the worst dz between the top and bottom of the dataset
        worst_offset=np.max(np.abs(topo[0,:]-topo[-1,:]))
        # find the minimum distance required to smoothly transition between those
        # two endpoints
        mindist=np.round(worst_offset/(dz/10.0))
        # calculate the mid point / average topography that will be the new borders
        aves=(topo[0,:]+topo[-1,:])/2.0
    else:
        # find the worst dz between the top and bottom of the dataset
        worst_offset=np.max([np.abs(topo[0,:]-edgevalue),np.abs(topo[-1,:]-edgevalue)])
        # find the minimum distance required to smoothly transition between those
        # two endpoints
        mindist=np.round(worst_offset/(dz/2.0))*2
        # calculate the mid point / average topography that will be the new borders
        aves=np.zeros(topo.shape[1])+edgevalue # (topo[0,:]+topo[-1,:])/2.0
    # create the output topography array
    sz=topo.shape
    if mindist%2==1:
        mindist+=1
    if mindist<2:
        mindist=2
    outputtopo=np.zeros((sz[0]+mindist,sz[1]))
    # fill in everything inside the padding zone with the original topography
    padding=np.round(mindist/2.0)
    outputtopo[padding:sz[0]+padding,:]=topo
    # now calculate a fading factor (0-1) between the actual topography and the new boundary
    fade=(np.arange(padding)/np.float(padding)).reshape((padding,1)).repeat(sz[1],axis=1)
    # finally calculate the topography over the padding zones by linearly combining the 
    # aves / mid-point with the actual topography fading from 100% ave to 100% real topo
    border=aves.reshape((1,-1)).repeat(padding,axis=0)
    topo_top=topo[0,:].reshape((1,-1)).repeat(padding,axis=0)
    topo_bottom=topo[-1,:].reshape((1,-1)).repeat(padding,axis=0)

    outputtopo[:padding,:]=(1-fade)*border + fade*topo_top
    outputtopo[-padding:,:]=fade*border + (1-fade)*topo_bottom
    return outputtopo,padding

    
def topo_preprocess(topo,test=False):
    acceptable_topo_diff=np.max(np.diff(topo))/2.0
    vpad=0.0
    hpad=0.0
    # This method of finding the edge might allow a slightly smoother interpolation, but with mroe negative values
    #  on the interior of the domain.  Not sure that actually matters or not though...
    # edgevalue=(topo[0,:].mean()+topo[-1,:].mean()+topo[:,0].mean()+topo[:,-1].mean())/4.0
    
    edge=np.min([topo[0,:].min(),topo[-1,:].min(),topo[:,0].min(),topo[:,-1].min()])
    if test:
        topo=topo.copy() #dont change the topo that is outside this routine
        t1=np.arange(16.0,-1,-1.0)[np.newaxis,:].repeat(160,axis=0)/18.0
        zt=topo[91:251,:17]
        topo[91:251,:17]=1500*t1+zt*(1-t1)
    
        t1=np.arange(13.0)[:,np.newaxis].repeat(40,axis=1)/12.0
        zt=topo[250:,78:118]
        topo[250:,78:118]=2200*t1+zt*(1-t1)
    
        topo[245:,:50]=1500
        topo[topo<1500]=1500
    
    (topo,vpad)=fix_top_bottom(topo,acceptable_topo_diff,edge)
    (topo,hpad)=fix_top_bottom(topo.T,acceptable_topo_diff,edge)
    topo=topo.T
    #this might matter more than I thought... and maybe all edges should get exactly to 0?
    topo-=edge #1500# (topo[0,:].mean()+topo[-1,:].mean()+topo[:,0].mean()+topo[:,-1].mean())/4.0

    if (topo.shape[0]+30)<512:
        xtra_vpad=(512-topo.shape[0])/2
    elif (topo.shape[0]+30)<1024:
        xtra_vpad=(1024-topo.shape[0])/2
    else:
        xtra_vpad=50
    if (topo.shape[1]+30)<512:
        xtra_hpad=(512-topo.shape[1])/2
    elif (topo.shape[1]+30)<1024:
        xtra_hpad=(1024-topo.shape[1])/2
    else:
        xtra_hpad=50
    xtra_vpad=int(xtra_vpad)
    xtra_hpad=int(xtra_hpad)
    newtopo=np.zeros((topo.shape[0]+xtra_vpad*2,topo.shape[1]+xtra_hpad*2),dtype=topo.dtype)
    newtopo[xtra_vpad:xtra_vpad+topo.shape[0],xtra_hpad:xtra_hpad+topo.shape[1]]=topo
    return newtopo,xtra_vpad+vpad,xtra_hpad+hpad
    
        
def update_weather(newatm,oldatm):
    """Copies internal weather information from the old hi res atmosphere to the new one
        while maintaining the updated boundary conditions of the new atmosphere"""
    for k in newatm.keys():
        if (k!="hgt") and (k!="date") and (k!="sfc") and (k!="p") and (k!="u") and (k!="v") and (k!="w"):
            newatm[k][:,1:-1,1:-1]=oldatm[k][:,1:-1,1:-1]
            # tmp=np.where(np.isfinite(oldatm[k][:,1:-1,1:-1]))
            # newatm[k][:,1:-1,1:-1][tmp]=oldatm[k][:,1:-1,1:-1][tmp]
            # tmp=np.where(~np.isfinite(oldatm[k][:,1:-1,1:-1]))
            # if len(tmp[0])>0:
            #     print("BROKEN:"+str(oldatm.date))
            #     oldatm[k][:,1:-1,1:-1][tmp]=newatm[k][:,1:-1,1:-1][tmp]

# Force vertical convergence/divergence to balance horizontal convergence/divergence
def adjust_winds(u,v,w,prestaggered=False):
    if not prestaggered:
        ustag=(u[:,:,1:]+u[:,:,:-1])/2
        vstag=(v[:,1:,:]+v[:,:-1,:])/2
    else:
        ustag=u
        vstag=v
    divergence=(ustag[:,1:-1,1:]-ustag[:,1:-1,:-1])+(vstag[:,1:,1:-1]-vstag[:,:-1,1:-1])
    w[0,1:-1,1:-1]=-divergence[0,:,:]
    for i in range(w.shape[0]-1):
        w[i+1,1:-1,1:-1]=w[i,1:-1,1:-1]-divergence[i+1,:,:]
    return(ustag,vstag,w)
    
def rotate_winds(u,v,hgt,dx=4000.0,rotation=None):
    if rotation==None:
        urot=np.cos(np.arctan(np.abs(hgt[:,:,1:]-hgt[:,:,:-1])/dx))
        vrot=np.cos(np.arctan(np.abs(hgt[:,1:,:]-hgt[:,:-1,:])/dx))
        rotation=[urot,vrot]
    else:
        urot,vrot=rotation
    u/=urot
    v/=vrot
    return rotation

def calc_ndsq(weather,base):
    return 6.37e-5
    # return 1.0e-6
    R  = 287.0
    Rv = 461.0
    cp = 1004.0
    L   = 2.5e6
    g  = 9.81
    ratio = 18.015/28.964

    t0 = 273.15
    p0=weather.p[0,...]
    pii=1.0/((100000.0/weather.p)**(R/cp))
    T2m=weather.th[0,...]*pii[0,...]
    
    es = 611.21*np.exp(17.502*(T2m-t0)/(T2m-32.19))
    qs0 = ratio * es/(p0-es)

    cap_gamma = -(g * (1.+(L*qs0)/(R*T2m)) / (cp + (L**2 * qs0*ratio) / (R*T2m**2)))
    env_gamma = np.mean(np.diff(weather.th*pii,axis=0)/np.diff(base.hgt3d,axis=0),axis=0)
    
    dry_gamma=np.mean(env_gamma-cap_gamma)
    print(T2m.mean(),env_gamma.mean(),cap_gamma.mean(),dry_gamma)
    ndsq=(g/T2m.mean())*(dry_gamma)
    ndsq=max(min(1e-5,ndsq),1e-8)
    return ndsq
    

def simul_next(base,wrfwinds,q,r_matrix,Fzs,padx,pady,options=None,forcing=None):
    weather=base.next()
    topo_adjust_weather(base.hgt3d, weather.hgt, weather)
    
    # weather.qv/=1.05
    if options!=None:
        use_linear_winds=options.use_linear_winds
        use_wrf_winds=options.use_wrf_winds
        wrfres=options.wrfres
        if options.const_U!=None:
            weather.u=weather.u*0+options.const_U
        if options.const_V!=None:
            weather.v=weather.v*0+options.const_V
        if options.const_p!=None:
            if np.iterable(options.const_p):
                weather.p=weather.p*0+np.array(options.const_p)[:,np.newaxis,np.newaxis]
                topo_adjust_weather(base.hgt3d,0*weather.hgt+base.hgt3d[:,0,0][:,np.newaxis,np.newaxis],weather)
            else:
                weather.p=weather.p*0+options.const_p
                topo_adjust_weather(base.hgt3d,0*weather.hgt,weather)
        if options.const_qv!=None:
            weather.qv=weather.qv*0+np.array(options.const_qv)[:,np.newaxis,np.newaxis]
            weather.qc*=0
        if options.const_t!=None:
            weather.th=weather.th*0+np.array(options.const_t)[:,np.newaxis,np.newaxis]
    if forcing!=None:
        Fzs=forcing.Fzs
        padx=forcing.padx
        pady=forcing.pady
        r_matrix=forcing.r_matrix
    
    if use_linear_winds:
        (junku,junkv,weather.w)=adjust_winds(weather.u,weather.v,weather.w)
        # calculate the (squared) brunt vaisala frequency, better off with the dry value (subtract cap_gamma)
        ndsq=calc_ndsq(weather,base)
        lt_winds.update_winds(base.hgt3d,Fzs,weather.u,weather.v,weather.w,dx=wrfres,
                              Ndsq=ndsq,r_matrix=r_matrix,padx=padx,pady=pady,rotation=False)
        # weather.u=fast_mean.fast_smooth(weather.u,5)
        # weather.v=fast_mean.fast_smooth(weather.v,5)
    elif use_wrf_winds:
        print("LOADING WRF WINDS")
        wind=wrfwinds.next()
        # wind.u=fast_mean.fast_smooth(wind.u,2)
        # wind.v=fast_mean.fast_smooth(wind.v,2)
        weather.u=wind.u
        weather.v=wind.v
        # weather.p=wind.p
        # weather.th=wind.th
        # weather.qv=wind.qv

    weather.u=(weather.u[:,:,1:]+weather.u[:,:,:-1])/2
    weather.v=(weather.v[:,1:,:]+weather.v[:,:-1,:])/2
    if r_matrix !=None:
        r_matrix=rotate_winds(weather.u,weather.v,base.hgt3d,dx=wrfres,rotation=r_matrix)
    (weather.u,weather.v,weather.w)=adjust_winds(weather.u,weather.v,weather.w,prestaggered=True)
    
    q.put(weather)



class Forcing_Reader(object):
    old=None
    new=None
    
    base=None
    q=None
    reader_process=None
    # rotation matrix to rotate winds into the real domain
    r_matrix=None
    
    # specific to useing WRF winds (and other variables) directly
    wrfwinds=None

    # specific to linear winds
    Fzs=None #FFT of topography with smoothed borders
    padx=None #padding required to smooth topography
    pady=None 
    # for constant forcing, we need to store the start date
    #  so all outputfiles can have a number added to that date
    init_date=""
    curdate=0
    # values to use to make U and V fields constant in space and time
    const_U=None
    const_V=None
    
    def __init__(self,options,domain):
        self.q=Queue()
        self.base=WRF_Reader(options.file_search,sfc=options.sfc_file_search,
                            bilin=True,nn=False,domain=domain,options=options)
        
        self.constant_forcing=options.constant_forcing
        
        if options.use_linear_winds:
            (self.lt_topo,self.pady,self.padx)=topo_preprocess(domain.topo)
            Ny,Nx=domain.topo.shape
            self.Fzs=fft.fftshift(fft.fft2(self.lt_topo))/((Nx+self.padx*2)*(Ny+self.pady*2))
            
        if options.use_wrf_winds:
            self.wrfwinds=WRF_Reader(options.wind_files,bilin=False,nn=True, 
                                     windonly=True,domain=domain,options=options)
            self.base.hgt3d=self.wrfwinds.hgt3d #if we are using winds from wrf, we need to use the 3d domain from wrf too
            
        self.reader_process=Process(target=simul_next,
                args=(self.base,self.wrfwinds,self.q,self.r_matrix,self.Fzs,self.padx,self.pady,options))
        self.reader_process.start()
        self.new=self.q.get()
        if self.constant_forcing:
            self.old=self.new
        # if base is copied... can I time_inc() and spawn off the next process before the last one finishes?
        self.init_date=self.new.date
        if not self.constant_forcing:
            self.base.time_inc()
            self.reader_process=Process(target=simul_next,
                    args=(self.base,self.wrfwinds,self.q,self.r_matrix,self.Fzs,self.padx,self.pady,options))
            self.reader_process.start()

    def next(self):
        self.old=self.new
        if self.constant_forcing:
            self.old.date=self.init_date[:-6]+"_{0:05d}".format(self.curdate)+self.init_date[-6:]
            self.curdate+=1
            return Bunch(old=self.old,new=self.new,deltas=None)
        self.new=self.q.get()
        self.base.time_inc() #note this seems to be necessary because putting it in its own thread copies the reader object(?)

        self.reader_process=Process(target=simul_next,
                args=(self.base,self.wrfwinds,self.q,self.r_matrix,self.Fzs,self.padx,self.pady,options))
        self.reader_process.start()

        update_weather(self.new,self.old)
         
        # deltas=self.calc_deltas(self.old,self.new)
        deltas=None
        return Bunch(old=self.old,new=self.new,deltas=deltas)

    def __iter__(self):
        return self
    def close(self):
        if self.base!=None:
            self.base.close()
        if self.wrfwinds!=None:
            self.wrwinds.close()
    def __enter__(self):
        return self
    def __exit__(self):
        self.close()
    def __del__(self):
        self.close()
        # super(Forcing_Reader,self).__del__()



def parallel_output(filename,varname,data,curdate):
    """Write an outputfile, self-contained so it can be called in its own thread"""
    outputcurdate=''.join(''.join(''.join(curdate.split('-')).split('_')).split(':'))
    sz=data.shape
    if len(sz)==3:
        Nz,Ny,Nx=sz
        ncout=swim_io.NC_writer(filename,Nx,Ny,Nz,var=varname)
    else:
        Ny,Nx=sz
        ncout=swim_io.NC_writer(filename,Nx,Ny,var=varname)
    ncout.appendToVar(data, date=outputcurdate)
    ncout.close()
    
def write_output(outputlist,weather,options):
    """Write output files in independant threads.
    
    Loops over output specified in the options and spawns a new thread for each. 
    Because this is parallelized, .copy() variables into new arrays to prevent the next input possibly overwriting these. 
    (NOTE I'm pretty sure this shouldn't actually happen, so we could try removing that copy,
    but then you are technically in a race condition. )"""
    for i,v in enumerate(options.output):
        if i>len(outputlist):
            outputlist.append(None)
        else:
            if outputlist[i]:
                outputlist[i].join()
        outputlist[i]=Process(target=parallel_output,
                              args=(options.output_base+"swim_"+v+"_"+str(weather.old.date)[:-6],
                                    v,weather.new[v].copy(),str(weather.old.date[:-6])))
        outputlist[i].start()

def setup_domain(options):
    """Get topography info"""
    try:
        topovarname=options.topovarname
    except KeyError,AttributeError:
        topovarname="HGT"
    sub_y0,sub_y1,sub_x0,sub_x1=options.subset
    halfk=options.halfk    
    topoinfo=swim_io.read_nc(options.topofile,var=topovarname)
    if sub_x1==None:
        xendpt=-halfk
    else:
        xendpt=sub_x1-halfk
    if sub_y1==None:
        yendpt=-halfk
    else:
        yendpt=sub_y1-halfk
    topo=topoinfo.data[0,sub_y0+halfk:yendpt,sub_x0+halfk:xendpt]
    # (Ny,Nx)=topo.shape
    return Bunch(topo=topo)


def write_parameter_file(options):
    with open(options.output_base+"swim_parameters.txt",'w') as f:
        for k in options.keys():
            f.write(k+" = "+repr(options[k])+"\n")
            # note in the old (non-ast.literal_eval) style, 
            # the characters: (,),[,],",', etc. must all be replaced with " "
        
def initialize(options):
    """Initialization"""
    write_parameter_file(options)
    # initialze physics engine as necessary
    swim_lib.swim_step.init(options.physics)
    # remove old output files if specified
    if options.clearold:
        # find old output files
        oldfiles=glob.glob(options.output_base+"*.nc")
        # loop over files, deleting them
        for o in oldfiles:
            os.remove(o)

    # set up an empty list for the output processes
    outputlist=[None]*len(options.output)
    return outputlist

def main(options):
    """Main program: initialization,setup domain/forcing, and loop over input data"""
    # domain could/should include linear wind topography?
    domain=setup_domain(options)
    
    # forcing object supplies methods and datastructures to loop over all forcing
    # data contained in files specified by options. 
    forcing=Forcing_Reader(options,domain)

    # any initialization that has to occur and returns a list for output process management
    #  perform after setting up forcing so initialization of e.g. microphysics can happen
    #  in parallel with the first/second forcing read that is set up in the the forcing reader init
    outputlist=initialize(options)
    
    # timeing information
    t0=time.time()
    # loop over all boundary conditions
    for weather in forcing:
        print(" ")
        print("-------------------------------------------------")
        print(weather.old.date)
        print("-------------------------------------------------")
        t1=time.time()
        # set up data for and call Fortran physics time stepping library (in swim_lib)
        weather.new.precip=swim.swim2d(domain,weather,swim_lib,options)
        t2=time.time()
        write_output(outputlist,weather,options)
        t3=time.time()
        if options.verbose:
            print("Finished Timestep:"+str(weather.old.date))
            print("Total Time:"+str(t3-t0))
            print("     input :"+str(t1-t0))
            print("     phyics:"+str(t2-t1))
            print("     output:"+str(t3-t2))
        t0=time.time()
    forcing.close()
    print("Exiting")

def default_options():
    """Return a structure of default options"""
    return Bunch(constant_forcing=False,
                 const_U =None,
                 const_V =None,
                 const_qv=None,
                 const_t =None,
                 const_p =None,
                 forcing_dir="forcing",
                 file_search="forcing/wrfout_d01_200*00",
                 wind_files="winds/wrfout_d01_200*00",
                 sfc_file_search="forcing/wrfout_d01_200*00",
                 domain_dir="baseline_info", 
                 topofile="baseline_info/4km_wrf_output.nc",
                 topovarname="HGT",
                 baseline_file="baseline_info/4km_wrf_output.nc",
                 use_linear_winds=False, 
                 use_wrf_winds=False,
                 wrfres=4000.0,
                 io_ratio=36.0/4.0, #ratio of low res driver to high res model grids
                 nlevels=20,
                 halfk=8,
                 subset=(30,-30,30,-30),
                 timestep=1.0*60*60, #time step between forcing changes [seconds] (1hr)
                 start_position=0,
                 verbose=True,
                 physics=int(0),
                 clearold=True,
                 output_base="output/",
                 output=["precip","th","qv","qc","p","qi","qs","qr","u","v","w"])
                
    
def read_options_file(filename):
    """Simple procedure to read through a file for key-value pairs
    
    lines should be key=value
    white space is stripped
    Anything after a # is treated as a comment. 
    """
    output=default_options()
    with open(filename,'r') as f:
        for l in f:
            l=l.split("#")[0]
            if re.match(".*=.*",l):
                inputlist=l.split("=")
                k=inputlist[0]
                v="=".join(inputlist[1:])
                output[k.strip()]=ast.literal_eval(v.strip())
            elif re.match(".*:.*",l):
                inputlist=l.split("=")
                k=inputlist[0]
                v="=".join(inputlist[1:])
                output[k.strip()]=ast.literal_eval(v.strip())
    return output
    

def convert_iterable(inputdata, default):
    if type(inputdata)==str:
        temp_value=inputdata.split(",")
        if len(temp_value)==1:
            temp_value=inputdata.strip().split()
        outputlist=[]
        default_type=type(default[0])
        for cur in temp_value:
            cur=cur.strip()
            if (type(default[0])!=str) and (np.iterable(default[0])):
                outputlist.append(convert_iterable(cur,default[0]))
            else:
                outputlist.append(default_type(cur))
        if type(default)==np.ndarray:
            # for some readon ndarray's new/init doesn't convert an iterable, it treats it as the shape...
            output=np.array(outputlist,dtype=default.dtype)
        else:
            output=type(default)(outputlist)
    else:
        output=inputdata
    return output

    
def setup_options(args):
    """Set up options for this model run. 
    
    Default options are set up by default_options()
    If args.inputfilename exists, options are read from file as key value pairs, overriding defaults
    if options are specified on the commandline, they override values in the input file
    """
    doptions=default_options()
    
    if args.inputfile==None:
        #note simple copy does not work obj.copy returns a dict and copy.copy() 
        # doesn't duplicate internal elements so the type casting below will not work properly
        options=copy.deepcopy(doptions) 
    else:
        options=read_options_file(args.inputfile)
    
    argdict=args.__dict__
    for k in argdict.keys():
        if k!="inputfile":
            if argdict[k]!=None:
                print(k,argdict[k])
                options[k]=ast.literal_eval(argdict[k].strip())
    
    # NOTE THE USE of ast.literal_eval makes the following unnecessary, left in for now. 
    # if for some reason ast doesn't get the type right (e.g. int instead of float?)
    #  this code could still be useful, albeit it is not quite as generic as ast. 
    # in case we read these options in from either an input file or the commandline
    #  we need to make sure the types are what the default options sets up
    # for k in doptions.keys():
    #     if (type(doptions[k])!=str) and (np.iterable(doptions[k])):
    #         # if the default option is actually an iterable (e.g. tuple,list) recursively convert 
    #         # all elements in it (i.e. lists of lists are allowed, but all sub-elements must be lists)
    #         options[k]=convert_iterable(options[k],doptions[k])
    #     else:
    #         if type(doptions[k])==type(True):
    #             if str(options[k]).strip().lower()=="false":
    #                 options[k]=False
    #             else:
    #                 options[k]=True
    #         else:
    #             options[k]=type(doptions[k])(options[k])
    return options
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='A high resolution Simple Weather (interpolation) Model (SWIM)')
        parser.add_argument('inputfile', default=None,nargs="?",action='store',
                            help="Input file specifying all other options, specified command line options will override")
        parser.add_argument('forcing_dir', default=None,nargs="?",action='store',
                            help="Directory containing coarse resolution WRF model data")
        parser.add_argument('domain_dir', default=None,nargs="?",action='store',
                            help="Directory containing high resolution model domain information \n"
                            +"(as from a single high resolution WRF input or output file)")
        parser.add_argument('-v', '--version',action='version',version='SWIM v0.6.1')
        parser.add_argument ('--verbose', action='store_true',
                default=None, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        options=setup_options(args)
        exit_code = main(options)

        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
