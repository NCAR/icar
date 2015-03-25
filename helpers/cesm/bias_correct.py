#!/usr/bin/env python
import sys
import glob
import time
from multiprocessing import Pool

import numpy as np
import matplotlib.pyplot as plt


from bunch import Bunch
from atm import vertical_interp
import units
import mygis

MAX_NUMBER_PROCESSES=16
# file_search="cesm_*_*_1990*.nc"
# file_search="cesm_*_*_*.nc"
file_search="cesm_002_20TR_*.nc"

days_per_month=[31,28,31,30,31,30,31,31,30,31,30,31]
start_day_per_month=[0]
start_day_per_month.extend(np.cumsum(days_per_month))

times_per_day=4
days_per_year=365 #no leap calendar

zoutputfile="month{0:02}_mean_z.nc"
qoutputfile="month{0:02}_mean_{1}.nc"

def compute_mean_z():
    """docstring for compute_mean_z"""
    print("Computing mean levels")
    print("loading files:")
    print(glob.glob(file_search))
    z=mygis.read_files(file_search,"z",axis=0)
    mygis.write("annual_mean_z.nc",z.mean(axis=0),varname="z")
    
    for month in range(12):
        print("Month "+str(month+1))
        curz=z[0]*0
        nyears=0
        month_start=start_day_per_month[month]
        start_point=month_start*times_per_day
        month_end=start_day_per_month[month+1]
        end_point=month_end*times_per_day
        
        while start_point<z.shape[0]:
            curz+=z[start_point:end_point].mean(axis=0)
            nyears+=1
            start_point+=times_per_day*days_per_year
            end_point+=times_per_day*days_per_year
            
        print(nyears)
        curz/=nyears
        mygis.write(zoutputfile.format(month+1),curz,varname="z")

def write_interpolated_6hrly(data,z,file_name=None):
    """docstring for write_interpolated_6hrly"""
    outputq_file="interpolated_{}.nc"
    if file_name!=None:
        outputq_file="{}_"+file_name
        
    for k in data.keys():
        if k!="z":
            mygis.write(outputq_file.format(k),data[k],varname=k)

def compute_rh(data,temperature):
    """docstring for compute_rh"""
    rh=units.sh2rh(temperature,data.p,data.qv)
    return rh

def vinterp_q(file_name=None):
    """interpolate all variables to a common z level"""
    print("Interpolating quantities")
    if file_name!=None:
        file_search=file_name
    
    qvarlist=["qv","theta","p","u","v"]
    qdata=Bunch()
    print("Loading data")
    t0=time.time()
    for varname in qvarlist:
        print(varname)
        qdata[varname]=mygis.read_files(file_search,varname,axis=0)
    
    print("z")
    qdata["z"]=mygis.read_files(file_search,"z",axis=0)
    mean_z=mygis.read_nc("annual_mean_z.nc","z").data
    print("  Reading data took {0:6.2f} seconds".format((time.time()-t0)))
    
    t0=time.time()
    temperature=qdata.theta*units.exner(qdata.p)
    print("  Temperature conversion took {0:6.2f} seconds".format((time.time()-t0)))
    t0=time.time()
    print("Computing RH")
    qdata["rh"]=compute_rh(qdata,temperature)
    qvarlist.append("rh")
    print("  RH calculation took {0:6.2f} seconds".format((time.time()-t0)))
    
    print("Interpolating")
    qout=Bunch()
    for varname in qvarlist:
        qout[varname]=np.zeros(qdata[varname].shape)
    
    t0=time.time()
    nsteps=qdata.z.shape[0]
    for i in range(nsteps):
        vertical_interp.interp_multivar(qdata,qout,i, qdata.z[i], mean_z,vartype=varname,
                                        inputt=temperature[i],inputp=qdata.p[i])
    print("  Interpolation took {0:6.2f} seconds per timestep".format((time.time()-t0)/nsteps))
                
    print("Writing data")
    write_interpolated_6hrly(qout,mean_z,file_name)

def compute_mean_q(wind_option):

    print("Computing Monthly Mean fields")
    qvarlist=["qv","theta","p","u","v","rh"]
    if (wind_option == "nowind"):
        qvarlist=["qv","theta","p","rh"]
        
    full_data=[]
    for varname in qvarlist:
        print(varname)
        data=mygis.read_files(varname+"_cesm_*_*.nc",varname,axis=0)
        for month in range(12):
            print("Month "+str(month+1))
            meanq=np.zeros(data.shape[1:])
            nyears=0

            month_start=start_day_per_month[month]
            start_point=month_start*times_per_day
            month_end=start_day_per_month[month+1]
            end_point=month_end*times_per_day

            while start_point<data.shape[0]:
                print(start_point,end_point)
                meanq+=data[start_point:end_point,...].mean(axis=0)
                nyears+=1
                start_point+=times_per_day*days_per_year
                end_point+=times_per_day*days_per_year

            if nyears>0:
                meanq/=nyears
                print(nyears)
                mygis.write(qoutputfile.format(month+1,varname),meanq,varname=varname)

def load_erai_means(wind_option):
    """docstring for load_erai_means"""
    eraid="erai/"
    varlist=["p","rh","ta","ua","va","z"]
    if (wind_option=="nowind"):
        varlist=["p","rh","ta","z"]
        
    outputdata=[]
    month_mid_point_doy=(start_day_per_month[1:]+np.array(start_day_per_month[:-1]))*0.5
    
    for month in range(1,13):
        curoutput=Bunch(doy=month_mid_point_doy[month-1])
        for v in varlist:
            if v=="ta":
                ta=mygis.read_nc(eraid+"regridded_ERAi_to_cesm_month{0:02}.nc".format(month),v).data
            else:
                curoutput[v]=mygis.read_nc(eraid+"regridded_ERAi_to_cesm_month{0:02}.nc".format(month),v).data
        curoutput["theta"] = ta / units.exner(curoutput.p)
        
        # erai is probably "upside down" so reverse the vertical dimension of all variables
        if curoutput.p[0,0,0]<curoutput.p[-1,0,0]:
            for v in curoutput.keys():
                if v!="doy":
                    curoutput[v]=curoutput[v][::-1,:,:]
        outputdata.append(curoutput)
    return outputdata

def load_cesm_means(wind_option):
    """docstring for load_cesm_means"""
    cesmd="means/"
    varlist=["p","rh","theta","u","v"]
    if wind_option=="nowind":
        varlist=["p","rh","theta"]
        
    outputdata=[]
    month_mid_point_doy=(start_day_per_month[1:]+np.array(start_day_per_month[:-1]))*0.5
    
    for month in range(1,13):
        curoutput=Bunch(doy=month_mid_point_doy[month-1])
        try:
            for v in varlist:
                curoutput[v]=mygis.read_nc(cesmd+"month{0:02}_mean_{1}.nc".format(month,v),v).data
            curoutput["z"]=mygis.read_nc(cesmd+"annual_mean_z.nc","z").data
        except:
            for v in varlist:
                curoutput[v]=mygis.read_nc("month{0:02}_mean_{1}.nc".format(month,v),v).data
            curoutput["z"]=mygis.read_nc("annual_mean_z.nc","z").data
            
        if curoutput["p"][0,0,0]<1300:
           curoutput["p"][0,0,0]*=100.0 # convert hPa to Pa 
        
        if (curoutput["z"][-1,0,0]<10000) and curoutput["p"][-1,0,0]<10000: # can't have 100hPa at 10km height
            print("WARNING!  Assuming that cesm z has been erroneously divided by 9.8")
            curoutput["z"]*=9.8
        
        outputdata.append(curoutput)
    
    return outputdata

def interp_era_to_cesm(erai,cesm):
    """docstring for interp_era_to_cesm"""
    print("Monthly vertical interpolation")
    varlist=erai[0].keys() #["p","rh","theta","ua","va"]
    mean_z=cesm[0].z
    inputp=np.zeros(erai[0].p.shape)
    for i in range(len(erai)):
        print("Month {}".format(i+1))
        inputp[:]=erai[i].p[:]
        temperature=erai[i].theta * units.exner(inputp)
        for v in varlist:
            if (v!="doy") and (v!="z"):
                erai[i][v]=vertical_interp.interp(erai[i][v], erai[i]["z"], mean_z,
                                                  vartype=v,  inputt=temperature, inputp=inputp)
        erai[i]["z"]=mean_z
    
    return erai

def convert_to_erai_name(varname):
    if varname=="u":
        return "ua"
    if varname=="v":
        return "va"
    else:
        return varname

def compute_biases(erai,cesm):
    """docstring for compute_biases"""
    biases=[]
    for e,c in zip(erai,cesm):
        curbias=Bunch(doy=e.doy)
        for v in c.keys():
            if (v!="doy") and (v!="z"):
                curbias[v]=c[v]-e[convert_to_erai_name(v)]
        curbias["z"]=c["z"]
        biases.append(curbias)
    return biases
    
def interpolate_monthly_to_daily(biases):
    """docstring for interpolate_monthly_to_daily"""
    daily=[]
    current=0
    for i in range(365):
        curbias=Bunch(doy=i)
        while (biases[current].doy<i):
            current+=1
            if current>=len(biases):
                current=len(biases)
                break
        
        if current==0:
            nextbias=0
            lastbias=-1
            weight=((i+365) - biases[lastbias].doy) / (365 - biases[lastbias].doy + biases[nextbias].doy)
        elif current==len(biases):
            current-=1
            nextbias=0
            lastbias=-1
            weight=(i - biases[lastbias].doy) / (365 - biases[lastbias].doy + biases[nextbias].doy)
        else:
            nextbias=current
            lastbias=current-1
            weight=(i - biases[lastbias].doy) / (biases[nextbias].doy - biases[lastbias].doy)
            
        for v in biases[current].keys():
            if (v!="doy") and (v!="z"):
                curbias[v]=biases[nextbias][v]*weight + biases[lastbias][v]*(1-weight)
        curbias["z"]=biases[nextbias]["z"]
        daily.append(curbias)
        
    daily.append(curbias) # duplicate the last day to prevent wrap around issues
    return daily
    

def load_biases(wind_option):
    """docstring for load_biases"""
    erai_data=load_erai_means(wind_option)
    cesm_data=load_cesm_means(wind_option)
    erai_data=interp_era_to_cesm(erai_data,cesm_data)
    monthly_biases=compute_biases(erai_data,cesm_data)
    daily_biases=interpolate_monthly_to_daily(monthly_biases)
    
    return daily_biases

def apply_bias_correction(wind_option):
    """docstring for  apply_bias_correction"""
    varlist=["theta","p","u","v","rh"]
    if wind_option=="nowind": print("Skipping U and V bias correction.")
    biases=load_biases(wind_option)
    current_files=glob.glob("cesm_*.nc")
    current_files.sort()
    for f in current_files:
        print("Bias Correcting : "+f)
        output_data=mygis.Dataset(f,mode="a")
        for v in varlist:
            print("  Variable : "+v)
            try:
                d=mygis.read_nc("vinterpolated/"+v+"_"+f,v).data
            except:
                d=mygis.read_nc(v+"_"+f,v).data
            
            if (wind_option=="nowind") and ((v=="u") or (v=="v")):
                pass
            if (wind_option=="meanwind") and ((v=="u") or (v=="v")):
                for i in range(d.shape[0]):
                    windbias=biases[int(np.round(i/float(times_per_day)))][v]
                    for j in range(d.shape[1]):
                        d[i,j,...]-=np.mean(windbias[j,...])
            else:
                for i in range(d.shape[0]):
                    d[i,...]-=biases[int(np.round(i/float(times_per_day)))][v]
                
            
            if v=="rh":
                d[d<1e-10]=1e-10
                d[d>100]=100
                rh=d
            else:
                output_data.variables[v][:]=d
        
        p=output_data.variables["p"][:]/100.0 # convert to mb
        t=output_data.variables["theta"][:] * units.exner(p)
        output_data.variables["qv"][:]=units.rh2sh(t,p,rh)
        output_data.variables["z"]=biases[0]["z"][np.newaxis,:,:,:]
        output_data.close()
    

def mean_z_available():
    """docstring for mean_z_available"""
    if glob.glob("annual_mean_z.nc"): 
        return True
    return False

def interp_q_available():
    """docstring for interp_q_available"""
    basefile=glob.glob(file_search)[0]
    if glob.glob("p_"+basefile): 
        return True
    elif glob.glob("vinterpolated/p_"+basefile):
        return True
    return False

def mean_q_available():
    """docstring for mean_q_available"""
    if glob.glob(qoutputfile.format(1,"p")): 
        return True
    elif glob.glob("means/"+qoutputfile.format(1,"p")):
        return True
    return False

def main(wind_option):
    """docstring for main"""
    if not mean_z_available():
        compute_mean_z()
    if not interp_q_available():
        print("WARNING, vertical interpolation is slow, parallelizing over CESM files")
        files=glob.glob(file_search)
        process_pool=Pool(min(MAX_NUMBER_PROCESSES,len(files)))
        process_pool.map(vinterp_q,files)
    if not mean_q_available():
        compute_mean_q(wind_option)
    
    apply_bias_correction(wind_option)

if __name__ == '__main__':
    print(sys.argv)
    if len(sys.argv)>2:
        cesm_file=sys.argv[2]
    else:
        cesm_file=None 
    if len(sys.argv)>1:
        option=sys.argv[1]
    else:
        option=""   
    if len(sys.argv)>2:
        wind_option = sys.argv[2]
    else:
        wind_option = ""
    
    if option == "z":
        compute_mean_z()
    elif option == "iq":
        vinterp_q(cesm_file)
    elif option == "q":
        compute_mean_q(wind_option)
    elif option == "bc":
        apply_bias_correction(wind_option)
    elif option == "dryrun":
        print("Will create mean z:" + str(not mean_z_available()) )
        print("Will interpolate q:" + str(not interp_q_available()) )
        print("Will create mean q:" + str(not mean_q_available()) )
    else:
        main(wind_option)
