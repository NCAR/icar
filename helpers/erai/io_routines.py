import subprocess,os
import numpy as np
import mygis
from bunch import Bunch


sfcvarlist=["SSHF_GDS4_SFC","SLHF_GDS4_SFC","Z_GDS4_SFC","BLH_GDS4_SFC","SSRD_GDS4_SFC","STRD_GDS4_SFC", "SSTK_GDS4_SFC", "CP_GDS4_SFC", "LSM_GDS4_SFC"]
icar_sfc_var=["sensible_heat","latent_heat","hgt_98","PBL_height","sw","lw", "tskin", "cp","landmask"]

atmvarlist=["Z_GDS4_HYBL","T_GDS4_HYBL","Q_GDS4_HYBL","LNSP_GDS4_HYBL","CLWC_GDS4_HYBL","CIWC_GDS4_HYBL","lv_HYBL2_a","lv_HYBL2_b","P0"]
icar_atm_var=["gph","t","qv","ln_p_sfc","cloud","ice","sigma_a","sigma_b","P0"]

atmuvlist=["U_GDS4_HYBL","V_GDS4_HYBL"]
icar_uv_var=["u","v"]

converted_sfc_files=[]
sfc_ncfiles=dict()

def grib2nc(erai_file,varlist,output_dir):
    """convert a grib file to a netcdf file"""
    print("Converting: "+erai_file.split("/")[-1])
    print("ncl_convert2nc "+erai_file+" -e grb -L -v "+",".join(varlist)+" -o "+output_dir)
    outputfile=output_dir+erai_file.split("/")[-1]+".nc"
    if not os.path.isfile(outputfile):
        try:
            os.system("ncl_convert2nc "+erai_file+" -e grb -L -v "+",".join(varlist)+" -o "+output_dir +"&> /dev/null")
        except:
            print("ERROR: ncl_convert2nc is not available in your $PATH (?)")

    return outputfile

def find_sfc_file(time,info):
    file_base= info.sfcdir+info.sfcfile
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    file_base= file_base.replace("_D_","{0:02}".format(time.day))
    if (time.hour>=12):
        hour="12"
    else:
        hour="00"
    file_base= file_base.replace("_h_",hour)
    offset=round(time.hour/3)
    return file_base,offset

def find_atm_file(time,info):
    file_base= info.atmdir+info.atmfile
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    file_base= file_base.replace("_D_","{0:02}".format(time.day))
    atm_file = file_base.replace("_h_","{0:02}".format(time.hour))

    file_base= info.atmdir+info.uvfile
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    file_base= file_base.replace("_D_","{0:02}".format(time.day))
    uv_file  = file_base.replace("_h_","{0:02}".format(time.hour))

    return uv_file,atm_file

def load_sfc(time,info):
    """load surface forcing from a grib file (or netcdf file if it has been converted previously)"""
    inputfile,offset=find_sfc_file(time,info)
    try:
        nc_file=sfc_ncfiles[inputfile]
    except KeyError:
        nc_file=grib2nc(inputfile,sfcvarlist,info.nc_file_dir)
        sfc_ncfiles[inputfile]=nc_file

    outputdata=Bunch()
    for s,v in zip(icar_sfc_var,sfcvarlist):
        nc_data=mygis.read_nc(nc_file,v,returnNCvar=True)
        input_data=nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
        if ((s=="latent_heat") or (s=="sensible_heat") or (s=="sw") or (s=="lw")):
            # if offset>=3:
            #     input_data[offset,...]-=input_data[offset-3,...]
            #     input_data[offset,...]/3.0
            if offset>=2:
                input_data[offset,...]-=input_data[offset-2,...]
                input_data[offset,...]/2.0
            elif offset>=1:
                input_data[offset,...]-=input_data[offset-1,...]

        if (s=="cp"):
            input_data[1:,:,:]-=input_data[:-1,:,:]

        outputdata[s]=input_data[int(offset),:,:]
        nc_data.ncfile.close()

    return outputdata

def load_atm(time,info):
    """Load atmospheric variable from a GRIB file"""
    uvfile,scfile=find_atm_file(time,info)
    uvnc_file=grib2nc(uvfile,atmuvlist,info.nc_file_dir)
    scnc_file=grib2nc(scfile,atmvarlist+["g4_lat_0","g4_lon_1"],info.nc_file_dir)

    outputdata=Bunch()
    for s,v in zip(icar_uv_var,atmuvlist):
        nc_data=mygis.read_nc(uvnc_file,v,returnNCvar=True)
        if len(nc_data.data.shape)==3:
            outputdata[s]=nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
        else:
            outputdata[s]=nc_data.data[info.ymin:info.ymax,info.xmin:info.xmax]
        nc_data.ncfile.close()

    for s,v in zip(icar_atm_var,atmvarlist):
        nc_data = mygis.read_nc(scnc_file,v,returnNCvar=True)

        if len(nc_data.data.shape)==3:
            outputdata[s] = nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
        elif len(nc_data.data.shape)==2:
            outputdata[s] = nc_data.data[info.ymin:info.ymax,info.xmin:info.xmax]
        elif len(nc_data.data.shape)==1:
            outputdata[s] = nc_data.data[:]
        else:
            try:
                outputdata[s] = nc_data.data[:]
            except:
                outputdata[s] = nc_data.data.get_value()

        nc_data.ncfile.close()

    return outputdata


def load_data(time,info):
    """docstring for load_data"""
    sfc=load_sfc(time,info)
    atm=load_atm(time,info)
    return Bunch(sfc=sfc,atm=atm)
