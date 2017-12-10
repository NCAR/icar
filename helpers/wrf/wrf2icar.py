#!/usr/bin/env python

"""
SYNOPSIS

    wrf2icar.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    Convert a WRF output file into a slightly smaller ICAR input file
        assumes a common set of variable names.  If these need to be changed
        e.g. for a new version of WRF, edit the list of names in wrf_vars.py
    
    Request help (-h or --help).

EXAMPLES

    wrf2icar.py "wrfout_d01_2000-10-01_00:00:00" -o icar_input_2000-10-01.nc

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    0.1
    
"""

import sys
import os
import traceback
import argparse
import numbers
# import types

from wrf_vars import wrfvars, tsvar, landmask, steps_per_day

from bunch import Bunch
import mygis
import numpy as np

global verbose
verbose=False

def load_data(filename, varname):
    """loads a variable:varname from the netcdf file:filename
    
    if varname is actually a number, then just return that number in the data structure
    """
    
    if verbose:print("Loading:"+str(varname))
    
    if isinstance(varname, numbers.Number):
        if verbose:print("   "+str(varname)+" is a number!")
        return Bunch(name="n"+str(varname), data=varname, attributes="a number", dims=(), dtype='f')
    
    # load the data as a netcdf4 data structure first
    d=mygis.read_nc(filename,varname,returnNCvar=True)
    # copy the dimensions of the variable
    dims=d.data.dimensions
    # read the attributes of the variable
    atts=mygis.read_atts(filename,varname)
    if verbose:print("    "+str(dims))
    if verbose:print("    "+str(atts))
    
    # now read the actual data
    data=d.data[:]
    d.ncfile.close()
    
    return Bunch(name=varname,data=data,dims=dims,attributes=atts, dtype='f')
    

def load_vars(filename,varnames):
    """Loads variables defined in the varnames list structure"""
    outputdata=[]
    op=None
    
    # loop through vars to load
    for v in varnames:
        # if this var is a string or a number, "load" it
        if (type(v)==str) or isinstance(v, numbers.Number):
            if op==None:
                # if there is no operator pending then just load the data
                outputdata.append(load_data(filename,v))
            else:
                # other wise, load the data, then apply the operator to this and the last data
                newdata=load_data(filename,v)
                
                if verbose:print("  applying operator"+str(op)+" to:"+outputdata[-1].name+" and "+newdata.name)
                # try just passing in the variables
                try:
                    # first try applying the operator to the data structures themselves
                    op(outputdata[-1],newdata)
                # if that doesn't work, try passing in the data only
                except:
                    outputdata[-1].data=op(outputdata[-1].data,newdata.data)
                    
                # clear the operator so we don't use it again
                op=None
                
        elif type(v)==list:
            # we this var is actually a list, call load_vars on this list
            if verbose:print("Loading a list")
            outputdata.extend(load_vars(filename,v))
            
        elif (hasattr(v,__call__)):#(type(v)==types.BuiltinFunctionType) or (type(v)==types.FunctionType):
            # if this var is actually a function than store it as an operator to process
            op=v
            if verbose:print(op)
    
    return outputdata
        
def load_tskin(filename,tsvarname, maskvarname):
    """load skin temp data from wrf output file: filename
    
    tskin uses it's own procedure so that we can average over the past 10 days for land points
    ideally it will create a higher resolution output based on a high-resolution land-sea mask
    """
    if verbose:print("Loading :"+maskvarname)
    ncmask=mygis.read_nc(filename,maskvarname,returnNCvar=True)
    mask=ncmask.data[0]
    ncmask.ncfile.close()
    land=np.where(mask==1)
    
    if verbose:print("Loading :"+tsvarname)
    
    ncts=mygis.read_nc(filename,tsvarname,returnNCvar=True)
    atts=mygis.read_atts(filename,tsvarname)
    dims=ncts.data.dimensions
    
    if verbose:print("    "+str(dims))
    if verbose:print("    "+str(atts))
    
    data=ncts.data[:]
    ncts.ncfile.close()
    
    for i in range(data.shape[0]-1,0,-1):
        start=max(i-(10*steps_per_day),0)
        data[i][land]=data[start:i+1].mean(axis=0)[land]
    
    return Bunch(name=tsvarname, data=data, attributes=atts, dims=dims, dtype='f')


def load_wrf_data(filename):
    """Load required data form the WRF output file : filename"""
    base_data=load_vars(filename,wrfvars)
    skin_t=load_tskin(filename,tsvar,landmask)
    base_data.append(skin_t)
    
    atts=mygis.read_atts(filename,global_atts=True)
    
    return Bunch(data=base_data,global_atts=atts)


def write_icar_file(filename,data):
    """write the output"""
    base_var=data.data.pop()
    if verbose:print("Writing:"+filename)
    mygis.write(filename, base_var.data,
                varname=base_var.name, dims=base_var.dims, attributes=base_var.attributes, 
                global_attributes=data.global_atts, extravars=data.data)


def main (filename, outputfile):
    """Convert WRF output:filename into ICAR input:outputfile"""
    input_data=load_wrf_data(filename)
    
    write_icar_file(outputfile,input_data)
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Simple program to trim down a WRF output file for ICAR input',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename', action='store', default="wrfout", help="WRF output filename")
        parser.add_argument('-o', dest="outputfile", action='store', default=None, help="ICAR input file to create")
        parser.add_argument('-v', '--version',action='version',
                version='wrf2icar 0.1')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        verbose=args.verbose
        
        if args.outputfile==None:
            outputfile="icar_"+args.filename
        else:
            outputfile=args.outputfile
        
        exit_code = main(args.filename, outputfile)
        
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
