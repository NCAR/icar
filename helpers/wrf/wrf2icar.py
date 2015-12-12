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
import types

from wrf_vars import wrfvars, tsvar, landmask

from bunch import Bunch
import mygis
import numpy as np

global verbose
verbose=False

def load_data(filename, varname):
    """docstring for load_data"""
    
    if verbose:print("Loading:"+str(varname))
    
    if isinstance(varname, numbers.Number):
        if verbose:print("   "+str(varname)+" is a number!")
        return Bunch(name="n"+str(varname), data=varname, attributes="a number", dims=(), dtype='f')
    
    d=mygis.read_nc(filename,varname,returnNCvar=True)
    dims=d.data.dimensions
    atts=mygis.read_atts(filename,varname)
    if verbose:print("    "+str(dims))
    if verbose:print("    "+str(atts))
    
    data=d.data[:]
    d.ncfile.close()
    
    return Bunch(name=varname,data=data,dims=dims,attributes=atts, dtype='f')
    

def load_vars(filename,varnames):
    """docstring for load_vars(filename,varnames)"""
    outputdata=[]
    op=None
    
    for v in varnames:
        if (type(v)==str) or isinstance(v, numbers.Number):
            if op==None:
                outputdata.append(load_data(filename,v))
            else:
                newdata=load_data(filename,v)
                
                if verbose:print("  applying operator"+str(op)+" to:"+outputdata[-1].name+" and "+newdata.name)
                # try just passing in the variables
                try:
                    op(outputdata[-1],newdata)
                # if that doesn't work, try passing in the data only
                except:
                    outputdata[-1].data=op(outputdata[-1].data,newdata.data)
                op=None
        elif type(v)==list:
            if verbose:print("Loading a list")
            outputdata.extend(load_vars(filename,v))
        elif (type(v)==types.BuiltinFunctionType) or (type(v)==types.FunctionType):
            op=v
            if verbose:print(op)
    
    return outputdata
        
def load_tskin(filename,tsvarname, maskvarname):
    """docstring for load_tskin"""
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
        start=max(i-10,0)
        data[i][land]=data[start:i+1].mean(axis=0)[land]
    
    return Bunch(name=tsvarname, data=data, attributes=atts, dims=dims, dtype='f')


def load_wrf_data(filename):
    """docstring for load_wrf_data"""
    base_data=load_vars(filename,wrfvars)
    skin_t=load_tskin(filename,tsvar,landmask)
    base_data.append(skin_t)
    
    atts=mygis.read_atts(filename,global_atts=True)
    
    return Bunch(data=base_data,global_atts=atts)


def write_icar_file(filename,data):
    base_var=data.data.pop()
    if verbose:print("Writing:"+filename)
    mygis.write(filename, base_var.data,
                varname=base_var.name, dims=base_var.dims, attributes=base_var.attributes, 
                global_attributes=data.global_atts, extravars=data.data)


def main (filename, outputfile):

    input_data=load_wrf_data(filename)
    write_icar_file(outputfile,input_data)
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename', action='store', default="some_file_name")
        parser.add_argument('-o', dest="outputfile", action='store', default=None)
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
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
