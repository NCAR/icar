#!/usr/bin/env python

"""
SYNOPSIS

    template_argparse.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION

    
"""

import sys
import os
import traceback
import argparse
import glob

import numpy as np
import mygis

global verbose

def update_data(curdata,lastdata):
    """docstring for update_data"""
    lost_rain=np.zeros(lastdata.shape)
    if np.max(curdata[0,5:-5,5:-5]) < np.max(lastdata[5:-5,5:-5]):
        lost_rain[:]=lastdata
        start_point=0

    for i in range(1,curdata.shape[0]):
        if np.max(curdata[i,5:-5,5:-5]) < np.max(curdata[i-1,5:-5,5:-5]):
            lost_rain[:]=curdata[i-1]
            start_point=i
    
    curdata[start_point:]+=lost_rain
    return lost_rain
    

def main (filesearch, outputfile):

    files=glob.glob(filesearch)
    files.sort()
    d0=mygis.read_nc(files[0],"rain").data
    outputdata=np.zeros((len(files), d0.shape[1],d0.shape[2]))
    
    if verbose:print("Outputdata.shape : " + str(outputdata.shape))
    last_data=d0[0]*0
    lost_rain=d0[0]*0
    
    for i,f in enumerate(files):
        if verbose:print(f)
        d=mygis.read_nc(f,"rain").data+lost_rain
        if np.max(d[-1,5:-5,5:-5]) < np.max(last_data[5:-5,5:-5]):
            if verbose: print("Updating lost_rain")
            d-=lost_rain
            lost_rain = update_data(d,last_data)
        
        outputdata[i]=d[-1]-last_data
        last_data=d[-1]
    
    if verbose:print("Writing output file")
    mygis.write(outputfile, outputdata, varname="pr")

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Create a daily precipitation file from ICAR output files. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filesearch',nargs="?",         action='store',     default="icar_*.nc",
                            help="glob pattern to search for icar output files")
        parser.add_argument('-o',        dest="outputfile", action='store',     default="daily_precip.nc",
                            help="name of outputfile")
        parser.add_argument('-v','--version',               action='version',   version='extract_daily_precip 1.0')
        parser.add_argument ('--verbose', dest='verbose',   action='store_true',default=False, help='verbose output')
        args = parser.parse_args()

        # note: verbose is defined above as global
        verbose=args.verbose
        exit_code = main(args.filesearch, args.outputfile)
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
