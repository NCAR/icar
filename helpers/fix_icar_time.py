#!/usr/bin/python
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#   ICAR helper script to fix tinestamps when an exisitng output file has been overwritten after a restart.
#   This leads to incorrect timestamps in the file, as the time:units attribute is not overwritten. As a result, 
#   one or more hours after the restart have the original 'days since' attribute, but new 'time' values, 
#   leading to weird jumps in time. 
#   Best is to remove the output file(s) before restarting (i.e. the output files of the timestep one is restarting from).
#   But in case one forgets to do so, this script can (brute force) fix it. 
#   More elegant solutions are undoubtedly possible (using option B below) , but I'll leave that to more skilled programmers
#
#   Usage:
#   > python fix_icar_times.py -i <inputfile> -o <outputfile>
#   N.B. if no output file is specified, inputfile will be overwritten. 
#
#   Note: 
#       - currently (Option A) Modifies the time:units attribute from "days since"  to "hours since". Not ideal but works.
#       - only tested for output files with hourly resolution. 
#    
#
#   Bert Kruyt,  March 2022.  NCAR.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import pandas as pd
# from datetime import datetime
import xarray as xr
# import numpy as np
import sys, getopt



def main(argv):

    """overwrite the times in an icar file based on its filename. Only tested for hourly output files"""

    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print(  'usage: fix_icar_time.py -i <inputfile> -o <outputfile>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print( 'use: fix_icar_time.py -i <inputfile> -o <outputfile>' )
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print( 'Input file is ', inputfile)
    print( 'Output file is ', outputfile)


    ####  Option A: less elegant, but more robust:   ####
    #_______  open the icar file: ______
    FIX = xr.open_dataset( inputfile )
    

    #_______ create the correct times: ______
    tstring = inputfile[inputfile.find('out_')+4:inputfile.rfind(".nc")]

    times2 = pd.date_range(tstring[:10], periods=len(FIX.time), freq='H')
    
    FIX['time'] = times2
    
    
    #_______ Write the fixed Dataset to nc file:  _________
   
    if outputfile == '':
        out_path = inputfile
    else:
        out_path = outputfile            
    
    FIX.to_netcdf( path=out_path, mode='w', encoding={'time': {'dtype': 'i4'}})  



    ###########   Option B: only modify units, but this doesnt always work as sometimes it is only the first hour, sometimes more hours.   ########
    # #_______  open the icar file: ______
    # FIX = xr.open_dataset( inputfile ,  decode_times=False)  
    
    # #_______ create the correct times: ______
    # units = FIX.time.units
    # tstring = inputfile[inputfile.find('out_')+4:inputfile.rfind(".nc")]
    # # create right time based on file name:
    # FIX['time'].attrs['units'] = units[:units.find('since')+6] + tstring[:10]
    


if __name__ == "__main__":
    main(sys.argv[1:])


