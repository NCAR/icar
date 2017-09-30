#!/usr/bin/env python

"""
SYNOPSIS

    plot_icar.py [-h] [--verbose] [-v, --version] [filename] [-o output_filename.png]

DESCRIPTION

    Given an ICAR input file, generate a map of the variables specified (e.g. rain)
    The output will be saved as an image (or vector) file as optionally specified 
    with -o outputfilename.png.  The output file type is defined by the extension 
    specified for the output filename. 
    
    If the user requests help (-h or --help) usage will be printed.
    If the user requests verbose output (--verbose), output will be printed at runtime

EXAMPLES

    plot_icar.py output/icar_out_2002_12_31_00-00.nc --verbose -o output_map.png

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    1.0     2016-05-31  edg     Initial Version

"""
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse

import matplotlib.pyplot as plt

import mygis
import map_vis
import custom_cmap

global verbose
verbose=False

# vars_to_map = ["rain"] #,"ta2m"] # more variables can be added to this list, this could also be made a commandline argument

def main (filename, vars_to_map, outputfile, cmin, cmax):
    
    #  reads latitude and longitude information
    if verbose: print("Reading lat/lon info")
    geo = mygis.read_geo(filename)
    map_bounds = [geo.lat.min(), geo.lat.max(), geo.lon.min(), geo.lon.max()]

    # use a subset of the reversed rainbow colormap (if not specified it will use "jet")
    cmap = custom_cmap.subset(cmap=plt.cm.rainbow_r,clim=(50,255))
    # alternatively a list of colormaps could be supplied so that each variable gets its own colormap

    if verbose: print("Looping over variables : "+str(vars_to_map))
    # optionally set the figure size with e.g. plt.figure(figsize=)
    for v in vars_to_map:
        
        # optionally put multiple plots on one figure with e.g. plt.subplot(nrows, ncols, current_plot_number)
        
        if verbose: print('  Reading variable: {} from file: {}'.format(v, filename))
        # read the last time entry [-1] in the file
        data = mygis.read_nc(filename,v).data[-1]
        if len(data.shape)>2:
            # if this is a 3D variable, just pull the first level
            data = data[0]
        
        # set up the output filename
        if len(vars_to_map)>1:
            # create a new output file by appending the variable name before the file extension
            ofile = outputfile[:-4] + "_"+v + outputfile[-4:]
        else:
            ofile = outputfile
        
        # to make a simple image with a colorbar this is all you need:
        # if verbose: print("  Plotting data")
        # plt.imshow(data, cmap=cmap)
        # plt.colorbar()
        # plt.clim(cmin, cmax) # to set the limits on the colorscale
        # plt.title(v)
        
        # To make a map instead use this code.  
        # Note that a lot of these keywords are optional, I've included them here
        #  to give a sense of what you can do with this routine : 
        if verbose: print("  Mapping data")
        m = map_vis.vis(data,geo=map_bounds,title=v,
                cmap=cmap, colorbar=True,latstep=2.0,lonstep=5.0,m=None,lat=geo.lat,lon=geo.lon,clim=(cmin,cmax),
                reproject=True,width=None,height=None,cbar_label=None, latlabels=[1,0,0,0], lonlabels=[0,0,0,1],
                statewidth=1.5, countrywidth=1.5, countrycolor="black", coastwidth=1.5, riverwidth=0.5,rivercolor="blue")

        # Stuffs more on a page, but can have a tendency to push labels off the page if the figure is not wide/tall enough
        # plt.tight_layout()
        if verbose: print("  Writing output file:"+ofile)
        plt.savefig(ofile)

        # in case we are looping over variables, we need to clear the figure for the next plot
        plt.clf()
        
        
            

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',nargs="?", action='store', default="output/*????_??_??_??-??.nc",
                            help = "Specifies the ICAR output file to read.")
        parser.add_argument('-v', dest="var", action='store', default="rain",
                            help = "Specifies the variable name to plot. ")
        parser.add_argument('-o', dest="output", action='store', default="output_map.png",
                            help = "Specifies the output plot file to create. ")
        parser.add_argument('-cmin', dest="cmin", action='store', default=0, type=float,
                            help = "Minimum value for the colorscale. ")
        parser.add_argument('-cmax', dest="cmax", action='store', default=2000, type=float,
                            help = "Maximum value for the colorscale. ")
        parser.add_argument('--version', action='version', version='Verion : 1.0')
        parser.add_argument('--verbose', action='store_true', default=False,
                            help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose
        
        #  could specify a list here
        vars_to_map = [args.var]
        
        exit_code = main(args.filename, vars_to_map, args.output, args.cmin, args.cmax)
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
