#!/usr/bin/env python

"""
SYNOPSIS

    setup_next_run.py [-h] [--verbose] [-v, --version] [options_filename] [template_filename] [-s #]

DESCRIPTION

    Updates an ICAR input options file to restart a simulation.  This will look for output files
    following the output filename prefix specified in an existing options file to find older output files.
    It will then find the last time step likely to be valid in that file and restart the simulation from there.

EXAMPLES

    setup_next_run.py icar_options.nml template.nml

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    1.1 Uses argparse to improve commandline handling, add -s flag to step backwards a number of output files

"""
from __future__ import absolute_import, print_function, division

import traceback
import argparse

global verbose
verbose=False

import glob,os,re,sys, fnmatch
from math import floor

import xarray as xr

def find_last_output(options_file, skip):
    """docstring for find_last_output"""
    with open(options_file,"r") as f:
        for l in f:
            ltest = l.split("!")[0] # only look at the line before any comment characters
            if fnmatch.fnmatch(ltest.strip(),"restart_file*=*"):
                outputfile = l.split("=")[1:]
                if len(outputfile)>1:
                    outputfile = "=".join(outputfile)
                else:
                    outputfile=outputfile[0]
                outputfile=outputfile.replace('"',"").replace("'","").replace(',',"").replace('\n',"")

    outputfile = outputfile.split()[0]
    print("Found outputfile prefix: "+outputfile)
    filelist=glob.glob(outputfile+"*00.nc")
    filelist.sort()

    return filelist[-1-skip], filelist[-2-skip], outputfile

def load_last_date(filename, prefix):
    """docstring for load_last_date"""
    time_data = xr.open_dataset(filename)["time"]
    last_hour = len(time_data)

    datetime_string = filename.replace(prefix,"")
    date_string = datetime_string.split("_")[-2]
    time_string = datetime_string.split("_")[-1][:-2]
    print(date_string + "  " + time_string)
    year    = date_string.split("-")[0]
    month   = date_string.split("-")[1]
    day     = date_string.split("-")[2]

    hour = int(time_string[:2])
    # hour    = last_hour - 2 # look back two time steps in case the last time steps were corrupted
                            # (not sure why more than 1 is / can be, but it happened)

    return "{0}, {1}, {2}, {3}, 0, 0".format(year,month,day,hour),hour


def main(options_file, template_file, skip):
    """docstring for main"""
    restart_file, backup_file, output_prefix = find_last_output(options_file, skip)
    try:
        restart_date, hour = load_last_date(restart_file, output_prefix)
    except:
        restart_file = backup_file
        restart_date, hour = load_last_date(restart_file, output_prefix)

    print("Using restart file: "+restart_file)
    print("      restart date: "+restart_date)

    # escape the filename for sed
    restart_file = restart_file.replace("/","\/")
    os.system("echo  's/__RESTART_DATE__/"+str(restart_date)   +"/g;"+
                      "' "+template_file+" #>"+options_file)
    os.system("sed -e 's/__RESTART_DATE__/"+str(restart_date)   +"/g;"+
                      "' "+template_file+" >"+options_file)


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Setup an ICAR options file to continue from a past run. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('options',  nargs="?",   action='store', default="icar_options.nml",
                            help="""ICAR options file name to read and write,
                                    this file will be overwritten by an updated file
                                    based off the template""")
        parser.add_argument('template', nargs="?",   action='store', default="icar_template.nml",
                            help="template file name to read")
        parser.add_argument('-s',       dest="skip", action='store', default=0,type=int,
                            help="number of output files to step backwards (e.g. skip the last output file -s 1)")
        parser.add_argument('-v', '--version',action='version',
                version='setup_next_run 2.0')
        parser.add_argument ('--verbose', action='store_true', default=False, dest='verbose',
                help='verbose output')
        args = parser.parse_args()

        verbose=args.verbose

        verbose = args.verbose

        exit_code = main(args.options, args.template, args.skip)
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
