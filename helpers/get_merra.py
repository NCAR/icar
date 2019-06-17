#!/usr/bin/env python
"""
SYNOPSIS

    get_merra.py [-h] [--verbose] [-v, --version] -s start_date -e end_date [North] [South] [East] [West]

DESCRIPTION

    Setup a shell script that downloads and prepares MERRA data to run ICAR

    For this shell script to work, wget and nco must be installed.  This script will generate a shall script that uses wget to download the data
    it then uses nco to flip the z dimension of MERRA data and concatenate the files.
    wget needs a .netrc file setup to authenticate the NASA download server, e.g. : machine urs.earthdata.nasa.gov login <username> password <password>
    See : https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget for more information


EXAMPLES

    get_merra.py -s 2010-01-01 -e 2011-12-31

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    1.0

"""
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse
import datetime

global verbose
verbose=False

def main (start_date, end_date, north, south, east, west):

    from datetime import timedelta

    wget_exe = 'wget --auth-no-challenge=on --content-disposition "'

    # south = "21.962"
    # west = "-131.133"
    # north = "55.009"
    # east = "-64.336"
    #
    # start_date = date(2010,1,1)
    # end_date = date(2011,1,1)

    URL="https://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2I3NVASM.5.12.4%2F{0.year:04}%2F{0.month:02}%2FMERRA2_{version}.inst3_3d_asm_Nv.{0.year:04}{0.month:02}{0.day:02}.nc4&FORMAT=bmM0Lw&BBOX={south}%2C{west}%2C{north}%2C-{east}&LABEL=MERRA2_{version}.inst3_3d_asm_Nv.{0.year:04}{0.month:02}{0.day:02}.SUB.nc&SHORTNAME=M2I3NVASM&SERVICE=SUBSET_MERRA2&LAYERS=LAYER_35%2C36%2C37%2C38%2C39%2C40%2C41%2C42%2C43%2C44%2C45%2C46%2C47%2C48%2C49%2C50%2C51%2C52%2C53%2C54%2C55%2C56%2C57%2C58%2C59%2C60%2C61%2C62%2C63%2C64%2C65%2C66%2C67%2C68%2C69%2C70%2C71%2C72&VERSION=1.02&DATASET_VERSION=5.12.4&VARIABLES=H%2CPHIS%2CPL%2CPS%2CQI%2CQL%2CQV%2CT%2CU%2CV"

    delta = timedelta(1)
    n_parallel_downloads = 4

    print("#!/bin/bash")
    print("")

    if verbose:
        print("# {} {} {} {} {} {}".format(args.startdate, args.enddate, args.north, args.south, args.east, args.west))
        print("")

        d = start_date+delta
        version = str(max(1,int((d.year - 1971)/10)))+"00"
        print("# Sample URL="+URL.format(d, version=version, north=north, south=south, east=east, west=west))
        print("")

    for i in range(int((end_date-start_date)/delta)):
        d = start_date + delta*i

        version = str(max(1,int((d.year - 1971)/10)))+"00"

        if ((i % n_parallel_downloads) < (n_parallel_downloads - 1)):
            print(wget_exe+URL.format(d, version=version, north=north, south=south, east=east, west=west)+'" &')
        else:
            print(wget_exe+URL.format(d, version=version, north=north, south=south, east=east, west=west)+'"')

    print("")
    print("wait")
    print("")
    print("for i in MERRA2_*.nc; do ncpdq -a -lev $i icar_${i}; done")
    print("")
    print("ncrcat icar_MERRA2_*.nc merra_forcing_{0.year:04}_{0.month:02}_{0.day:02}.nc".format(start_date))
    print("")


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This script generates a shell script to download and format MERRA data for ICAR. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-s', "--startdate",
            dest="startdate",
            help="The Start Date - format YYYY-MM-DD",
            required=True,
            # default=datetime.datetime(2010, 1, 1),
            type=lambda d: datetime.datetime.strptime(d, '%Y-%m-%d'))
            # type=datetime.date.fromisoformat)
        parser.add_argument('-e', "--enddate",
            dest="enddate",
            help="The Start Date - format YYYY-MM-DD",
            required=True,
            # default=datetime.datetime(2010, 1, 31),
            type=lambda d: datetime.datetime.strptime(d, '%Y-%m-%d'))
            # type=datetime.date.fromisoformat)
        # parser.add_argument('filename',nargs="?", action='store', default="some_file_name")
        parser.add_argument('--south',  dest="south",   action='store', default="21.962")
        parser.add_argument('--west',   dest="west",    action='store', default="-131.133")
        parser.add_argument('--north',  dest="north",   action='store', default="55.009")
        parser.add_argument('--east',   dest="east",    action='store', default="-64.336")
        parser.add_argument('-v', '--version',action='version',
                version='get_merra 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        # global verbose
        verbose = args.verbose

        exit_code = main(args.startdate, args.enddate, args.north, args.south, args.east, args.west)

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
