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
import re

def main (options_file, template_file):

    entered_restart_section=False

    with open(options_file,"r") as opt:
        with open(template_file,"w") as tmpl:
            for l in opt:

                if re.match("\s*&restart_info\s*",l):
                    entered_restart_section=True

                key = (l.split("=")[0]).strip().lower()
                if key == "restart_file":
                    if verbose: print("Writing restart_file line")
                    tmpl.write('    restart_file="__RESTART_FILE__"\n')
                elif key == "restart_date":
                    if verbose: print("Writing restart_date line")
                    tmpl.write('    restart_date= __RESTART_DATE__,\n')
                elif key == "restart_step":
                    if verbose: print("Removing restart_step line")
                    pass
                elif key == "restart":
                    if verbose: print("Writing restart=true line")
                    tmpl.write('    restart=true,\n')
                else:
                    tmpl.write(l)

            if not entered_restart_section:
                if verbose: print("Adding restart section")
                tmpl.write("\n")
                tmpl.write("&restart_info\n")
                tmpl.write('    restart_date= __RESTART_DATE__,\n')
                tmpl.write('    restart_file="__RESTART_FILE__"\n')
                tmpl.write("/\n")

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Simple script to make a setup template file based on an options file. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('option',   nargs="?", action='store', default="icar_options.nml")
        parser.add_argument('template', nargs="?", action='store', default="template.nml")
        parser.add_argument('-v', '--version',action='version',
                version='Make Template 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        global verbose
        verbose=args.verbose

        exit_code = main(args.option, args.template)
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
