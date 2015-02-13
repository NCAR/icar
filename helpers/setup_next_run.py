#!/usr/bin/env python

import glob,os,re,sys, fnmatch
from math import floor

import mygis

def find_last_output(options_file):
    """docstring for find_last_output"""
    with open(options_file,"rU") as f:
        for l in f:
            ltest=l.split("!")[0] # only look at the line before any comment characters
            if fnmatch.fnmatch(ltest.strip(),"output_file*=*"):
                outputfile=l.split("=")[1:]
                if len(outputfile)>1:
                    outputfile="=".join(outputfile)
                else:
                    outputfile=outputfile[0]
                outputfile=outputfile.replace('"',"").replace("'","").replace(',',"").replace('\n',"")
    
    print("Found outputfile prefix: "+outputfile)
    filelist=glob.glob(outputfile+"*00.nc")
    filelist.sort()

    return filelist[-1], filelist[-2]

def load_last_date(filename):
    """docstring for load_last_date"""
    time_data=mygis.read_nc(filename,"time").data
    last_hour=len(time_data)
    
    year=filename.split("_")[-4]
    month=filename.split("_")[-3]
    day=filename.split("_")[-2]
    hour=last_hour-1 # look back one time step incase it hadn't finished writing the last timestep
    
    return "{0}, {1}, {2}, {3}, 0, 0".format(year,month,day,hour),hour
    
    
def main(options_file):
    """docstring for main"""
    restart_file,backup_file=find_last_output(options_file)
    restart_date,hour=load_last_date(restart_file)
    if hour<2:
        restart_file=backup_file
        restart_date,hour=load_last_date(restart_file)
    
    print("Using restart file: "+restart_file)
    print("      restart date: "+restart_date)
    
    # escape the filename for sed
    restart_file=restart_file.replace("/","\/")
    os.system("sed -e 's/__RESTART_FILE__/"+restart_file   +"/g;"+
                      "s/__RESTART_DATE__/"+str(restart_date)   +"/g;"+
                      "' template.nml >"+options_file)
    
        
    
    

if __name__ == '__main__':
    if len(sys.argv)>1:
        options_file=sys.argv[1]
    else:
        options_file=glob.glob("icar_*opt.nml")[0]
        
    main(options_file)