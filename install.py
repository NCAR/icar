#!/usr/bin/env python
import os,shutil,glob,re
import sys,getopt,subprocess

def install(install_dir,copyout,diff):
    pyfiles=glob.glob("*.py")
    for curfile in pyfiles:
        if diff:
            print("Diff for: "+curfile)
            subprocess.call("diff "+curfile+" "+install_dir+curfile,shell=True)
        else:
            if copyout:
                shutil.copy2(install_dir+curfile,curfile)
            else:
                shutil.copy2(curfile,install_dir)
    
        
if __name__ == '__main__':
    opts,args=getopt.getopt(sys.argv[1:],"rd")
    
    host=os.uname()[1].split(".")[0]
    if re.match("mirage",host):
        host="mirage"
    install_dirs=dict(lobuche="/Users/gutmann/Dropbox/bin/python/statdown/",
                      mirage="/glade/home/gutmann/bin/python/",
                      hydrolab="/home/gutmann/bin/python/",
                      nomad="/Users/gutmann/Dropbox/bin/python/statdown/",
                      pakaldi="/Users/gutmann/Dropbox/bin/python/statdown/")
    # print(host,install_dirs[host])
    copyout=False
    diff=False
    if "r" in args:
        copyout=True
    if "d" in args:
        diff=True
    for opt in opts:
        if opt[0]=="-r":
            copyout=True
        if opt[0]=="-d":
            diff=True
    if copyout:
        print("Copying files out from "+install_dirs[host])
    if diff:
        print("Differencing files from "+install_dirs[host])
    install(install_dirs[host],copyout,diff)
    os.chdir("lib")
    install(install_dirs[host]+"lib/",copyout,diff)
    os.chdir("../")
                      