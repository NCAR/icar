#!/usr/bin/env python
import os,shutil,glob,re
import sys,getopt

def install(install_dir,copyout):
    pyfiles=glob.glob("*.py")
    for curfile in pyfiles:
        if copyout:
            shutil.copy2(install_dir+curfile,curfile)
        else:
            shutil.copy2(curfile,install_dir)
    
    
        
if __name__ == '__main__':
    opts,args=getopt.getopt(sys.argv[1:],"r")
    print(opts)
    print(args)
    
    host=os.uname()[1].split(".")[0]
    if re.match("mirage",host):
        host="mirage"
    install_dirs=dict(lobuche="/Users/gutmann/Dropbox/bin/python/statdown/",
                      mirage="/glade/home/gutmann/bin/python/",
                      hydrolab="/home/gutmann/bin/python/",
                      nomad="/Users/gutmann/Dropbox/bin/python/statdown/",
                      pakaldi="/Users/gutmann/Dropbox/bin/python/statdown/")
    copyout=False
    if len(args)>0:
        if args[0][0]=="r":
            copyout=True
    if len(opts)>0:
        print(opts[0][0])
        if opts[0][0]=="-r":
            copyout=True
    print(copyout)
    install(install_dirs[host],copyout)
    os.chdir("lib")
    install(install_dirs[host]+"lib/",copyout)
    os.chdir("../")
                      