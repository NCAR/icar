#!/usr/bin/env python
import os,shutil,glob,re

def install(install_dir):
    pyfiles=glob.glob("*.py")
    for curfile in pyfiles:
        shutil.copy2(curfile,install_dir)
    
    
        
if __name__ == '__main__':
    host=os.uname()[1].split(".")[0]
    if re.match("mirage",host):
        host="mirage"
    install_dirs=dict(lobuche="/Users/gutmann/Dropbox/bin/python/statdown/",
                      mirage="/glade/home/gutmann/bin/python/",
                      hydrolab="/home/gutmann/bin/python/",
                      nomad="/Users/gutmann/Dropbox/bin/python/statdown/",
                      pakaldi="/Users/gutmann/Dropbox/bin/python/statdown/")
    install(install_dirs[host])
    os.chdir("lib")
    install(install_dirs[host]+"lib/")
    os.chdir("../")
                      