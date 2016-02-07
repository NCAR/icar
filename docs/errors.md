##Common errors

1) **Segmentation fault**:
    Most likely due to your shell's stacksize limit (particularly with ifort). 
    
To *fix* it try the following: 

    in bash:
        ulimit -s unlimited
    in csh: 
        unlimit stacksize
        
**NB**: Some systems have a hard limit of 64MB (ulimit -s 65532), this *may* be enough depending on domain size. 
Reference: https://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors
