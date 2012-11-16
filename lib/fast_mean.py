#!/usr/bin/env python
from scipy.weave import converters
from scipy import weave
import numpy as np
# surprisingly, this only cuts runtime in half compared to three python for loops around an np.mean
# could be re-written to updated the mean each window j-step instead of recomputing it entirely...
def smooth(data,windowsize):
    data2=data.copy()
    nz,nx,ny = data.shape

    code = """
            #line 10 "fast_mean.py"
            double tmp, n;
            n=windowsize*windowsize;
            for (int k=0;k<nz;k++){
                 for (int i=windowsize; i<nx-windowsize; i++) {
                     for (int j=windowsize; j<ny-windowsize; j++) {
                         tmp=0.0;
                         n=0;
                         for (int a=-1*windowsize;a<=windowsize;a++){
                             for (int b=-1*windowsize;b<=windowsize;b++){
                                tmp=tmp+double(data2(k,i+a,j+b));
                                n++;
                             }
                         }
                         data(k,i,j)=tmp/n;
                     }
                 }
                // tmp=double(k);
            }
            return_val = n;
            """
    # compiler keyword only needed on windows with MSVC installed
    err = weave.inline(code,
                       ['data', 'data2', 'windowsize', 'nx', 'ny','nz'],
                       type_converters=converters.blitz,compiler="gcc")
    return data
    

# Use the blitz converters already used for weave inlines
# and add to start of list so we don't get the catchall
# this was listed as a fix for a bug that affected linux systems, 
# but it doesn't seem to matter... so I'm not using it now
# from scipy.weave import blitz_spec
# class memmap_converter(blitz_spec.array_converter):
#     def init_info(self):
#         blitz_spec.array_converter.init_info(self)
#         self.matching_types = [np.memmap]
# eqrm_converters = [memmap_converter()] + converters.blitz

# This version is MUCH faster, ~50x faster than above code, 100x faster than pure python+numpy
# and output is identical to the above code. 
def fast_smooth(data,windowsize):
    data2=data.copy()
    nz,nx,ny = data.shape

    code = """
            #line 58 "fast_mean.py"
            double n,cursum,lastsum,nextsum;
            int j;
            n=(windowsize*2+1)*(windowsize*2+1);
            for (int k=0;k<nz;k++){
                 for (int i=windowsize; i<nx-windowsize; i++) {
                     cursum=0.0;
                     lastsum=0.0;
                     nextsum=0.0;
                     j=windowsize;
                     for (int a=-1*windowsize;a<=windowsize;a++){
                         for (int b=-1*windowsize;b<=windowsize;b++){
                            cursum=cursum+data2(k,i+a,j+b);
                         }
                     }
                     data(k,i,j)=cursum/n;
                     for (j=windowsize+1; j<ny-windowsize; j++) {
                         for (int a=-1*windowsize;a<=windowsize;a++){
                            cursum=cursum-data2(k,i+a,j-windowsize-1);
                            cursum=cursum+data2(k,i+a,j+windowsize);
                         }
                         //data(k,i,j)=cursum/n;
                         data(k,i,j)=float(cursum/n);
                     }
                 }
            }
            return_val = n;
            """
    # compiler keyword only needed on windows with MSVC installed
    err = weave.inline(code,
                       ['data', 'data2', 'windowsize', 'nx', 'ny','nz'],
                       type_converters=converters.blitz,compiler="gcc")
    return data
    