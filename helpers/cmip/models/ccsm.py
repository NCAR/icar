import numpy as np

import mygis
from bunch import Bunch

def vcoord(filename):
    """compute the vertical coordinate in space and time for a given file"""
    na=np.newaxis
    a = mygis.read_nc(filename,"a").data[na,:,na,na]
    b = mygis.read_nc(filename,"b").data[na,:,na,na]
    p0= mygis.read_nc(filename,"p0").data
    ps= mygis.read_nc(filename,"ps").data[:,na,:,:]
    p= a*p0+b*ps
    return p

