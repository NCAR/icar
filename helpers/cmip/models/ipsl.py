import numpy as np

import mygis
from bunch import Bunch

def vcoord(filename):
    """compute the vertical coordinate in space and time for a given file"""
    na=np.newaxis
    ap= mygis.read_nc(filename,"ap").data[na,:,na,na]
    b = mygis.read_nc(filename,"b").data[na,:,na,na]
    ps= mygis.read_nc(filename,"ps").data[:,na,:,:]
    p= ap+b*ps
    return p

