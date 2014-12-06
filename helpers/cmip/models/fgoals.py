import numpy as np

import mygis
from bunch import Bunch

def vcoord(filename):
    """compute the vertical coordinate in space and time for a given file"""
    na=np.newaxis
    ptop = mygis.read_nc(filename,"ptop").data
    sigma = mygis.read_nc(filename,"lev").data[na,:,na,na]
    ps= mygis.read_nc(filename,"ps").data[:,na,:,:]
    p= ptop+sigma*(ps-ptop)
    return p

