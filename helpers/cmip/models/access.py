import numpy as np

import mygis
from bunch import Bunch

def vcoord(filename):
    """compute the vertical coordinate in space and time for a given file"""
    na=np.newaxis
    a = mygis.read_nc(filename,"lev").data[na,:,na,na]
    b = mygis.read_nc(filename,"b").data[na,:,na,na]
    orog= mygis.read_nc(filename,"orog").data
    z= a+b*orog
    return z

