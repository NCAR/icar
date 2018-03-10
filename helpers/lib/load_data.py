'''
Created on Mar 18, 2011

@author: gutmann

    This module is for programs used to load data from various file formats
    Stripped down slightly to remove additional dependencies for the ICAR helpers directory

'''

import numpy as np
from datetime import datetime
import re
import os
import glob

def cols_date(filename,dtype='f',year=-1,month=-1,day=-1,hour=-1,minute=-1,second=-1):
    '''
    cols: reads a text file with an arbitrary number of header lines
    starts at first line that has nothing but numbers on it.
    returns a numpy array of data using numpy.loadtxt

    defaults to float data type, can be made to use other datatypes with dtype='d'

    EXAMPLE:
    import load_data
    data=load_data.cols('filename')
    '''
    from julday import mjul_day

    data=cols(filename,dtype=dtype)
    years=data[:,year]
    months=data[:,month]
    days=data[:,day]
    if hour>-1:
        hours = data[:,hour]
    else: hours=np.zeros(len(days))
    if minute>-1:
        minutes = data[:,minute]
    else: minutes=np.zeros(len(days))
    if second>-1:
        seconds = data[:,second]
    else: seconds=np.zeros(len(days))
    times=[mjul_day(years[i], months[i], days[i], hours[i], minutes[i], seconds[i]) for i in range(len(years))]
    #    times=[datetime(result[:,year],result[:,month],result[:,day],result[:,hour],result[:,minute])]
    times=np.array(times)
    data=data[:,max((year,month,day,hour,minute))+1:]
    return (times, data)


from bunch import Bunch
def cols(filename,dtype='d',delimiter=None,readheader=False):
    '''
    cols: reads a text file with an arbitrary number of header lines
    starts at first line that has nothing but numbers on it.
    returns a numpy array of data using numpy.loadtxt

    defaults to float data type, can be made to use other datatypes with dtype='d'

    EXAMPLE:
    import load_data
    data=load_data.cols('filename')
    '''
    from numpy import loadtxt
    from is_number import is_number

    headerdata=''
    f=open(filename, 'r')

    inheader=True
    headerlength=-1
    while inheader:
        headerlength+=1
        line=f.readline()
        curdata=line.split()
        inheader=False
        for test in curdata:
            inheader= (inheader or (not is_number(test)))

        if inheader & readheader:
            headerdata+=line

    f.close()

    if readheader:
        return Bunch(data=loadtxt(filename,skiprows=headerlength,dtype=dtype,delimiter=delimiter),header=headerdata)
    else:
        return loadtxt(filename,skiprows=headerlength,dtype=dtype,delimiter=delimiter)


def readfirst(data):
    line=data.next()
    for i in range(len(line)):
        line[i]=re.sub('(?P<num>[0-9])-','\g<num> ',line[i])
        line[i]=line[i].replace('/',' ').replace(':',' ').replace('+',' ').replace('kg',' ').replace('\'',' ').replace('"',' ')
    time=np.array(line[0].split(),'i')
    if time[0]<1800: time[0]+=2000
    time=datetime(time[0],time[1],time[2],time[3],time[4])
    outputtime=[datetime(2000,1,1) for i in range(100)]
    outputtime[0]=time

    outputdata=np.zeros((100,len(line)-1),'d')
    outputdata[0,:]=np.array(line[1:],'d')
    return (outputtime,outputdata)
