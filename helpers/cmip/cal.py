import numpy as np

import datetime
import mygis
from bunch import Bunch

days_per_month=[31,28,31,30,31,30,31,31,30,31,30,31]
start_day_per_month=[0]
start_day_per_month.extend(np.cumsum(days_per_month))

def std_date(model_time,roundseconds=True,y0=1850,m0=1,filename=None):
    """docstring for std_date"""
    if y0==None:
        y0=find_start_year(filename)
    if m0==None:
        m0=find_start_month(filename)

    curdate=datetime.datetime(int(y0),int(m0),1,0,0)+datetime.timedelta(model_time)
    year=curdate.year
    month=curdate.month
    day=curdate.day
    hour=curdate.hour
    minute=curdate.minute
    second=curdate.second
    if roundseconds:
        if np.round(second)>np.floor(second):
            if np.round(second)==60:
                second=00
                year,month,day,hour,minute=add_to_minute(year,month,day,hour,minute)
            else:
                second=np.round(second)
        else:
            second=np.round(second)
    
    date=Bunch(year=year,month=month,day=day,hour=hour,minute=minute,second=second)
    date[0]=year
    date[1]=month
    date[2]=day
    date[3]=hour
    date[4]=minute
    date[5]=second
    return date


def noleap_month_from_doy(doy):
    """docstring for noleap_month_from_doy"""
    for m,s in enumerate(start_day_per_month):
        if s>np.floor(doy):
            return m # m is zero based, otherwise would return m-1
def noleap_day_from_month_doy(doy):
    """docstring for noleap_month_from_doy"""
    for m,s in enumerate(start_day_per_month):
        if s>np.floor(doy):
            return m # m is zero based, otherwise would return m-1

def add_to_month(year,month):
    month=month+1
    if month==13:
        month=1
        year=year+1
    return year,month

def add_to_day(year,month,day):
    day=day+1
    if day==days_per_month[month]+1:
        day=1
        year,month=add_to_month(year,month)
    return year,month,day

def add_to_hour(year,month,day,hour):
    hour=hour+1
    if hour==24:
        hour=0
        year,month,day=add_to_day(year,month,day)
    return year,month,day,hour

def add_to_minute(year,month,day,hour,minute):
    minute=minute+1
    if minute==60:
        minute=0
        add_to_hour(year,month,day,hour)
    return year,month,day,hour,minute

def find_start_year(filename):
    """Try to determine the starting year for the calendar in a given file"""
    units_string=mygis.read_attr(filename,"units",varname="time")
    # e.g. "days since 1960-01-01 00:00:00"
    date=units_string.split()[2]
    year=int(date[:4])
    return year

def find_start_month(filename):
    """Try to determine the starting month for the calendar in a given file"""
    units_string=mygis.read_attr(filename,"units",varname="time")
    # e.g. "days since 1960-01-01 00:00:00"
    date=units_string.split()[2]
    month=int(date.split("-")[1])
    return month

def noleap_date(model_time,roundseconds=True,y0=1850,m0=1,filename=None):
    """docstring for noleap_date"""
    if y0==None:
        y0=find_start_year(filename)
    if m0==None:
        m0=find_start_month(filename)
    if m0!=1:
        model_time+=start_day_per_month[m0-1]
    
    year=np.floor(model_time/365.0)+y0
    doy=(model_time%365)
    month=noleap_month_from_doy(doy)
    day=np.floor(doy-start_day_per_month[month-1])+1
    dayfrac=doy%1.0
    
    hour=np.floor(dayfrac*24)
    minute=np.floor(dayfrac*1440 - hour*60)
    second=dayfrac*86400-((hour*60)+minute)*60
    if roundseconds:
        if np.round(second)>np.floor(second):
            if np.round(second)==60:
                second=00
                year,month,day,hour,minute=add_to_minute(year,month,day,hour,minute)
            else:
                second=np.round(second)
        else:
            second=np.round(second)
    
    date=Bunch(year=year,month=month,day=day,hour=hour,minute=minute,second=second)
    date[0]=year
    date[1]=month
    date[2]=day
    date[3]=hour
    date[4]=minute
    date[5]=second
    return date
    