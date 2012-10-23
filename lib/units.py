#!/usr/bin/env python
# encoding: utf-8
"""
units.py

converts between various formulations of atmospheric humidity
    specific humidity (sh)
    relative humidity (rh)
    dewpoint (dp)

    xy2a = convert x,y point (or u/v vector) to the corresponding angle (0,0) to (x,y)

    for atmospheric (LOP/SWIM) purposes can compute lifting condensation level (LCL) 
    and the temperature at the LCL

    could also add
    vapor pressure (vp)

Created by Ethan Gutmann on 2011-08-15.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""
import numpy as np

def rtod(angle):
    return angle*360.0/(2*np.pi)
def dtor(angle):
    return angle/360.0*(2*np.pi)
    
def xy2a(x,y):
    '''
    convert an x,y coordinate to the angle to that coordinate (0-360)
    '''
    
    angle=rtod(np.arctan(np.cast['f'](x)/y))
    if np.array(x).size==1:
        if x>0:
            if y>0:
                pass
            else:
                angle=180+angle
        else:
            if y>0:
                angle=360+angle
            else:
                angle=180+angle
                
    else:
        for i in range(len(x)):
            if x[i]>0:
                if y[i]>0:
                    pass
                else:
                    angle[i]=180+angle[i]
            else:
                if y[i]>0:
                    angle[i]=360+angle[i]
                else:
                    angle[i]=180+angle[i]
            
            
    return angle


def t2vp(t):
    '''Saturated vapor pressure for a given temperature
    T in degrees C or K'''
    if np.array(t).max()>150:
        tIsInK=True
        t-=273.15
    else:
        tIsInK=False
    vp=6.112*np.exp(17.67*t/(t+243.5)) #*10**((7.5*t)/(237.7+t))
    if tIsInK:
        t+=273.15
    return vp

def sh2rh(t,p,sh):
    '''T(deg.C or K), p(mb), Specific Humidity(kg/kg)'''
    mr=sh/(1-sh)
    return mr2rh(t,p,mr)

def rh2sh(t,p,rh):
    '''T(deg.C or K), p(mb), Relative Humidity(0-1 or 0-100)'''
    mr=rh2mr(t,p,rh)
    return mr/(1+mr)

def sh2dp(t,p,sh):
    '''T(deg.C or K), p(mb), SH(kg/kg)'''
    return rh2dp(t,sh2rh(t,p,sh))
    
def mr2rh(t,p,mr):
    '''T(deg.C or K), p(mb), MR(kg/kg)'''
    # e=mr*p/(621.97+mr)
    e=mr*p/(0.62197+mr)
    rh=e/t2vp(t)
    return rh

def rh2mr(t,p,rh):
    '''T(deg.C or K), p(mb), rh(0-1 or 0-100 assumed if max is over 1.5)'''
    if np.array(rh).max()>1.5:
        rh/=100.0
        rhIsPercent
    e=t2vp(t)*rh
    mr=0.62197*e/(p-e)
    if rhIsPercent:
        rh*=100.0
    return mr

def dp2mr(t,p,dp):
    '''T(deg.C or K), p(mb), dp(deg.C or K)'''
    return rh2mr(t,p,dp2rh(t,dp))

def mr2dp(t,p,mr):
    '''T(deg.C or K), p(mb), mr(kg/kg)'''
    if np.array(t).max()>150:
        tIsInK=True
        t-=273.15
    else:
        tIsInK=False
    dp=rh2dp(t,mr2rh(t,p,mr))
    if tIsInK:
        t+=273.15
        dp+=273.15
    return dp

def rh2dp(t,rh):
    '''T in deg.C or K RH:0-1 or 0-100 if RH>1'''
    if np.array(rh).max()>1:rh/=100.0
    if np.array(t).max()>150:
        tIsInK=True
        t-=273.15
    else:
        tIsInK=False
    e=t2vp(t) * rh
    dp=243.5*np.log(e/6.112)/(17.67-np.log(e/6.112))
    if tIsInK: 
        t+=273.15
        dp+=273.15
    return dp

def dp2rh(t,dp):
    '''T and DP in deg.C or K (T and DP can be in different units)'''
    e=t2vp(dp)
    esat=t2vp(t)
    return np.max((np.min((0,e/esat)),1))

def find_lcl(p,t,dp):
    '''Find the lifting condensation level (m) p=[Pa,hPa],t[C,K],dp[C,K]
    
    Units are permitted to vary, 
        p can be hPa, or Pa
        t and dp can be C or K
            but both t and dp must be the same (C or K)
    
    Main code from Greg Thompson's FORTRAN77 (meteo_funcs.f) based on Bolton (1980)?
    '''
    if np.array(t).max()<150:
        tIsInC=True
        t+=273.15
        dp+=273.15
    else:
        tIsInC=False
    if np.array(p).max()<1500:
        pIsInHPa=True
        p*=100.0
    else:
        pIsInHPa=False
    
    rcp=287.04/1004.0
    cpr=1.0/rcp
    
    tlcl=t_lcl(t,dp)
    theta=t*(100000.0/p)**rcp
    plcl = 100000.0 * (tlcl/theta)**cpr
    lcl=(1-((plcl/101325.0)**(1.0/5.25588)))/2.25577E-5
    
    if pIsInHPa:
        p/=100
    if tIsInC:
        t-=273.15
        dp-=273.15
    return lcl
    

def t_lcl(t,dp):
    '''Convert T and dew point to T at lifting condensation level
    
    Units can be in C or K, if T<150 assumes C and returns C
    if T>=150 assumes K and returns K
    Code from Greg Thompson FORTRAN77 (meteo_funcs.f) based on Bolton (1980) eqn #15
    claims accuracy of 0.1K for -35<T<35 C
    '''
    tIsInC=False
    if np.array(t).max()<150:
        tIsInC=True
        t+=273.15
        dp+=273.15
    t_lcl=(1.0/(1.0/((dp-56.0))+(np.log(t/dp)/800.))) + 56.0
    if tIsInC:
        t-=273.15
        dp-=273.15
        t_lcl-=273.15
    return t_lcl
