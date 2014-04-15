#!/usr/bin/env python
import sys
from math import exp
import argparse,traceback,os

outputfile="sounding.txt"
z0=0
p_surf=1000.0
RH=0.95
dz=100
ztop=30000
U=10.0
V=0.0


"""gen_sounding.py [Surface T (K), [lapse_rate (K/km)]]
    
    Generate a sounding file for use by WRF ideal.exe
    Optionally supply the surface temperature in Kelvin
        default = 270 K
    if a surface T is supplied, you can optionally supply a lapse_rate
        Kelvin per Kilometer
        default = 5 K/km (stable but linear...)
"""

def t2vp(t):
    '''Saturated vapor pressure for a given temperature
    T in K'''
    # vp=6.112*exp(17.67*(t-273.15)/(t-29.65))
    if (t<273.15):
        a=21.8745584
        b=7.66
    else:
        a=17.2693882
        b=35.86
    
    vp = 610.78 * exp(a*(t-273.16)/(t-b)) #(Pa)
    
    return vp/100 #convert to mb / hPa


def rh2mr(t,p,rh):
    '''T(K), p(mb), rh(0-1)'''
    e=t2vp(t)*rh
    mr=0.62197*e/(p-e)
    return mr

def inverse_exner(theta,p):
    """convert potential temperature to temperature"""
    rcp=287.04/1004.0
    t= theta/(1000.0/p)**rcp
    return t

def exner(t,p):
    """compute potential temperature"""
    rcp=287.04/1004.0
    theta= t*(1000.0/p)**rcp
    return theta

def calc_p(z,z0,psurf):
    """Compute the pressure at a given level from some surface pressure"""
    slp=psurf/(1 - 2.25577E-5*z0)**5.25588  # if z0=0 then slp=psurf
    p = slp*(1 - 2.25577E-5*z)**5.25588
    return p

def main(t_surf=270.0,lapse_rate=5.0):
    """
    Generate a sounding for use by WRF ideal.exe
    
        t_surf = surface temperature (K) (not potential T, but they are "the same")
        lapse_rate = environmental lapse rate to use (K/km)
                NOTE:positive values means temperature decreases with height...
    """
    
    theta=t_surf
    t=inverse_exner(t_surf,p_surf)
    mr_surf=rh2mr(t,p_surf,RH)
    with open(outputfile,"wu") as f:
        f.write("{0} {1} {2}\n".format(p_surf,t_surf,mr_surf*1000))
        for z in range(z0,ztop,dz):
            p=calc_p(z,z0,p_surf)
            theta+=lapse_rate*(dz/1000.0)
            t=inverse_exner(theta,p)
            mr=rh2mr(t,p,RH)
            f.write("{0} {1} {2} {3} {4}\n".format(z,theta,mr*1000,U,V))
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Generate a moist sounding for WRF\'s ideal.exe')
        parser.add_argument('t_surf',nargs="?",action='store',help="Surface temperature [K] (not potential)",default=270)
        parser.add_argument('lapse_rate',nargs="?",action='store',help="Environmental lapse rate [K/km] (in potential T)",default=3)
        # parser.add_argument('p_surf',nargs="?",action='store',help="Surface pressure [mb]",default=1000)
        
        parser.add_argument('-v', '--version',action='version',
                version='gen_sounding.py 1.0')
        args = parser.parse_args()

        T_surface=float(args.t_surf)
        lapse_rate=float(args.lapse_rate)
        exit_code=main(T_surface, lapse_rate)
        
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
