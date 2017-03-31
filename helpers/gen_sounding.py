#!/usr/bin/env python
import sys
from math import exp
import numpy as np
import argparse,traceback,os
import units

global U, RH
outputfile="sounding.txt"
following_moist_adiabat=False
z0=0
p_surf=1000.0
RH=0.95
dz=100
ztop=25000
U=10.0
V=0.0

# useful ref:
# D Bolton, 1980: The Computation of Equivalent Potential Temperature. Mon. Wea. Rev., Vol. 108, pp.1046-1053.

"""gen_sounding.py [Surface T (K), [lapse_rate (K/km)]]
    
    Generate a sounding file for use by WRF ideal.exe
    Optionally supply the surface temperature in Kelvin
        default = 270 K
    if a surface T is supplied, you can optionally supply a lapse_rate
        Kelvin per Kilometer
        default = 5 K/km (stable but linear...)
"""

# From G. Thompson's meteo_funcs.f
#       function theta_e(pres_Pa,temp_K,w_non,tlcl_K)
def theta_e(pres,temp,mr,tlcl):
    """ The following code was based on Bolton (1980) eqn #43
    and claims to have 0.3 K maximum error within -35 < T < 35 C
    pres  = Pressure in Pascals
    temp  = Temperature in Kelvin
    mr    = mixing ratio (non-dimensional = kg/kg)
    tlcl  = Temperature at Lifting Condensation Level (K)"""
    
    mr = max(mr, 1.e-8)
    power=(0.2854*(1.0 - (0.28*mr) ))
    xx = temp * (100000.0/pres)**power
    
    p1 = (3.376/tlcl) - 0.00254
    p2 = (mr*1000.0) * (1.0 + 0.81*mr)
    
    return xx * np.exp(p1*p2)
    
# c
# c+---+-----------------------------------------------------------------+
# c
#       function t_lcl(temp_K,tdew_K)
# c..
# c..         The following code was based on Bolton (1980) eqn #15
# c..         and claims to have 0.1 K maximum error within -35 < T < 35 C
# c..            temp_K  = Temperature in Kelvin
# c..            tdew_K  = Dewpoint T at Lifting Condensation Level (K)
# c..
#       tt = temp_K
#       tttd= tdew_K
#       denom= ( 1.0/(tttd-56.0) ) + (log(tt/tttd)/800.)
#       t_lcl = ( 1.0 / denom ) + 56.0
#       return
#       end
# c
# c+---+-----------------------------------------------------------------+
# c
#       function theta_wetb(thetae_K)
def theta_wetb(thetae_K):
    """
        Eqn below was gotten from polynomial fit to data in
        Smithsonian Meteorological Tables showing Theta-e
        and Theta-w
    """
    c=np.array([-1.00922292e-10, -1.47945344e-8, -1.7303757e-6,
                -0.00012709,      1.15849867e-6, -3.518296861e-9,
                3.5741522e-12])
    d=np.array([0.00000000,   -3.5223513e-10, -5.7250807e-8,
                -5.83975422e-6, 4.72445163e-8, -1.13402845e-10,
                8.729580402e-14])
    
    x=min(475.0,thetae_K)
    if( x < 335.5 ) :
        answer = c[0]+x*(c[1]+x*(c[2]+x*(c[3]+x*(c[4]+x*(c[5]+
                    x*c[6] )))))
    else:
        answer = d[0]+x*(d[1]+x*(d[2]+x*(d[3]+x*(d[4]+x*(d[5]+
                    x*d[6] )))))
    
    return answer + 273.15
# c
# c+---+-----------------------------------------------------------------+
# c
#       function compT_fr_The(thelcl_K,pres_Pa)
def compT_fr_The(thelcl_K,pres_Pa):
    """
    Compute temperature at a given pressure from theta_e at LCL
    
    pres_Pa = Pressure in Pascals
    thelcl  = Theta-e at LCL (units in Kelvin)
    
    Temperature (K) is returned given Theta-e at LCL
    and a pressure.  This describes a moist-adiabat.
    This temperature is the parcel temp at level Pres
    along moist adiabat described by theta-e.
    """
    guess= (thelcl_K - 0.5 * ( max(thelcl_K-270., 0.))**1.05) * (pres_Pa/100000.0)**.2
    epsilon=0.01
    for i in range(100):
        w1 = t2mr(pres_Pa/100.,guess)
        w2 = t2mr(pres_Pa/100.,guess+1.)
        tenu = theta_e(pres_Pa,guess,w1,guess)
        tenup = theta_e(pres_Pa,guess+1,w2,guess+1.)
        cor = (thelcl_K - tenu) / (tenup - tenu)
        guess = guess + cor
        if (abs(cor)<epsilon):
            return guess

    print('  convergence not reached ')
    thwlcl_K=theta_wetb(thelcl_K)
    return thwlcl_K*((pres_Pa/100000.0)**0.286)

# end G. Thompson routines
# c+---+-----------------------------------------------------------------+
# c

def t2mr(pres,t):
    """saturated mixing ratio at a given T and P
    T in K
    P in hPa(mb)"""
    vpsat=t2vp(t)
    return 0.622*vpsat / (pres-vpsat)

def t2vp(t):
    '''Saturated vapor pressure for a given temperature
    T in K'''
    # vp=6.112*exp(17.67*(t-273.15)/(t-29.65))
    # saturated vapor pressure with respect to ice
    if (t<273.15):
        a=21.8745584
        b=7.66
    else:
    # saturated vapor pressure with respect to liquid water
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
    
        t_surf = surface temperature (K) (not potential T, but they are "the same" p=p0)
        lapse_rate = environmental lapse rate to use (K/km)
                NOTE:positive values means temperature decreases with height...
    """
    
    theta=t_surf
    tlcl=t_surf
    t=inverse_exner(t_surf,p_surf)
    mr_surf=rh2mr(t,p_surf,RH)
    mr=mr_surf
    if following_moist_adiabat:
        print("Following the moist adiabat")
    else:
        print("Following a constant lapse rate {} dT/dz (K/km)".format(lapse_rate))
    with open(outputfile,"wu") as f:
        f.write("{0} {1} {2}\n".format(p_surf,t_surf,mr_surf*1000))
        # th_e_lcl = theta_e(p_surf*100,t_surf,mr_surf,t_surf)
        lastz=z0
        lastp=p_surf
        for z in range(z0,ztop,dz):
            # p=calc_p(z,z0,p_surf)
            p=units.zt2p(z-lastz,p0=lastp,t0=t,dtdz= (lapse_rate-9.8)/1000.0)
            lastz=z
            lastp=p
            if following_moist_adiabat:
                t=compT_fr_The(th_e_lcl,p*100)
                theta=exner(t,p)
                mr=rh2mr(t,p,RH)
            else:
                theta+=lapse_rate*(dz/1000.0)
                t=inverse_exner(theta,p)
                mr=rh2mr(t,p,RH)
            
            f.write("{0} {1} {2} {3} {4}\n".format(z,theta,mr*1000,U,V))
    

if __name__ == '__main__':
    # global U,RH
    try:
        parser= argparse.ArgumentParser(description='Generate a moist sounding for WRF\'s ideal.exe')
        parser.add_argument('t_surf',nargs="?",action='store',help="Surface temperature [K] (not potential) <270>",default=270)
        parser.add_argument('lapse_rate',nargs="?",action='store',help="Environmental lapse rate [K/km] (in potential T) <3> ",default=3)
        parser.add_argument('U',nargs="?",action='store',help="Wind Speed [m/s] <10>",default=10)
        parser.add_argument('RH',nargs="?",action='store',help="Relative Humidity [0-1] <0.95>",default=0.95)
        # parser.add_argument('p_surf',nargs="?",action='store',help="Surface pressure [mb]",default=1000)
        
        parser.add_argument('-v', '--version',action='version',
                version='gen_sounding.py 1.0')
        args = parser.parse_args()
        
        U=float(args.U)
        RH=float(args.RH)

        T_surface=float(args.t_surf)
        lapse_rate=float(args.lapse_rate)
        exit_code=main(T_surface, lapse_rate)
        
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
