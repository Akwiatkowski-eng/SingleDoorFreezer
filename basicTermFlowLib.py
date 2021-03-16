# This library is created to contain all basic functions which are usefull in 
# thermodynamic calculations
# 
# creation date 24.10.2020
# last mod 24.10.2020

import math as m
import CoolProp as CP


# -------------Reynolds number----------------
#this function calculates reynolds number
def Re(u,d,v): 
    """
    this function calculates reynolds number

    Parameters
    ----------
    u : characteristic velocity, m/s
    
    d : characteristic dimension, m
    
    v : kinematic viscosity, m2/s 
        

    Returns
    -------
    Re : Reynolds number.

    """
    
    Re = u*d/v
    return Re

#-------------pressure drop coef----------
#this function determines pressuredrop 
def pdCoeff(Re):
    """
    this function determines pressuredrop 

    Parameters
    ----------
    Re : Reynolds number from 0 to 10e5.

    Returns
    -------
    lambda1 : pressure drop coefficient.

    """
    if Re<=2400.:
        lambda1 = 64./Re
    if Re<1000000. and  Re>2400.:
        lambda1 = 0.316/Re**(0.25)
    #else: 
        #lambda1 = print('I cannot calculate for Re>10e5')
    return lambda1


    




# -----------Pressure drop in straight line -----------
def pressDrop(lambda2, u, d, ro):
    """
    this function calculate pressure drop related to tube length

    Parameters
    ----------
    lambda1 : pressure drop coefficent
        
    u : flow velocity, m/s
        
    d : hydraulic diameter, m

    ro : fluid density kg/m3
        .

    Returns
    -------
    pd : pressure drop per unit length, Pa/m

    """
    pd = lambda2 * (ro/2.) * pow(u,2) /d
    return pd


#-----------Velocity from mass flow ---------------
#this function calculates velocity from mass flow in circular channel
def velFromMas(m,d,ro):
    """
    this function calculates velocity from mass flow in circular channel

    Parameters
    ----------
    m : mass flow, kg/s
    
    d : tube diameter, m
    
    ro : fluid density, kg/m3

    Returns
    -------
    u : velocity m/s.

    """
    
    u = 4*m/(3.14*d**2*ro)
    return u

#----------Conduction heat transfer coeff--------------
def condHT(case,l1,l2,l3,lambda1,lambda2,lambda3,):
    """
    This function calculates conduction heat transfer coeff for 3 walls 

    Parameters
    ----------
    case : type 'wall' for flat wall conduction or 'cylinder' for cylindrical conductio.
    
    l1 : thickness of first wall 
        .
    l2 : thickness of second wall 
        . 
    l3 : thickness of third wall 
        .
    lambda1 : heat conduction coeff for third wall material
        
    lambda2 : heat conduction coeff for second wall material
        
    lambda3 : heat conduction coeff for third wall material
        

    Returns
    -------
    k : overall conduction heat transfer.

    """
    if case == 'wall':
        k = 1./(l1/lambda1 + l2/lambda2 + l3/lambda3)
    if case == 'cylinder':
        k = 1 # is going to be defined later
    return k

#-----------Average Temp-----------------------
def tempAvg(T1,T2):
    """
    this function calculates average temp out of two

    Parameters
    ----------
    T1 : first temperature
    .
    T2 : second temperature
    .

    Returns
    -------
    tavg : average temperature
    .

    """
    tavg = (T1 + T2)/2.
    return tavg

#---------Logarythmic mean temperature difference
def LMTD(dT1,dT2):
    """
    This function calculates logarythmic mean temp difference

    Parameters
    ----------
    dT1 : temp difference at inlet
    .
    dT2 : temp difference ot outlet.

    Returns
    -------
    LMTD : logarythmic mean temperatue difference.

    """
    
    LMTD = (dT1 - dT2)/m.log(dT1/dT2)
    return LMTD



#Re1 = Re(10, 1, 0.0005)    

#print(condHT('wall',0.05,0.03,0.1,200,100,20))