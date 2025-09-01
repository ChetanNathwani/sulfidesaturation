import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from oxygen_fugacity import OxygenFugacity
from fugacity import Fugacity
from sulfur_partition_coefficients import PartitionCoefficient
from Iacono_Marziano_COH import IaconoMarziano
from melt_composition import MeltComposition
from VC_COH import VolatileCalc
from newvariables import NewVariables
from degassingrun import COHS_degassing
from S_Fe import Sulfur_Iron
from SCSS_model import Sulfur_Saturation
from scipy.special import erf
import math as math
import PySulfSat as ss
import matplotlib as mpl
import periodictable as pt
import seaborn as sb

def S6_total_OandM(composition, T_C, del_FMQ):
    T_K = T_C + 273.15
    """o2: log10"""
    wtsio2 = composition["SiO2_Liq"]
    wttio2 = composition["TiO2_Liq"]
    wtal2o3 = composition["Al2O3_Liq"]
    wtfeo = composition["FeOt_Liq"]
    wtmno = composition["MnO_Liq"]
    wtmgo = composition["MgO_Liq"]
    wtcao = composition["CaO_Liq"]
    wtna2o = composition["Na2O_Liq"]
    wtk2o = composition["K2O_Liq"]
    wtp2o5 = composition["P2O5_Liq"]

    xtot = wtsio2/60.08+wttio2/79.9+wtal2o3/50.98+wtfeo/71.85+wtmgo/40.32+ wtcao/56.08+wtna2o/30.99+wtk2o/47.1+wtmno/70.94
    
    xna = (wtna2o/30.99)/xtot
    xmg = (wtmgo/40.32)/xtot
    xal = (wtal2o3/50.98)/xtot
    xsi = (wtsio2/60.08)/xtot
    xk = (wtk2o/47.1)/xtot
    xca = (wtcao/56.08)/xtot
    xti = (wttio2/79.9)/xtot
    xmn = (wtmno/70.94)/xtot
    xfet = (wtfeo/71.85)/xtot
    fe2_fetotal = 1/(1 + 10**(0.25*del_FMQ - 1.36 + 2.4*xca + 2*xna + 3.7*xk))
    xferrous = xfet * fe2_fetotal
    fo2 = del_FMQ - 25050/T_K+8.5
    
    c_sulfide = 8.77-23590/T_K+(1673/T_K)*(6.7*(xna+xk)+4.9*xmg+8.1*xca+8.9*(xfet+xmn)+5*xti+1.8*xal
                                                   -22.2*xti*(xfet+xmn)+7.2*((xfet+xmn)*xsi))-2.06*erf(-7.2*(xfet+xmn))
    c_sulfate = (-8.02) +(21100+44000*xna+18700*xmg+4300*xal+35600*xca+44200*xk+16500*xferrous+12600*xmn)/T_K
    lnk = (-55921)/T_K+25.07-0.6465*np.log(T_K) # SO3/S
    lnrs =(c_sulfate - lnk - c_sulfide) + 2 * np.log(10)*fo2
    rs =1-1/(1+np.exp(lnrs))

    return rs

def calc_fe2_fetotal(composition, T_C, del_FMQ):
    T_K = T_C + 273.15
    """o2: log10"""
    wtsio2 = composition["SiO2"]
    wttio2 = composition["TiO2"]
    wtal2o3 = composition["Al2O3"]
    wtfeo = composition["FeOT"]
    wtmno = composition["MnO"]
    wtmgo = composition["MgO"]
    wtcao = composition["CaO"]
    wtna2o = composition["Na2O"]
    wtk2o = composition["K2O"]
    wtp2o5 = composition["P2O5"]

    xtot = wtsio2/60.08+wttio2/79.9+wtal2o3/50.98+wtfeo/71.85+wtmgo/40.32+ wtcao/56.08+wtna2o/30.99+wtk2o/47.1+wtmno/70.94
    
    xna = (wtna2o/30.99)/xtot
    xmg = (wtmgo/40.32)/xtot
    xal = (wtal2o3/50.98)/xtot
    xsi = (wtsio2/60.08)/xtot
    xk = (wtk2o/47.1)/xtot
    xca = (wtcao/56.08)/xtot
    xti = (wttio2/79.9)/xtot
    xmn = (wtmno/70.94)/xtot
    xfet = (wtfeo/71.85)/xtot
    
    fe2_fetotal = 1/(1 + 10**(0.25*del_FMQ - 1.36 + 2.4*xca + 2*xna + 3.7*xk))
    
    return(fe2_fetotal)

def calc_logfo2(T_C, del_FMQ):
    T_K = T_C + 273.15
    return(del_FMQ - 25050/T_K+8.5)

# Calculate paritioning between MSS and silicate melt for Cu
def calc_DMSSSM_Cu(T, del_FMQ, FeOmelt):
    d = 1.18
    a = 0.28
    b = -0.05
    c = -0.66
    logd = d + a * (10000/T) + b*del_FMQ + c*np.log10(FeOmelt)
    D = 10**logd
    return D

def calc_H2O(SiO2):
    return a * SiO2**2 + b * SiO2 + c

def calc_DSLSM_Cu(T, del_FMQ, FeOmelt):
    d = 1.13
    a = 0.39
    b = 0.02
    c = -0.86
    logd = d + a * (10000/T) + b*del_FMQ + c*np.log10(FeOmelt)
    D = 10**logd
    return D

def calc_DMSSSM_Cu(T, del_FMQ, FeOmelt):
    d = 1.18
    a = 0.28
    b = -0.05
    c = -0.66
    logd = d + a * (10000/T) + b*del_FMQ + c*np.log10(FeOmelt)
    D = 10**logd
    return D

def calc_D(T, del_FMQ, FeOmelt, element, sulfur_phase):
    if element == 'Cu':
        if sulfur_phase == 'MSS':
            d = 1.18
            a = 0.28
            b = -0.05
            c = -0.66
        if sulfur_phase == 'SL':
            d = 1.13
            a = 0.39
            b = 0.02
            c = -0.86
    if element == 'Au':
        if sulfur_phase == 'MSS':
            d = 0.77
            a = 0.21
            b = -0.17
            c = -0.19
        if sulfur_phase == 'SL':
            d = 3.18
            a = 0.21
            b = 0.08
            c = -0.82
    if element == 'Ag': # Note, the R2 values for these are quite low
        if sulfur_phase == 'MSS':
            d = 0.73
            a = 0.16
            b = -0.06
            c = -0.34
        if sulfur_phase == 'SL':
            d = 2.01
            a = 0.16
            b = 0.00
            c = -0.16


    logd = d + a * (10000/T) + b*del_FMQ + c*np.log10(FeOmelt)
    D = 10**logd
    return D

def calc_melt_comp(Co, D, F, type):
    '''
        Calculate the melt composition based on initial concentration and melting model.

    Parameters
    ----------
    Co : float
        Initial concentration of the element in the source material.
    D : float
        Bulk distribution coefficient.
    F : float
        Melt fraction (in percent, e.g., 10 for 10%).
    type : str
        Type of melting model to use. Must be one of:
        - 'EQ' for equilibrium crystallisation.
        - 'FC' for fractional crystallization.

    Returns
    -------
    Cl : float
        Concentration of the element in the melt.

    '''
    if type == 'EQ':
        Cl = (Co/(D+F/100*(1-D)))
    elif type == 'FC':
        Cl = Co*(F/100)**(D-1)
    else:
        raise ValueError("Invalid type. Must be 'EQ' or 'FC'.")
    return Cl