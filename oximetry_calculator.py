from __future__ import division
import numpy as np
import math

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:15:11 2016

@author: Josh
"""

# Biophysical parameters from the description of umbilical cord relaxaton rates
R1_plas = 0.44 # 1/s
R1_ery_0 = 0.95 # 1/s
r_prime_dHB = 0.29 # 1/s
R2_plas = 1.60 # 1/s
R2_dia_plus_R2_oxy = 5.1 # 1/s
R2_deoxy_minus_R2_oxy = -3.1 # 1/s
w_dia_plus_w_oxy = -0.071 # ppm
w_deo_minus_w_oxy = -0.62 # ppm
tau = 2.27 # ms

def intermediate_constants(tau_180):
    mu = tau * (1 - 2*tau*math.tanh(tau_180/(2*tau))/tau_180)
    K0 = R2_plas
    K1 = R2_dia_plus_R2_oxy + mu*math.pow(w_dia_plus_w_oxy + w_deo_minus_w_oxy, 2)
    K2 = -mu*math.pow((w_dia_plus_w_oxy + w_deo_minus_w_oxy), 2)
    K3 = -R2_deoxy_minus_R2_oxy - 2*mu* w_deo_minus_w_oxy*(w_dia_plus_w_oxy + w_deo_minus_w_oxy)
    K4 = 2*mu*w_deo_minus_w_oxy * (w_dia_plus_w_oxy + w_deo_minus_w_oxy)
    K5 = mu*math.pow(w_deo_minus_w_oxy, 2)
    M0 = R1_plas
    M1 = R1_ery_0 - R1_plas + r_prime_dHB
    M2 = -r_prime_dHB
    
    return (K0, K1, K2, K3, K4, K5, M0, M1, M2)

def SO2_from_T1_Hct(T1, Hct):
    # equation 9
    R1 = 1/T1
    print R1

    sO2 = (R1 - R1_plas - Hct * (R1_ery_0 - R1_plas + r_prime_dHB) ) / (-r_prime_dHB*Hct)
    return sO2
    
def Hct_from_T1_sO2(T1, sO2):
    # equation 10
    R1 = 1/T1
    
    Hct = (R1 - R1_plas) / (R1_ery_0 - R1_plas + r_prime_dHB*(1-sO2))
    return Hct

def Hct_sO2_from_T1_T2(T1, T2, tau_180):
    # using equations 7 and 9
    R1 = 1/T1 # 1/s (to be set)
    R2 = 1/T2 # 1/s (to be set)

    (K0, K1, K2, K3, K4, K5, M0, M1, M2) = intermediate_constants(tau_180)
    
    A1 = K2 * math.pow(M2, 2) - K4*M1*M2 - K5*math.pow(M1, 2)
    B1 = K1 * math.pow(M2, 2) - K3*M1*M2 + K4*(R1*M2 - M0*M2) + K5*(math.pow(M1, 2) + 2*R1*M1 - 2*M0*M1)
    C1 = K0 * math.pow(M2, 2) + K3 * (R1*M2 - M0*M2) + K5 * (2*M0*M1 - 2*M1*R1 - math.pow(R1-M0, 2)) - R2 * math.pow(M2, 2)
    D1 = K5 * math.pow(R1-M0, 2)
    
    # equation 7
    polynomial_coefficients = [A1, B1, C1, D1]    
    roots = np.roots(polynomial_coefficients)
    Hct_vals = np.real(roots)
    
    SO2_vals = []
    for hct_root in Hct_vals:
        SO2_vals.append(SO2_from_T1_Hct(T1, hct_root))
    
    return Hct_vals, SO2_vals
    
def alt_Hct_sO2_from_T1_T2(T1, T2, tau_180):
    # using equations 8 and 10
    R1 = 1/T1 # 1/s (to be set)
    R2 = 1/T2 # 1/s (to be set)
    
    (K0, K1, K2, K3, K4, K5, M0, M1, M2) = intermediate_constants(tau_180)
    
    A2 = K5 *(R1*M2 - M0*M2)
    B2 = K5 * (R1*M1 - M0*M1 - math.pow(M0-R1, 2)) + K3 * (R1*M2 - M0*M2) + K0*math.pow(M2, 2) - R2*math.pow(M2, 2)
    C2 = K4 * math.pow(M0-R1, 2) + K3*(R1*M1 - M0*M1) + K1 * (R1*M2 - M0*M2) + 2*K0*M1*M2 - 2*R2*M1*M2
    D2 = K2 * math.pow(M0-R1, 2) + K1 * (R1*M1 - M0*M1) + K0*math.pow(M1, 2) - R2*math.pow(M1, 2)
    
    # equation 7
    polynomial_coefficients = [A2, B2, C2, D2]  
    roots = np.roots(polynomial_coefficients)
    SO2_vals = np.real(roots)
    
    Hct_vals = []
    for so2_root in SO2_vals:
        Hct_vals.append(Hct_from_T1_sO2(T1, so2_root))
    return Hct_vals, SO2_vals
    
def sO2_from_T2_Hct(T2, Hct, tau_180):
    # equation 11
    R2 = 1/T2
    (K0, K1, K2, K3, K4, K5, M0, M1, M2) = intermediate_constants(tau_180)    
    
    polynomial_coeff = [K5 *Hct - K5*math.pow(Hct, 2), K3*Hct + K4 * math.pow(Hct, 2), K0 + K1 * Hct + K2 * math.pow(Hct, 2) - R2]
    sO2_roots = np.real(np.roots(polynomial_coeff))
    return sO2_roots

def test_with_T1_and_T2():
     T1 = 1.3
     T2 = 0.15
     tau_180 = 20
     
     Hct_vals, SO2_vals = Hct_sO2_from_T1_T2(T1, T2, tau_180)
     print Hct_vals, SO2_vals
    
test_with_T1_and_T2()