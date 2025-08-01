import numpy as np
import pandas as pd
import math
import cmath
import time
import re
import os
from scipy.fft import fft, fftfreq
from scipy import constants as cnst
from matplotlib import pyplot as plt
from lmfit import Model
# is this needed, if importing through circuit_reader?

# A file that keeps all calculation methods

# definition of constants
phi0 = cnst.h/(2*cnst.e)

# probably do a use_calc_val data format transformation

# general circuit/jj functions
def calc_lc_rfreq(l, c):
    # frequency, not omega!
    return 1/(2*np.pi*np.sqrt(l*c))

def calc_c_from_rfreq(rfreq, l):
    return 1/(4*np.pi**2*l*rfreq**2)

def calc_lc_q(r, l, c, mode="series"):
    if mode=="series" or mode=="s":
        qfactor = 1/r*np.sqrt(l/c)
    elif mode=="parallel" or mode=="p":
        qfactor = r*np.sqrt(c/l)
    else:
        qfactor = 0
    return qfactor

def calc_k(l1, l2, m):
    return m/np.sqrt(l1*l2)

def calc_lj(ic, phi=0):
    return phi0/(2*np.pi*ic)/np.cos(phi)

def calc_lj_ibias(ic, ibias):  # incomplete
    ibias = ibias % (2*ic)
    phi = np.arcsin(ibias/ic)
    return phi0/(2*np.pi*ic)/np.cos(phi)

def calc_fj(ic):
    return ic/(4*np.pi*cnst.e)
    
def calc_fj_jjparams(ic, cpic, phi=0):
    return 1/(2*np.pi*np.sqrt(np.abs(calc_lj(ic, phi))*cpic*ic))

# for sinusoidal fns.
def amplitude(data):
    return (np.max(data) - np.min(data))/2

def offset(data):
    return np.max(data) - (np.max(data) - np.min(data))/2
    
def dB2power(x):  # for vna
    return pow(10, x / 10.)

def power2dB(x):  # for vna
    return 10 * np.log10(x)

# rifkin circuit specific functions
def calc_leff(lt, omega, r, l, cur_l, k=1, ls=0, stabilize_tank=True):
    m2 = k**2*lt*l
    if r == 0 and stabilize_tank:
        leff = lt - m2/(l + cur_l + ls)
    elif not stabilize_tank:
        leff = lt + ls - m2/(l + cur_l)
    else:
        z_0r = omega**2*cur_l**2*r/(r**2 + omega**2*cur_l**2)
        z_0i = omega*cur_l*r**2/(r**2 + omega**2*cur_l**2) + omega*ls
        leff = lt - omega*m2*(z_0i+omega*l)/(z_0r**2 + (z_0i+omega*l)**2)
    # leff = lt - (r**2*m2*(l+cur_l) + omega**2*m2*cur_l**2*l)/(r**2*(cur_l+l)**2+omega**2*cur_l**2*l**2)
    return leff

def calc_reff(rt, lt, omega, r, l, cur_l, k=1, ls=0):
    m2 = k**2*lt*l
    if r == 0:
        try: reff = rt * np.ones(len(cur_l))
        except TypeError: reff = rt
    else:
        z_0r = omega**2*cur_l**2*r/(r**2 + omega**2*cur_l**2)
        z_0i = omega*cur_l*r**2/(r**2 + omega**2*cur_l**2) + ls
        reff = rt + omega**2*m2*z_0r/(z_0r**2 + (z_0i+omega*l)**2)
    # reff = rt + (r*omega**2*m2*cur_l**2)/(r**2*(cur_l+l)**2+omega**2*cur_l**2*l**2)
    return reff
    
def calc_coupled_z(rt, lt, omega, r, l, cur_l, k=1, ls=0, stabilize_tank=True):
    m2 = k**2*lt*l
    m = np.sqrt(m2)
    z_0r = omega**2*cur_l**2*r/(r**2 + omega**2*cur_l**2)
    if stabilize_tank: z_0i = omega*cur_l*r**2/(r**2 + omega**2*cur_l**2)
    else: z_0i = omega*cur_l*r**2/(r**2 + omega**2*cur_l**2) + ls
    # z_0 = z_0r + 1j*z_0i  # complex(z_0r, z_0i)
    z_coupled = 1 / (1/(1j*omega*m) + 1/(z_0r + 1j*(omega*(l-m) + omega*ls + z_0i)))
    # print((z_coupled + 1j*omega*(lt - m)).imag/omega)
    return z_coupled
    
def calc_zeff(reff, leff, ct, omega, mode="series"):
    zl = omega*leff
    zc = 1/(ct*omega)
    if mode == "series":
        real = reff
        imag = zl - zc
    else:  # if parallel
        denom = reff**2 + (zl - zc)**2
        real = reff*zl**2 / denom
        imag = zl*(reff**2 + zc*(zc + zl)) / denom
    z_mag = np.sqrt(real**2 + imag**2)
    return (real, imag, z_mag)

def calc_cur_l(reff, rt, leff, lt, m, omega, l, ls, no_r=False, stabilize_tank=True):
    # only supports ls>0 cases if no_r
    rdiff = reff - rt
    ldiff = leff - lt
    if no_r and not stabilize_tank: cur_l = - m**2 / ldiff - ls - l
    else:
        num = (rdiff)**2*l**2 + omega**2*(m**2+ldiff*l)**2
        denom = -omega**2*ldiff*(m**2 + ldiff*l) - rdiff**2*l
        cur_l = num/denom
    return cur_l
    
def calc_r(reff, rt, leff, lt, m, omega, l, ls, no_r=False):
    # only supports ls>0 cases for no_r
    rdiff = reff - rt
    ldiff = leff - lt
    num = (rdiff)**2*l**2 + omega**2*(m**2+ldiff*l)**2
    if no_r:
        r = rdiff
    else:
        num = (rdiff)**2*l**2 + omega**2*(m**2+ldiff*l)**2
        denom = rdiff*m**2
        r = num/denom
    return r
    
def calc_flat_wire(l, w, h, sigma=5.8e7):
    # Get rectangular wire inductance and resistance
    # all in SI units! copper conductance = 5.8e7 S*m
    R = l/(w*h)/sigma
    L = 2e-7*l*(np.log(2*l/(w+h)) + 0.5 + 0.2235*(w+h)/l)
    return R, L

def calc_k_report(z1, z2, zm, freq):
    omega = 2*np.pi*freq
    l1, l2 = z1/omega, z2/omega
    m = zm/omega
    k = calc_k(l1, l2, m)
    print(f"l1: {l1} H, l2: {l2} H, k: {k}")
    
def calc_capacitance(A, d, k, w=0, h=0):
    if w*h != 0: A = w*h
    return k*cnst.epsilon_0*A/d
    
# fit fns (lmfit)  # model_name.fit(y_vals, x=x_vals, a=0, ...)
def exp_decay(t, a, b, c):
    return a*np.exp(-t/b)+c

def sine_wave(t, a, freq, phase):
    return a*np.sin(2*np.pi*freq*t + phase)

def resistance_phase(phi, r, alpha):
    return r/(1+alpha*np.cos(phi))

def lorentzian(x, par0, par1, par2):  # linear case
    # par0: amplitude, par1: res freq, par2: bw
    return par0 / (1. + pow(((x - par1) / (par2 / 2.)), 2.))

def inverse(x, a, b, c):
    return a/(x-b)+c

def exp_decay_alt(x, a, b, c):
    return np.exp(-a*(x-b)) + c

def exp_inverse(x, a, b, c, d):
    return a*np.exp(-c / (x - b)) + d

exp_decay_model = Model(exp_decay)
exp_alt_model = Model(exp_decay_alt)
sine_model = Model(sine_wave)
rp_model = Model(resistance_phase)
lorentz_model = Model(lorentzian)
inverse_model = Model(inverse)
exp_inverse_model = Model(exp_inverse)

def do_fft(time_, data, round_stepsize=True):
    stepsize = time_[1] - time_[0]  # higher maxfreq when data is finer
    # duration = time_[-1] - time_[0]  # better rbw when data is longer
    if round_stepsize:  # only works when stepsize is 10^n
        round_digit = round(-np.log10(time_[1]-time_[0]))
        stepsize = round(stepsize, round_digit)
    npts = len(time_)  # round(duration/stepsize)+1  # must be integer
    
    yf = np.abs(fft(data)[0:npts//2])  # fftshift(fft(y))
    if npts % 2 == 0: xf = fftfreq(npts, stepsize)[0:npts//2]  # why do i need to do this separately???
    else: xf = fftfreq(npts-1, stepsize)[0:npts//2]
    
    # fft will remove last point in fft because of how fftfreq works
    #if npts == len(data):
    #    xf = fftfreq(npts - 1, stepsize)  # [:-1]
    #else: xf = fftfreq(len(data) - 1, stepsize)
    #yf = fft(data[:-1])  # [:npts//2]
    return xf, yf

'''
def round_value(value):  #, nine_count=3):  # don't do for > 1 yet
    decimals_to_round = 0
    value_mult = value
    # nine_counts = 0  # rounds values like 9.99995
    first_nonzero = False
    finished_rounding = False
    while not finished_rounding:
        # current_digit_value = int(value_mult)%10
        current_digit_value = round(value_mult)%10
        #print("NextDigit", current_digit_value, value_mult, decimals_to_round)
        #if nine_counts > nine_count:
        #    print("TooManyNines")
        #    decimals_to_round -= nine_count
        #    finished_rounding = True
        decimals_to_round += 1
        value_mult = value_mult*10
        if current_digit_value != 0:
            first_nonzero = True
        #    print("FirstNonzero")
        #if current_digit_value == 9:
        #    print("AddNineCount")
        #    nine_counts += 1
        if first_nonzero and current_digit_value == 0:
        #    print("FinishedRounding")
            finished_rounding = True
        #input("...")
    return round(value, decimals_to_round)
'''

'''
l = 5e-10
l0 = 2e-9
k0 = 0.2
l1 = 5e-10
k1 = 0.253144227664784
l2 = 2.1178582977171555e-09
lj = calc_lj(1.6e-7, np.pi)

l_to_cancel = l*(1-k0**2)
ltot_jj = l2 * (1 - k1**2*l1/(l1 + l0*(1 - k0**2)))
ltot_v = l * (1 - k0**2*l0/(l0 + l1*(1 - k1**2 * l2 / (l2 + lj))))
l_to_cancel, ltot_jj, ltot_v

'''
