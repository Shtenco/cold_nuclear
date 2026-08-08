#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D-D solid-state / screened-fusion consistency model.

Purpose
-------
1) Reproduce UKAEA PROCESS Bosch-Hale D-D thermal reactivity.
2) Calibrate a simple low-energy Gamow/S-factor cross-section model to Bosch-Hale.
3) Apply the experimentally reported effective screening potentials from:
   Karahadian et al., Nature Communications (2026), DOI 10.1038/s41467-026-74421-1.
4) Check whether the model reproduces the order of magnitude of sub-keV enhancement.
5) Stress-test the (invalid) extrapolation of the beam-derived screening potential to
   room-temperature lattice deuterons.
6) Quantify THz phonon energies using exact SI constants.
"""

import math
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

E_CHARGE = 1.602176634e-19
H_PLANCK = 6.62607015e-34
K_B = 1.380649e-23
C_LIGHT = 299792458.0
KEV_J = 1e3 * E_CHARGE
BARN = 1e-28

DD1 = dict(bg=31.397, mrc2=937814.0, cc1=5.4336e-12, cc2=0.00585778, cc3=0.00768222, cc4=0.0, cc5=-2.964e-06, cc6=0.0, cc7=0.0)
DD2 = dict(bg=31.397, mrc2=937814.0, cc1=5.65718e-12, cc2=0.00341267, cc3=0.00199167, cc4=0.0, cc5=1.0506e-05, cc6=0.0, cc7=0.0)
BG = DD1["bg"]
MU = DD1["mrc2"] * KEV_J / C_LIGHT**2

def bosch_hale_reactivity(T_keV, r):
    T = np.asarray(T_keV, dtype=float)
    theta = T / (1.0 - T * (r["cc2"] + T * (r["cc4"] + T * r["cc6"])) / (1.0 + T * (r["cc3"] + T * (r["cc5"] + T * r["cc7"]))))
    xi = (r["bg"]**2 / (4.0 * theta)) ** (1.0 / 3.0)
    return 1e-6 * r["cc1"] * theta * np.sqrt(xi / (r["mrc2"] * T**3)) * np.exp(-3.0 * xi)

def bh_total(T_keV):
    return bosch_hale_reactivity(T_keV, DD1) + bosch_hale_reactivity(T_keV, DD2)

def sigma_gamow(E_keV, S0_keV_b, Ue_keV=0.0):
    E = np.asarray(E_keV, dtype=float)
    return (S0_keV_b / E) * np.exp(-BG / np.sqrt(E + Ue_keV)) * BARN

def maxwell_reactivity_gamow(T_keV, S0_keV_b):
    kT_J = T_keV * KEV_J
    pref = math.sqrt(8.0 / (math.pi * MU)) / kT_J**1.5
    E0 = (BG * T_keV / 2.0) ** (2.0 / 3.0)
    logE0 = math.log(E0)
    def phi(E): return -BG / math.sqrt(E) - E / T_keV
    phi0 = phi(E0)
    def scaled_integrand_logx(x):
        E = math.exp(x)
        return math.exp(phi(E) - phi0) * E
    val = quad(scaled_integrand_logx, logE0 - 30.0, logE0 + 30.0, epsabs=0.0, epsrel=1e-9, limit=500)[0]
    integral = S0_keV_b * BARN * math.exp(phi0) * val
    return pref * KEV_J**2 * integral

def room_temp_log_reactivity(T_K, S0_keV_b):
    T_keV = K_B * T_K / KEV_J
    E0 = (BG * T_keV / 2.0) ** (2.0 / 3.0)
    logE0 = math.log(E0)
    def phi(E): return -BG / math.sqrt(E) - E / T_keV
    phi0 = phi(E0)
    def scaled_integrand_logx(x):
        E = math.exp(x)
        return math.exp(phi(E) - phi0) * E
    val = quad(scaled_integrand_logx, logE0 - 35.0, logE0 + 35.0, epsabs=0.0, epsrel=1e-9, limit=500)[0]
    log_integral = math.log(S0_keV_b * BARN) + phi0 + math.log(val)
    kT_J = T_keV * KEV_J
    log_pref = 0.5 * math.log(8.0 / (math.pi * MU)) - 1.5 * math.log(kT_J) + 2.0 * math.log(KEV_J)
    return log_pref + log_integral, T_keV, E0

def naive_static_screening_log_reactivity(T_K, S0_keV_b, Ue_keV):
    T_keV = K_B * T_K / KEV_J
    kT_J = T_keV * KEV_J
    pref_log = 0.5 * math.log(8.0 / (math.pi * MU)) - 1.5 * math.log(kT_J) + 2.0 * math.log(KEV_J)
    def psi(E): return -BG / math.sqrt(E + Ue_keV) - E / T_keV
    def scaled_y(y):
        E = T_keV * y
        return math.exp(psi(E) - psi(0.0))
    val = quad(scaled_y, 0.0, 100.0, epsabs=0.0, epsrel=1e-10)[0] * T_keV
    log_integral = math.log(S0_keV_b * BARN) + psi(0.0) + math.log(val)
    return pref_log + log_integral

def main():
    fit_T = np.array([0.20, 0.25, 0.50, 1.00])
    unit = np.array([maxwell_reactivity_gamow(T, 1.0) for T in fit_T])
    target = bh_total(fit_T)
    S0 = math.exp(np.mean(np.log(target / unit)))
    print(f"Calibrated total low-energy S0 = {S0:.3f} keV*b")
    for T in [0.2, 0.25, 0.5, 1, 2, 5, 10]:
        bh = float(bh_total(T)); gm = maxwell_reactivity_gamow(T, S0)
        print(T, bh, gm, gm/bh)
    for E in [0.25,0.5,1,2,2.5,3,6.5]:
        bare=float(sigma_gamow(E,S0,0.0))
        print("screen",E,float(sigma_gamow(E,S0,1.7)/bare),float(sigma_gamow(E,S0,0.9)/bare))
    nD=6.8e28; Qavg_J=3.65e6*E_CHARGE
    log_sv,T_keV_300,E0_300=room_temp_log_reactivity(300.0,S0)
    log_R=math.log(0.5*nD**2)+log_sv; log_P=log_R+math.log(Qavg_J)
    print("bare 300K log10 P",log_P/math.log(10))
    for U in [0.9,1.7]:
        log_sv_s=naive_static_screening_log_reactivity(300.0,S0,U)
        log_R_s=math.log(0.5*nD**2)+log_sv_s; log_P_s=log_R_s+math.log(Qavg_J)
        print("UNPHYSICAL stress",U,log_P_s/math.log(10))

if __name__ == "__main__":
    main()
