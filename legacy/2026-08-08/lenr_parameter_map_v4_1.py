#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LENR parameter map v4.1 — corrected energy-balance treatment.

Differences from v4:
- External drive is NOT randomly sampled.
- The model reports P_out and therefore the break-even ceiling:
      P_in,max(Q>=1) = P_out
- A real Q is only obtained after inserting measured external input power.
- Standard D-D neutron production is calculated independently.
- Generic neutron upper-limit benchmarks are converted into the MAXIMUM
  standard-D-D heat power compatible with those limits.

No apparatus construction settings are provided.
"""

import numpy as np
import pandas as pd
import math
from pathlib import Path

E_CHARGE = 1.602176634e-19
BG = 31.397
UE_REF = 1.7
P_FULL_REF = 5.19996e16
Q_DD_AVG_J = 3.65e6 * E_CHARGE
LN10 = math.log(10.0)
C_POWER = math.log10(P_FULL_REF) + BG / (math.sqrt(UE_REF) * LN10)

def sigmoid(z):
    z = np.clip(z, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-z))

def p_full_screened(Ue_keV):
    U = np.maximum(np.asarray(Ue_keV, dtype=float), 1e-8)
    log10p = C_POWER - BG / (np.sqrt(U) * LN10)
    return 10.0 ** np.clip(log10p, -300.0, 100.0)

def state_gate(x_load, defect, detune, Qmode=10.0):
    gx = sigmoid((x_load - 0.78) / 0.035)
    gd = sigmoid((defect - 0.28) / 0.045)
    resonance = 1.0 / np.sqrt(1.0 + (2.0 * Qmode * detune)**2)
    gr = sigmoid((resonance - 0.42) / 0.08)
    return gx * gd * gr, resonance

def branchB_power(x_load, defect, detune, Ue_max, fV, ft):
    G, resonance = state_gate(x_load, defect, detune)
    Ue_floor = 0.03
    Ue = Ue_floor + (Ue_max - Ue_floor) * G
    local = p_full_screened(Ue)
    Pavg_m3 = fV * ft * G * local
    focus = 1.0 / np.maximum(fV * ft, 1e-300)
    return Pavg_m3, Ue, G, resonance, focus

def neutron_rate_per_cm3(P_W_cm3):
    return 0.5 * P_W_cm3 / Q_DD_AVG_J

def scan(N=300_000, seed=20260808):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.55, 0.995, N)
    defect = rng.uniform(0.02, 0.80, N)
    detune = rng.uniform(0.0, 0.15, N)
    Umax = rng.uniform(0.5, 1.9, N)
    fV = 10.0 ** rng.uniform(-14.0, -2.0, N)
    ft = 10.0 ** rng.uniform(-12.0, -1.0, N)
    Pm3, Ue, G, resonance, focus = branchB_power(x, defect, detune, Umax, fV, ft)
    Pcm3 = Pm3 / 1e6
    Nstd = neutron_rate_per_cm3(Pcm3)
    limits = [1e3, 1e6, 1e9, 1e12]
    beta = {lim: np.minimum(1.0, lim / np.maximum(Nstd, 1e-300)) for lim in limits}
    df = pd.DataFrame({
        "x_D_loading": x,
        "defect_state": defect,
        "normalized_detuning": detune,
        "resonance_response": resonance,
        "Ue_max_keV": Umax,
        "Ue_effective_keV": Ue,
        "state_gate_G": G,
        "active_volume_fraction_fV": fV,
        "active_time_fraction_ft": ft,
        "effective_active_fraction_fVftG": fV*ft*G,
        "focus_factor_1_over_fVft": focus,
        "predicted_avg_power_W_cm3": Pcm3,
        "max_avg_input_for_Q_ge_1_W_cm3": Pcm3,
        "standard_DD_neutron_rate_n_s_cm3": Nstd,
        "required_n_branch_fraction_if_Nlim_1e3": beta[1e3],
        "required_n_branch_fraction_if_Nlim_1e6": beta[1e6],
        "required_n_branch_fraction_if_Nlim_1e9": beta[1e9],
        "required_n_branch_fraction_if_Nlim_1e12": beta[1e12],
    })
    cand = df[(df.predicted_avg_power_W_cm3 >= 1e-3) & (df.predicted_avg_power_W_cm3 <= 10.0)].copy()
    out = Path(__file__).resolve().parents[2] / "results" / "master"
    out.mkdir(parents=True, exist_ok=True)
    df.sample(min(60_000, len(df)), random_state=1).to_csv(out/"lenr_parameter_scan_v4_1.csv", index=False)
    cand.sort_values("predicted_avg_power_W_cm3").head(20_000).to_csv(out/"lenr_candidate_region_v4_1.csv", index=False)
    nb=[]
    n_per_W=0.5/Q_DD_AVG_J
    for lim in limits:
        pmax=lim/n_per_W
        nb.append({"neutron_upper_limit_n_s_cm3":lim,"max_standard_DD_heat_W_cm3":pmax,"max_standard_DD_heat_mW_cm3":pmax*1e3})
    pd.DataFrame(nb).to_csv(out/"lenr_neutron_budget_v4_1.csv",index=False)
    print(f"Total points: {len(df):,}; 1mW..10W/cm3 candidates: {len(cand):,}")
    return df,cand

if __name__ == "__main__":
    scan()
