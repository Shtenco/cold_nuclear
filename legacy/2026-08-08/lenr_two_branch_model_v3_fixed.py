#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-branch LENR / materials-enhanced D-D model.

BRANCH A preserves the earlier active-fraction scaling.
BRANCH B is a phenomenological transient state model.

IMPORTANT: Branch B is NOT a validated microscopic theory. It is retained as a
falsifiable state-space hypothesis. This archived copy fixes the missing
`matplotlib.pyplot` import found in the original saved source.
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

E_CHARGE = 1.602176634e-19
BG = 31.397
P_FULL_REF = 5.19996e16
UE_REF = 1.7
Q_DD_AVG_J = 3.65e6 * E_CHARGE
LOG10E = math.log(10.0)
C_POWER = math.log10(P_FULL_REF) + BG / (math.sqrt(UE_REF) * LOG10E)

def sigmoid(z):
    z = np.clip(z, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-z))

def p_full_screened(Ue_keV):
    U = np.maximum(np.asarray(Ue_keV, dtype=float), 1e-8)
    log10p = C_POWER - BG / (np.sqrt(U) * LOG10E)
    return 10.0 ** np.clip(log10p, -300.0, 100.0)

def branch_A(eta, Ue_keV=1.7):
    return eta * float(p_full_screened(Ue_keV))

def pulse_train(t, start=30.0, stop=80.0, period=0.50, width=0.08):
    if t < start or t >= stop:
        return 0.0
    return 1.0 if ((t - start) % period) < width else 0.0

def equal_energy_continuous_drive(t, start=30.0, stop=80.0, duty=0.08/0.50):
    return duty if start <= t < stop else 0.0

def simulate(pulsed=True, dt=0.01, t_end=120.0):
    t = np.arange(0.0, t_end + dt, dt)
    n = len(t)
    x = np.zeros(n); d = np.zeros(n); a = np.zeros(n); Ue = np.zeros(n)
    eta = np.zeros(n); P = np.zeros(n); drive_arr = np.zeros(n)
    x[0] = 0.20; d[0] = 0.05

    k_load = 0.090; k_unload = 0.004; k_pulse_desorb = 0.012
    k_def_transport = 0.70; k_def_drive = 0.055; tau_heal = 55.0
    tau_a = 0.12; resonance_gain = 5.0
    x_c, wx = 0.78, 0.035
    d_c, wd = 0.28, 0.045
    a_c, wa = 2.00, 0.25
    Ue_floor = 0.03; Ue_max = 1.70
    eta_max = 3.0e-11

    for i in range(1, n):
        ti = t[i-1]
        load_on = 1.0 if ti < 80.0 else 0.0
        drive = pulse_train(ti) if pulsed else equal_energy_continuous_drive(ti)
        drive_arr[i-1] = drive
        dx = (k_load * load_on * (1.0 - x[i-1])
              - k_unload * (1.0 - load_on) * x[i-1]
              - k_pulse_desorb * drive * x[i-1])
        dd = (k_def_transport * abs(dx) * (1.0 - d[i-1])
              + k_def_drive * drive * (1.0 - d[i-1])
              - d[i-1] / tau_heal)
        da = (resonance_gain * drive - a[i-1]) / tau_a
        x[i] = np.clip(x[i-1] + dt * dx, 0.0, 1.0)
        d[i] = np.clip(d[i-1] + dt * dd, 0.0, 1.0)
        a[i] = max(0.0, a[i-1] + dt * da)
        gx = sigmoid((x[i] - x_c) / wx)
        gd = sigmoid((d[i] - d_c) / wd)
        ga = sigmoid((a[i] - a_c) / wa)
        overlap = gx * gd * ga
        Ue[i] = Ue_floor + (Ue_max - Ue_floor) * overlap
        eta[i] = eta_max * overlap
        P[i] = eta[i] * p_full_screened(Ue[i])

    drive_arr[-1] = drive_arr[-2]
    rxn_rate = P / Q_DD_AVG_J
    neutron_rate = 0.5 * rxn_rate
    return pd.DataFrame({
        "t_s": t, "drive": drive_arr, "D_loading_x": x,
        "defect_state_d": d, "mode_envelope_a": a,
        "Ue_keV": Ue, "active_fraction_eta": eta,
        "power_W_m3": P, "reaction_rate_m3_s": rxn_rate,
        "neutron_rate_m3_s": neutron_rate,
    })

def summarize(df, label):
    dt = float(df.t_s.iloc[1] - df.t_s.iloc[0])
    peak_idx = int(df.power_W_m3.values.argmax())
    peak = df.iloc[peak_idx]
    E_J_m3 = float(np.trapezoid(df.power_W_m3.values, df.t_s.values))
    active_time = float(np.sum(df.active_fraction_eta.values > 1e-12) * dt)
    return {
        "branch": label,
        "peak_time_s": peak.t_s,
        "peak_power_W_m3": peak.power_W_m3,
        "peak_power_W_cm3": peak.power_W_m3 / 1e6,
        "peak_Ue_keV": peak.Ue_keV,
        "peak_eta": peak.active_fraction_eta,
        "peak_D_loading": peak.D_loading_x,
        "peak_defect_state": peak.defect_state_d,
        "peak_mode_envelope": peak.mode_envelope_a,
        "integrated_energy_J_m3": E_J_m3,
        "time_eta_gt_1e-12_s": active_time,
        "peak_neutron_rate_m3_s_if_standard_DD": peak.neutron_rate_m3_s,
    }

def main():
    branchA_rows=[]
    for label,eta in [("A_eta_1e-14",1e-14),("A_eta_1e-12",1e-12),("A_eta_1e-11",1e-11),("A_eta_1e-10",1e-10)]:
        p=branch_A(eta,1.7)
        branchA_rows.append({"branch":label,"peak_power_W_m3":p,"peak_power_W_cm3":p/1e6,"peak_eta":eta})
    pulsed=simulate(True); control=simulate(False)
    summary=pd.DataFrame(branchA_rows+[summarize(pulsed,"B_pulsed_staged"),summarize(control,"B_equal_energy_continuous")])
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
