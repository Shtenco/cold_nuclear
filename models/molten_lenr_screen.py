#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Molten precious / rare-earth metal extension to LENR MASTER MODEL.

Purpose
-------
This is an exploratory *falsification / sensitivity* model. It does NOT claim
that molten metals cause LENR or D-D fusion. It asks a narrower question:

  If the low-energy D-D screening phenomenology used in the 2026-08-08 master
  model were transplanted from a solid lattice into a liquid metal, what
  combinations of D loading, transient local order, effective screening and
  rare active fraction would be required to produce a specified heat density?

The solid-state "defect" gate is replaced with a liquid-state gate:
  D loading x_D × short-range-order SRO × electronic-response gate.

Material identity enters only through approximate atomic number density and an
EXPERIMENT-PRIORITY heuristic (hydrogen affinity, direct screening evidence,
and temperature penalty). The heuristic is NOT a fusion-rate measurement.

No apparatus geometry, pressure, neutron-source recipe, or operating procedure
is included.
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd

E_CHARGE = 1.602176634e-19
N_A = 6.02214076e23
LN10 = math.log(10.0)
BG = 31.397
P_FULL_REF_W_M3 = 5.19996e16
UE_REF_KEV = 1.7
C_POWER = math.log10(P_FULL_REF_W_M3) + BG / (math.sqrt(UE_REF_KEV) * LN10)
N_D_REF = 6.8e28
Q_DD_AVG_J = 3.65e6 * E_CHARGE

def sigmoid(z):
    z = np.clip(z, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-z))

def p_full_screened(Ue_keV):
    U = np.maximum(np.asarray(Ue_keV, dtype=float), 1e-8)
    log10p = C_POWER - BG / (np.sqrt(U) * LN10)
    return 10.0 ** np.clip(log10p, -300.0, 100.0)

def melt_state_gate(xD, sro, electronic_response):
    gx = sigmoid((np.log10(np.maximum(xD, 1e-12)) + 2.0) / 0.60)
    gs = sigmoid((sro - 0.35) / 0.10)
    ge = sigmoid((electronic_response - 0.45) / 0.10)
    return gx * gs * ge

def melt_power_W_m3(n_metal_m3, xD, sro, electronic_response, Ue_max_keV, f_active):
    G = melt_state_gate(xD, sro, electronic_response)
    Ue_floor = 0.03
    Ue = Ue_floor + (Ue_max_keV - Ue_floor) * G
    nD = n_metal_m3 * xD
    density_scale = (nD / N_D_REF) ** 2
    P_local_ref = p_full_screened(Ue)
    Pavg = f_active * G * density_scale * P_local_ref
    return Pavg, Ue, G, nD

MATERIALS = [
    ("Pd","precious",106.42,12.02,1828,0.75,1.00),
    ("Pt","precious",195.084,21.45,2041,0.25,0.45),
    ("Rh","precious",102.906,12.41,2237,0.25,0.40),
    ("Ir","precious",192.217,22.56,2719,0.18,0.35),
    ("Ru","precious",101.07,12.37,2607,0.25,0.35),
    ("Os","precious",190.23,22.59,3306,0.12,0.25),
    ("Au","precious",196.967,19.32,1337,0.04,0.20),
    ("Ag","precious",107.868,10.49,1235,0.06,0.20),
    ("Sc","rare_earth",44.956,2.985,1814,0.90,0.30),
    ("Y","rare_earth",88.906,4.472,1799,0.98,0.35),
    ("La","rare_earth",138.905,6.145,1193,1.00,0.35),
    ("Ce","rare_earth",140.116,6.770,1068,0.98,0.35),
    ("Pr","rare_earth",140.908,6.773,1208,0.96,0.33),
    ("Nd","rare_earth",144.242,7.007,1297,0.96,0.33),
    ("Sm","rare_earth",150.36,7.52,1345,0.94,0.32),
    ("Eu","rare_earth",151.964,5.244,1099,0.90,0.30),
    ("Gd","rare_earth",157.25,7.90,1585,0.96,0.34),
    ("Tb","rare_earth",158.925,8.23,1629,0.95,0.33),
    ("Dy","rare_earth",162.50,8.55,1680,0.95,0.33),
    ("Ho","rare_earth",164.930,8.795,1734,0.94,0.32),
    ("Er","rare_earth",167.259,9.066,1802,0.94,0.32),
    ("Tm","rare_earth",168.934,9.321,1818,0.92,0.31),
    ("Yb","rare_earth",173.045,6.90,1097,0.88,0.30),
    ("Lu","rare_earth",174.967,9.841,1925,0.92,0.31),
]

def material_df():
    cols=["symbol","family","atomic_mass_g_mol","density_g_cm3","Tm_K","H_affinity_heuristic","screen_evidence_heuristic"]
    df=pd.DataFrame(MATERIALS,columns=cols)
    df["liquid_density_proxy_g_cm3"] = 0.90*df.density_g_cm3
    df["metal_atom_density_m3"] = df.liquid_density_proxy_g_cm3*1e6/df.atomic_mass_g_mol*N_A
    df["temperature_retention_proxy"] = np.exp(-(df.Tm_K-1000.0)/2500.0).clip(0.20,1.0)
    df["pure_material_priority_score"] = (
        0.50*df.H_affinity_heuristic + 0.35*df.screen_evidence_heuristic + 0.15*df.temperature_retention_proxy
    )
    return df

def required_factive(n_metal, xD, sro, er, Ue_max, target_W_cm3):
    G=melt_state_gate(xD,sro,er)
    Ue=0.03+(Ue_max-0.03)*G
    nD=n_metal*xD
    base=G*(nD/N_D_REF)**2*float(p_full_screened(Ue))
    target=target_W_cm3*1e6
    return np.inf if base<=0 else target/base

def pure_threshold_table(mdf):
    rows=[]
    for _,m in mdf.iterrows():
        for xD in [1e-4,1e-3,1e-2,1e-1]:
            for Ue in [0.3,0.5,0.9,1.7]:
                fa_1mW=required_factive(m.metal_atom_density_m3,xD,0.60,0.70,Ue,1e-3)
                fa_1W=required_factive(m.metal_atom_density_m3,xD,0.60,0.70,Ue,1.0)
                rows.append(dict(symbol=m.symbol,family=m.family,Tm_K=m.Tm_K,x_D=xD,Ue_max_keV=Ue,
                                 f_active_required_1mW_cm3=fa_1mW,
                                 f_active_required_1W_cm3=fa_1W,
                                 pure_material_priority_score=m.pure_material_priority_score))
    return pd.DataFrame(rows)

def ideal_binary(a,b,xB):
    xA=1-xB
    VA=a.atomic_mass_g_mol/a.density_g_cm3/0.90
    VB=b.atomic_mass_g_mol/b.density_g_cm3/0.90
    Vmix=xA*VA+xB*VB
    Mmix=xA*a.atomic_mass_g_mol+xB*b.atomic_mass_g_mol
    rho=Mmix/Vmix
    n=N_A/(Vmix*1e-6)
    Tproxy=xA*a.Tm_K+xB*b.Tm_K
    return rho,n,Tproxy

def binary_screen(mdf):
    precious=mdf[mdf.family=="precious"]
    rare=mdf[mdf.family=="rare_earth"]
    rows=[]
    for _,p in precious.iterrows():
        for _,r in rare.iterrows():
            for xr in [0.25,0.50,0.75]:
                rho,n,T=ideal_binary(p,r,xr)
                comp = math.sqrt(max(r.H_affinity_heuristic,1e-9)*max(p.screen_evidence_heuristic,1e-9))
                sro_bonus = math.exp(-((xr-0.50)/0.24)**2)
                temp = max(0.20,min(1.0,math.exp(-(T-1000)/2500)))
                priority = 0.45*r.H_affinity_heuristic + 0.35*p.screen_evidence_heuristic + 0.10*comp*sro_bonus + 0.10*temp
                name=f"{p.symbol}{int(round((1-xr)*100))}-{r.symbol}{int(round(xr*100))}"
                fa=required_factive(n,1e-2,0.60,0.70,0.9,1e-3)
                rows.append(dict(alloy=name,precious=p.symbol,rare_earth=r.symbol,x_RE=xr,
                                 ideal_density_proxy_g_cm3=rho,ideal_Tmix_proxy_K=T,
                                 metal_atom_density_m3=n,experimental_priority_score=priority,
                                 f_active_required_for_1mW_cm3_at_xD1pct_Ue0p9=fa))
    return pd.DataFrame(rows).sort_values(["experimental_priority_score"],ascending=False)

def generic_melt_scan(N=400_000,seed=20260808):
    rng=np.random.default_rng(seed)
    xD=10**rng.uniform(-6,-0.3,N)
    sro=rng.uniform(0.02,0.95,N)
    er=rng.uniform(0.05,0.95,N)
    Umax=rng.uniform(0.1,1.9,N)
    f=10**rng.uniform(-16,-2,N)
    n=6.0e28
    P,U,G,nD=melt_power_W_m3(n,xD,sro,er,Umax,f)
    Pcm=P/1e6
    Nn=0.5*Pcm/Q_DD_AVG_J
    df=pd.DataFrame(dict(x_D=xD,short_range_order=sro,electronic_response=er,Ue_max_keV=Umax,
                         Ue_effective_keV=U,state_gate_G=G,f_active=f,
                         power_W_cm3=Pcm,standard_DD_neutron_rate_n_s_cm3=Nn))
    cand=df[(df.power_W_cm3>=1e-3)&(df.power_W_cm3<=10)].copy()
    return df,cand

def main():
    out=Path(__file__).resolve().parent.parent/"results"/"molten"
    out.mkdir(parents=True,exist_ok=True)
    m=material_df()
    m.to_csv(out/"material_properties_and_priority.csv",index=False)
    pure=pure_threshold_table(m)
    pure.to_csv(out/"pure_metal_thresholds.csv",index=False)
    alloys=binary_screen(m)
    alloys.to_csv(out/"precious_rareearth_binary_screen.csv",index=False)
    full,cand=generic_melt_scan()
    full.sample(min(80000,len(full)),random_state=7).to_csv(out/"melt_parameter_scan_sample.csv",index=False)
    cand.sort_values("power_W_cm3").head(30000).to_csv(out/"melt_candidate_region.csv",index=False)
    top_pure=m.sort_values("pure_material_priority_score",ascending=False).head(12)
    top_alloy=alloys.head(30)
    top_pure[["symbol","family","Tm_K","pure_material_priority_score"]].to_csv(out/"TOP_pure_materials.csv",index=False)
    top_alloy.to_csv(out/"TOP_30_binary_alloys.csv",index=False)
    print("MOLTEN LENR SCREEN")
    print("scan points:",len(full),"candidate 1mW..10W/cm3:",len(cand))
    print(top_alloy.head(20)[["alloy","experimental_priority_score","ideal_Tmix_proxy_K","f_active_required_for_1mW_cm3_at_xD1pct_Ue0p9"]].to_string(index=False))

if __name__=="__main__":
    main()
