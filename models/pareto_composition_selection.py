#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pareto selection across Pd–RE compositions (25/75, 50/50, 75/25).

Objectives are maximized independently; no single weighted score determines the
Pareto front. A geometric balance score is included only as a secondary display
metric and must not be interpreted as physics.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MOLTEN = ROOT / "results" / "molten" / "TOP_30_binary_alloys.csv"
TRANS_ELEMENT = ROOT / "results" / "transmutation" / "element_transmutation_efficiency.csv"
OUT = ROOT / "results" / "selection"

PARTNERS = ["La","Ce","Y","Nd","Gd"]
OBJECTIVES = ["molten_priority","sigma_element_change_b_per_atom","change_fraction"]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mol = pd.read_csv(MOLTEN)
    el = pd.read_csv(TRANS_ELEMENT).set_index("element")
    sub = mol[(mol.precious=="Pd") & (mol.rare_earth.isin(PARTNERS))].copy()
    rows=[]
    for _,r in sub.iterrows():
        xre=float(r.x_RE); xpd=1.0-xre; re=r.rare_earth
        sa=xpd*el.loc["Pd","sigma_abs_natural_b"] + xre*el.loc[re,"sigma_abs_natural_b"]
        sc=xpd*el.loc["Pd","sigma_element_change_b"] + xre*el.loc[re,"sigma_element_change_b"]
        rows.append({
            "alloy":r.alloy,"partner":re,"x_Pd":xpd,"x_RE":xre,
            "molten_priority":r.experimental_priority_score,
            "sigma_abs_b_per_atom":sa,
            "sigma_element_change_b_per_atom":sc,
            "change_fraction":sc/sa,
            "change_percent":100.0*sc/sa,
            "f_active_req_1mW_xD1pct_U0p9":r.f_active_required_for_1mW_cm3_at_xD1pct_Ue0p9,
            "ideal_Tmix_proxy_K":r.ideal_Tmix_proxy_K,
        })
    df=pd.DataFrame(rows)
    front=[]; doms=[]
    for i,row in df.iterrows():
        ds=[]
        for j,o in df.iterrows():
            if i==j: continue
            if all(o[c]>=row[c] for c in OBJECTIVES) and any(o[c]>row[c] for c in OBJECTIVES):
                ds.append(o.alloy)
        front.append(len(ds)==0); doms.append(";".join(ds))
    df["pareto_front"]=front; df["dominated_by"]=doms
    normcols=[]
    for c in OBJECTIVES:
        mn,mx=df[c].min(),df[c].max()
        nc=c+"_norm"; normcols.append(nc)
        df[nc]=(df[c]-mn)/(mx-mn) if mx>mn else 1.0
    df["geometric_balance_score"]=(df[normcols].clip(lower=1e-12).prod(axis=1))**(1/3)
    keep=["alloy","partner","x_Pd","x_RE","molten_priority","change_percent",
          "sigma_element_change_b_per_atom","sigma_abs_b_per_atom",
          "f_active_req_1mW_xD1pct_U0p9","ideal_Tmix_proxy_K",
          "pareto_front","dominated_by","geometric_balance_score"]
    out=df[keep].sort_values(["pareto_front","geometric_balance_score"],ascending=[False,False])
    out.to_csv(OUT/"pd_re_composition_pareto.csv",index=False)
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
