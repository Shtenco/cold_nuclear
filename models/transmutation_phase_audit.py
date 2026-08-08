#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRANSMUTATION AUDIT: Pd–La / Pd–Ce / Pd–Y / Pd–Nd / Pd–Gd

Purpose
-------
Separate two questions that were previously mixed:
1) If a neutron capture occurs, does the nucleus eventually become a DIFFERENT
   chemical element (Z changes) on a ~3-year diagnostic horizon?
2) Can crossing a solid/liquid phase boundary change the transmutation RATE?

Important:
- This is a nuclear-data accounting model, not evidence for LENR.
- No neutron-source construction or operating recipe is specified.
- Nuclear pathways are determined by isotope + reaction + decay.
  A liquid phase can alter target density, neutron spectrum, or a hypothetical
  event-source term, but it does not turn an otherwise forbidden A->B path
  into an allowed one merely by melting the metal.

Data anchors
------------
NIST NCNR: neutron absorption cross sections at 2200 m/s and natural abundances.
NNDC ENSDF / Nuclear Wallet Cards: decay products and half-lives.
"""

from pathlib import Path
import pandas as pd, math

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / "results" / "transmutation"

def load():
    return pd.read_csv(OUT/"isotope_pathways.csv")

def element_summary(df):
    g = df.groupby("element").agg(
        sigma_abs_natural_b=("capture_weight_b","sum"),
        sigma_element_change_b=("element_change_weight_b","sum"),
    )
    g["change_fraction"] = g.sigma_element_change_b/g.sigma_abs_natural_b
    return g

def alloy_5050(g, partner):
    sa = 0.5*g.loc["Pd","sigma_abs_natural_b"] + 0.5*g.loc[partner,"sigma_abs_natural_b"]
    sc = 0.5*g.loc["Pd","sigma_element_change_b"] + 0.5*g.loc[partner,"sigma_element_change_b"]
    return {
        "alloy":f"Pd50-{partner}50",
        "sigma_abs_b_per_atom":sa,
        "sigma_change_b_per_atom":sc,
        "change_fraction":sc/sa,
        "pd_capture_share":0.5*g.loc["Pd","sigma_abs_natural_b"]/sa,
    }

def phase_ratio(number_density_ratio=0.97, Tsolid_over_TL=0.99,
                Tliquid_over_TL=1.01, source_ratio=1.0):
    """
    Simple sensitivity only.
    For a 1/v-like absorption cross section:
       rate_liquid/rate_solid
       ~= (Nliq/Nsol)*sqrt(Tsol/Tliq)*(Sliq/Ssol)
    The last term S is exactly where a hypothetical phase-dependent LENR
    event source would enter. Standard nuclear data do not provide S.
    """
    return (number_density_ratio
            * math.sqrt(Tsolid_over_TL/Tliquid_over_TL)
            * source_ratio)

if __name__ == "__main__":
    df = load()
    g = element_summary(df)
    print("\nELEMENT-CHANGE EFFICIENCY PER ABSORPTION")
    print((100*g.change_fraction).sort_values(ascending=False).round(6).astype(str)+" %")
    print("\n50:50 Pd-RE ALLOYS")
    rows = [alloy_5050(g,x) for x in ["La","Ce","Y","Nd","Gd"]]
    print(pd.DataFrame(rows).sort_values("change_fraction",ascending=False).to_string(index=False))
    print("\nREFERENCE PHASE RATIO (same source, Nliq/Nsol=0.97):",
          phase_ratio())
