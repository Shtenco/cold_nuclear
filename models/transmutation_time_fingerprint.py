#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time-fingerprint audit for selected Pd–rare-earth transmutation chains.

This is a normalized nuclear bookkeeping model, not evidence for LENR and not
an apparatus recipe. A reference neutron/event flux is used only to make the
capture -> radioactive-intermediate -> daughter timing explicit.

ODE for each chain:
    dP/dt = -a P
    dI/dt =  a P - b I
    dD/dt =        b I
where a = Phi * sigma and b = ln(2)/half_life.

The 'phase multiplier' is deliberately separated from the nuclear pathway.
A standard near-liquidus reference uses the previously audited factor 0.960348
from number-density and 1/v temperature scaling. An optional 10x hypothesis
scenario is included only as sensitivity, not as a claim.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from pathlib import Path
import pandas as pd

BARN_CM2 = 1e-24
PHI_REF = 1e12  # normalized reference only, n/cm^2/s
STANDARD_LIQUID_OVER_SOLID = 0.9603480186014252

@dataclass(frozen=True)
class Chain:
    name: str
    parent: str
    intermediate: str
    daughter: str
    sigma_b: float
    half_life_h: float
    natural_fraction: float = 1.0

CHAINS = [
    Chain("La_to_Ce", "139La", "140La", "140Ce", 8.93, 40.284, 0.9991),
    Chain("Y_to_Zr", "89Y", "90Y", "90Zr", 1.28, 64.05, 1.0),
    Chain("Ce_to_Pr", "140Ce", "141Ce", "141Pr", 0.57, 32.504*24.0, 0.8848),
    Chain("Pd_to_Ag", "108Pd", "109Pd", "109Ag", 8.55, 13.437, 0.2646),
]

SCENARIOS = {
    "solid_reference": 1.0,
    "liquid_standard_same_source": STANDARD_LIQUID_OVER_SOLID,
    "liquid_hypothesis_source_x10": STANDARD_LIQUID_OVER_SOLID * 10.0,
}

TIMES_H = [0, 6, 12, 24, 48, 72, 7*24, 14*24, 30*24, 60*24, 120*24, 365*24]

def populations(t_s: float, capture_rate_s: float, decay_rate_s: float):
    a = capture_rate_s
    b = decay_rate_s
    P = math.exp(-a*t_s)
    if abs(b-a) < 1e-30:
        I = a*t_s*math.exp(-a*t_s)
    else:
        I = a/(b-a) * (math.exp(-a*t_s) - math.exp(-b*t_s))
    D = max(0.0, 1.0 - P - I)
    return P, I, D

def peak_time_s(a: float, b: float):
    if a <= 0 or b <= 0 or abs(a-b) < 1e-30:
        return math.nan
    return math.log(b/a)/(b-a)

def run(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    summary = []
    for c in CHAINS:
        b = math.log(2.0)/(c.half_life_h*3600.0)
        for scenario, phase_mult in SCENARIOS.items():
            a = PHI_REF * phase_mult * c.sigma_b * BARN_CM2
            tpk = peak_time_s(a, b)
            Ppk, Ipk, Dpk = populations(tpk, a, b)
            summary.append({
                "chain": c.name,
                "parent": c.parent,
                "intermediate": c.intermediate,
                "daughter": c.daughter,
                "scenario": scenario,
                "phase_multiplier": phase_mult,
                "reference_flux_n_cm2_s": PHI_REF,
                "sigma_b": c.sigma_b,
                "half_life_h": c.half_life_h,
                "intermediate_peak_time_days": tpk/86400.0,
                "intermediate_peak_fraction_per_initial_parent": Ipk,
                "daughter_fraction_at_intermediate_peak": Dpk,
            })
            for th in TIMES_H:
                P,I,D = populations(th*3600.0, a, b)
                rows.append({
                    "chain": c.name,
                    "scenario": scenario,
                    "time_h": th,
                    "parent_fraction": P,
                    "intermediate_fraction": I,
                    "daughter_fraction": D,
                    "daughter_fraction_per_natural_element_atom": D*c.natural_fraction,
                })
    pd.DataFrame(rows).to_csv(outdir/"time_fingerprint_curves.csv", index=False)
    sdf = pd.DataFrame(summary)
    sdf.to_csv(outdir/"time_fingerprint_summary.csv", index=False)
    print(sdf.to_string(index=False))

if __name__ == "__main__":
    run(Path(__file__).resolve().parent.parent / "results" / "time_fingerprints")
