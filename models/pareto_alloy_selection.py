#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pareto selection across molten-material and transmutation objectives.

No arbitrary weighted average is used. A candidate is Pareto-dominated if
another candidate is at least as good in every objective and strictly better
in at least one.

Objectives (all maximized):
1) molten experimental-priority heuristic;
2) percent of absorptions that change chemical element within the audit horizon;
3) effective element-change cross section per alloy atom.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MOLTEN = ROOT / "results" / "molten" / "TOP_30_binary_alloys.csv"
TRANS = ROOT / "results" / "transmutation" / "alloy_50_50_ranking.csv"
OUT = ROOT / "results" / "selection"

OBJECTIVES = [
    "experimental_priority_score",
    "percent_absorptions_that_change_element_3y",
    "effective_element_change_b_per_atom",
]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mol = pd.read_csv(MOLTEN)
    tr = pd.read_csv(TRANS)
    df = mol[mol.alloy.isin(tr.alloy)].merge(tr, on="alloy")
    front=[]; dominators=[]
    for i,row in df.iterrows():
        dom=[]
        for j,other in df.iterrows():
            if i == j:
                continue
            if all(other[c] >= row[c] for c in OBJECTIVES) and any(other[c] > row[c] for c in OBJECTIVES):
                dom.append(other.alloy)
        front.append(len(dom)==0)
        dominators.append(";".join(dom))
    df["pareto_front"] = front
    df["dominated_by"] = dominators
    keep = ["alloy"] + OBJECTIVES + [
        "effective_absorption_b_per_atom",
        "fraction_captures_on_Pd",
        "fraction_captures_on_partner",
        "pareto_front",
        "dominated_by",
    ]
    result = df[keep].sort_values(OBJECTIVES, ascending=False)
    result.to_csv(OUT/"pd_re_pareto.csv", index=False)
    print(result.to_string(index=False))

if __name__ == "__main__":
    main()
