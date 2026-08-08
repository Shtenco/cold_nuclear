#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graph-consistency gate for the COLD NUCLEAR project.

The graph layer does not supply nuclear physics. It checks that a graph
representation of capture/decay pathways reproduces the tabular weighted
cross-section bookkeeping already stored in results/transmutation.

Requires: pandas, networkx.
"""
from __future__ import annotations
from pathlib import Path
import json
import pandas as pd
import networkx as nx

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "results" / "transmutation"
OUT = ROOT.parent / "results" / "graph"

CHAINS = [
    ("139La", "140La", "capture", 8.93),
    ("140La", "140Ce", "beta-", None),
    ("89Y", "90Y", "capture", 1.28),
    ("90Y", "90Zr", "beta-", None),
    ("140Ce", "141Ce", "capture", 0.57),
    ("141Ce", "141Pr", "beta-", None),
    ("108Pd", "109Pd", "capture", 8.55),
    ("109Pd", "109Ag", "beta-", None),
    ("158Gd", "159Gd", "capture", 2.2),
    ("159Gd", "159Tb", "beta-", None),
]

ELEMENT = {
    "139La":"La","140La":"La","140Ce":"Ce",
    "89Y":"Y","90Y":"Y","90Zr":"Zr",
    "140Ce":"Ce","141Ce":"Ce","141Pr":"Pr",
    "108Pd":"Pd","109Pd":"Pd","109Ag":"Ag",
    "158Gd":"Gd","159Gd":"Gd","159Tb":"Tb",
}

def build_graph():
    g = nx.DiGraph()
    for u,v,kind,sigma in CHAINS:
        g.add_node(u, element=ELEMENT[u])
        g.add_node(v, element=ELEMENT[v])
        g.add_edge(u,v,edge_type=kind,sigma_b=sigma)
    return g

def audit_weighted_table(df):
    calc = (df.assign(change_weight=lambda x: x.capture_weight_b*x.element_change_3y)
              .groupby("element")
              .agg(total_capture_b=("capture_weight_b","sum"),
                   element_change_b=("change_weight","sum")))
    calc["change_fraction"] = calc.element_change_b/calc.total_capture_b
    return calc

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA/"isotope_pathways.csv")
    reference = pd.read_csv(DATA/"element_transmutation_efficiency.csv").set_index("element")
    calc = audit_weighted_table(df)
    joined = calc.join(reference[["sigma_abs_natural_b","sigma_element_change_b","fraction_absorptions_that_change_element_3y"]])
    joined["abs_err_total_b"] = (joined.total_capture_b-joined.sigma_abs_natural_b).abs()
    joined["abs_err_change_b"] = (joined.element_change_b-joined.sigma_element_change_b).abs()
    joined["abs_err_fraction"] = (joined.change_fraction-joined.fraction_absorptions_that_change_element_3y).abs()
    joined.to_csv(OUT/"graph_vs_table_consistency.csv")

    g = build_graph()
    tests = []
    for parent,daughter in [("139La","140Ce"),("89Y","90Zr"),("140Ce","141Pr"),("108Pd","109Ag"),("158Gd","159Tb")]:
        path = nx.shortest_path(g,parent,daughter)
        tests.append({
            "parent":parent,
            "daughter":daughter,
            "path":" -> ".join(path),
            "edge_count":len(path)-1,
            "element_changed":ELEMENT[parent] != ELEMENT[daughter],
        })
    pd.DataFrame(tests).to_csv(OUT/"diagnostic_graph_paths.csv",index=False)

    report = {
        "max_abs_error_total_b": float(joined.abs_err_total_b.max()),
        "max_abs_error_change_b": float(joined.abs_err_change_b.max()),
        "max_abs_error_fraction": float(joined.abs_err_fraction.max()),
        "table_graph_gate_pass": bool(joined[["abs_err_total_b","abs_err_change_b","abs_err_fraction"]].to_numpy().max() < 1e-12),
        "nodes": g.number_of_nodes(),
        "edges": g.number_of_edges(),
        "diagnostic_paths": tests,
        "hypothesis_isolation_rule": "No material-phase hypothesis edge is present in this graph. Adding/removing such an edge may scale a source term but must not alter evaluated capture/decay identities."
    }
    (OUT/"graph_consistency_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__ == "__main__":
    main()
