# COLD NUCLEAR — reproducible computation repository

**Status:** computational research programme. **Not proof of LENR, cold fusion, or anomalous transmutation.**

This repository collects reproducible calculations developed around:

1. low-energy D–D fusion / materials-enhanced screening consistency checks;
2. solid vs molten metal hypotheses;
3. Pd–rare-earth alloys (Pd–La, Pd–Ce, Pd–Y, Pd–Nd, Pd–Gd);
4. isotope reaction / decay networks and element-to-element transmutation;
5. a causal / multiway graph layer inspired by `Shtenco/info_graph_theory`;
6. replacement of hand-entered nuclear constants by evaluated nuclear-data and transport toolchains.

## Scientific rule

Every quantity should be labelled as one of:

- **EVALUATED / MEASURED** — imported from an external evaluated dataset or experiment;
- **DERIVED** — follows from stated equations and inputs;
- **STRESS TEST** — deliberate extrapolation used to expose sensitivity;
- **HYPOTHESIS** — phenomenological parameter or unvalidated mechanism;
- **LEGACY** — historical calculation retained for audit, not endorsed as physical truth.

A numerical model is never treated as evidence that a nuclear process occurs in nature.

## Repository layout

```text
models/
  transmutation_phase_audit.py
  lenr_master_core.py
  molten_alloy_core.py

results/
  master/
  molten/
  transmutation/

research/
  COMPUTATIONAL_STACK.md
  GRAPH_ARCHITECTURE.md
  SCIENTIFIC_BOUNDARIES.md
```

## Current findings

### D–D / screening

The master audit reproduces a Bosch–Hale-style D–D baseline and separately evaluates a **naive low-temperature screening extrapolation**. The latter is explicitly marked invalid as a reactor prediction; it is retained only as a sensitivity/stress test.

### Molten metals

The molten branch does **not** copy crystalline defects into the liquid. It replaces them with D loading, short-range order, electronic response, local coordination and an active space-time fraction.

The main unresolved variable is whether strong effective low-energy reaction enhancement survives the solid → liquid transition. This is not assumed.

### Transmutation

The transmutation audit separates reaction identity / radioactive decay from the source strength of nuclear events and from phase-dependent materials effects.

A phase transition does not by itself change nuclear identity rules. Any discontinuity in an A→B transmutation rate across liquidus must enter through source strength, spectrum, local particle transport or another measurable state-dependent quantity.

For the first natural-isotope 50:50 Pd–rare-earth audit, the most diagnostically useful systems were Pd–La, Pd–Y and Pd–Ce. Pd–Gd is an extremely strong absorber but most captures remain within Gd isotopes rather than changing chemical element.

## Reproducibility policy

Large Monte-Carlo tables are **generated, not hand-edited**. The repository stores source code, seeds, scan bounds, compact CSV/JSON summaries and run logs. Very large deterministic raw tables may be omitted from Git history when they can be regenerated exactly from committed code and seed.

## Target computational stack

```text
CALPHAD / phase thermodynamics
          ↓
DFT / electronic structure
          ↓
MD / liquid short-range structure and D distribution
          ↓
material-state graph
          ↓
ENDF/GNDS nuclear data → transport / reaction rates
          ↓
decay network / isotope inventory
          ↓
causal multiway verification graph
```

Candidate open-source backends under evaluation: OpenMC, LLNL FUDGE, NJOY, radioactivedecay, pycalphad, Quantum ESPRESSO, LAMMPS, DeePMD/DP-GEN, pymatgen and pynucastro.

## Initial consolidated calculations

- D–D Bosch–Hale/Gamow consistency table;
- screening enhancement sensitivity table;
- branch-A / branch-B phenomenological summaries;
- neutron-budget consistency table;
- corrected one-group Hg→Au / Ir→Pt / Ru→Rh audit;
- molten precious × rare-earth alloy ranking;
- Pd–La / Pd–Y / Pd–Ce / Pd–Nd / Pd–Gd natural-isotope transmutation ranking;
- below/above-liquidus phase sensitivity table.

Consolidation date: **2026-08-08**.
