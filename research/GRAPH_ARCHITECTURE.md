# Nuclear-material multiway graph architecture

Date: 2026-08-08

This document defines how ideas from `Shtenco/info_graph_theory` are reused without claiming that CIMFIG changes nuclear physics.

## State

A model state is represented as

```text
Omega = (
  nuclide_inventory,
  alloy_composition,
  phase_fractions,
  temperature,
  pressure,
  D_loading,
  local_structure_descriptors,
  electronic_descriptors,
  particle_spectrum
)
```

## Edge classes

### A. Evaluated nuclear edges

Examples:

```text
139La + n -> 140La
140La -> 140Ce + beta-
89Y + n -> 90Y
90Y -> 90Zr + beta-
108Pd + n -> 109Pd
109Pd -> 109Ag + beta-
```

Weights come from evaluated nuclear data and decay data.

### B. Material-state edges

Examples:

```text
solid -> solid+liquid
solid+liquid -> liquid
D_bulk -> D_interface
local_cluster_A -> local_cluster_B
```

Weights/probabilities come from CALPHAD / DFT / MD, not nuclear tables.

### C. Hypothesis edges

These are allowed only when explicitly labelled.

Example:

```text
local_metal_state -> anomalous_nuclear_source
```

No free edge is allowed to masquerade as evaluated physics. Its strength must be inferred from data or scanned as a falsifiable parameter.

## Two linked graphs

The project uses two graphs rather than one overloaded graph.

### Materials graph

```text
composition + T + D
        ↓
solid / interface / liquid
        ↓
electronic structure
        ↓
D-D local statistics
```

### Nuclear graph

```text
nuclide + particle spectrum
        ↓
reaction products
        ↓
radioactive intermediates
        ↓
stable / long-lived products
```

The graphs interact only through measurable coupling variables such as particle source spectrum, local number density, and explicitly tested material-state descriptors.

## Multiway histories

A sample may have many competing histories:

```text
Pd-La-D state
  ├─ no nuclear event
  ├─ Pd capture -> Ag chain
  ├─ La capture -> Ce chain
  ├─ D-D branch -> He3+n
  └─ background / contamination path
```

The observable is the final isotope inventory plus its time dependence, not a visually attractive graph path.

## Coarse graining

Microscopic MD configurations are too numerous for the nuclear layer. They are coarse-grained into descriptors such as:

- phase fraction;
- D atomic fraction;
- D-D pair probability below selected radii;
- coordination-number distribution;
- electronic density / DOS descriptors;
- interface fraction;
- diffusion coefficient.

The coarse map must be tested for stability: changing the microscopic sampling should not arbitrarily change the inferred nuclear prediction.

## Verification gates

### G1 — conservation

Mass number, charge and isotope inventory accounting must close for each evaluated reaction path.

### G2 — solver agreement

For a simplified benchmark, Bateman/ODE, graph propagation and transport/depletion solutions must agree within stated tolerance.

### G3 — phase-null test

If the nuclear source is held fixed, crossing liquidus should not create a fictitious discontinuity solely from the graph implementation.

### G4 — delayed fingerprint

A proposed transmutation must predict the correct intermediate half-life and daughter growth curve.

### G5 — background path

The graph must include contamination / impurity / measurement-background alternatives. A new-element observation is not assigned to transmutation unless those alternatives are quantitatively weaker.

### G6 — hypothesis isolation

Removing all HYPOTHESIS edges must reproduce standard nuclear physics. If it does not, the graph implementation is invalid.

## Main experimental-computation targets

### Pd-La

Primary diagnostic pair:

```text
139La(n,g)140La -> 140Ce
108Pd(n,g)109Pd -> 109Ag
```

### Pd-Y

```text
89Y(n,g)90Y -> 90Zr
108Pd(n,g)109Pd -> 109Ag
```

### Pd-Ce

```text
140Ce(n,g)141Ce -> 141Pr
142Ce(n,g)143Ce -> 143Pr -> 143Nd
```

## Key falsifier

If an alleged liquid-phase anomaly is real, the model should predict **both**:

1. a phase-correlated change in source/event rate; and
2. the correct isotope/decay fingerprint afterward.

Heat without the corresponding nuclear products falsifies an ordinary D-D / neutron-capture interpretation at the claimed rate.
