# Computational stack for COLD NUCLEAR

Date: 2026-08-08

The project is split into independent physics layers. No single home-grown formula is allowed to stand in for all of them.

## 1. Evaluated nuclear data

### LLNL FUDGE
Repository: `llnl/fudge`

Role:
- read/process GNDS nuclear data;
- translate ENDF-6 ↔ GNDS;
- provide evaluated reaction data and covariance information to the nuclear layer.

Use in this project:
- replace hand-entered reaction tables where possible;
- generate isotope/reaction edges programmatically;
- retain provenance for every cross section.

## 2. Particle transport / activation

### OpenMC
Repository: `openmc-dev/openmc`

Role:
- continuous-energy Monte-Carlo particle transport;
- spectrum-aware reaction rates;
- depletion / activation workflows.

Use in this project:
- replace one-group `R = N Phi sigma` estimates by energy-integrated rates;
- test self-shielding and resonance effects, especially for Gd and Nd;
- compute isotope inventories from a specified external source spectrum.

Important boundary: OpenMC can calculate consequences of a specified particle source. It does not prove that a hypothetical LENR source exists.

## 3. Radioactive decay

### radioactivedecay
Repository: `radioactivedecay/radioactivedecay`

Role:
- radioactive decay chains;
- branching;
- metastable states;
- time-dependent inventories.

Use in this project:
- convert activation products into delayed chemical-element signatures;
- generate predicted time fingerprints such as `La-140 -> Ce-140` and `Pd-109 -> Ag-109`.

## 4. Reaction-network graph

### pynucastro
Repository: `pynucastro/pynucastro`

Role:
- nuclear reaction network construction and ODE generation;
- graph representation of reaction networks.

Use in this project:
- bridge nuclear reaction data into graph/network diagnostics;
- compare graph evolution with independent Bateman/transport solutions.

## 5. Phase thermodynamics

### pycalphad
Repository: `pycalphad/pycalphad`

Role:
- CALPHAD thermodynamics;
- phase equilibrium calculations;
- solid / liquid / multiphase regions.

Use in this project:
- stop treating `T_liquidus` as a single guessed threshold;
- identify solid, liquid, and solid+liquid windows in Pd–RE alloys;
- generate phase fractions used by the material-state graph.

## 6. Electronic structure

### Quantum ESPRESSO
Repository: `QEF/q-e`

Role:
- plane-wave DFT / electronic-structure calculations.

Use in this project:
- calculate electronic density and local environments in Pd–RE–D snapshots;
- test whether any assumed effective screening proxy correlates with actual electronic observables;
- replace free `Ue` parameters with measurable/calculable descriptors where possible.

Important boundary: DFT electronic screening is not automatically equivalent to the phenomenological nuclear `Ue` used in low-energy fusion fits.

## 7. Atomistic liquid dynamics

### LAMMPS
Repository: `lammps/lammps`

Role:
- classical molecular dynamics;
- liquid structure and diffusion;
- radial distribution functions and local coordination.

Use in this project:
- sample solid / interface / liquid configurations;
- estimate D clustering, residence times, pair-distance distributions and short-range order.

## 8. Machine-learned interatomic potentials

### DeePMD-kit
Repository: `deepmodeling/deepmd-kit`

Role:
- machine-learned interatomic potentials.

Use in this project:
- scale first-principles-quality Pd–RE–D sampling to large atom counts and long MD trajectories once a validated training dataset exists.

## 9. Causal / multiway verification layer

### info_graph_theory
Repository: `Shtenco/info_graph_theory`

Role here is **not** to alter nuclear laws.

We reuse its architectural ideas:
- finite graph/hypergraph state;
- local rewrite rules;
- multiway competing histories;
- coarse graining;
- independent consistency gates;
- explicit falsifiers.

The graph layer consumes outputs from the physics codes above. It must not invent cross sections, decay constants, or nuclear energy releases.

## Target pipeline

```text
pycalphad
   ↓ phase fractions / liquidus window
Q-E / DFT
   ↓ electronic descriptors
LAMMPS / DeePMD
   ↓ local atomic configurations, D-D pair statistics
material-state graph
   ↓ state-dependent candidate source term (hypothesis only)
FUDGE nuclear data
   ↓
OpenMC transport / activation
   ↓ isotope production
radioactivedecay / Bateman network
   ↓ delayed products
pynucastro / NetworkX
   ↓
CIMFIG-style causal verification graph
```

## Rule for acceptance

A claimed anomaly must survive at least four independent checks:

1. isotope and charge/mass bookkeeping;
2. spectrum-aware reaction/transport calculation;
3. time-dependent decay fingerprint;
4. control/background model.

Only after those pass should a materials-phase correlation be discussed as evidence for a new source mechanism.
