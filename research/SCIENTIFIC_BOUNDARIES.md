# Scientific boundaries and falsification rules

Date: 2026-08-08

## What the current project does establish computationally

- standard reaction-rate bookkeeping can be reproduced;
- Bosch–Hale / Gamow low-energy D–D consistency can be checked numerically;
- isotope capture + decay networks can be propagated in time;
- solid/liquid material-state hypotheses can be parameterized without treating them as established nuclear mechanisms;
- candidate alloys can be ranked by diagnostic usefulness;
- historical unit errors and circular parameters can be exposed.

## What it does NOT establish

The repository does not establish that:

- cold fusion produces net energy;
- a molten metal creates nuclear reactions;
- a particular effective screening parameter measured in one regime applies to a thermal liquid metal;
- heat without nuclear products is evidence of fusion;
- the graph architecture changes nuclear conservation laws;
- a calculated transmutation path implies that the required particle source exists.

## Required labels

### EVALUATED

Value taken from a recognized evaluated dataset or directly from a cited experiment.

### DERIVED

Value obtained by deterministic calculation from EVALUATED inputs.

### HYPOTHESIS

Unvalidated mechanism or state-dependent coupling.

### STRESS_TEST

Intentional extrapolation beyond validated regime to measure sensitivity.

### LEGACY

Historical project value retained so earlier conclusions can be audited.

## Hard falsifiers for a standard D–D explanation

If the model claims a D–D reaction rate corresponding to measurable heat, then standard D–D branch products must be consistent with that rate. A sufficiently low neutron / proton / tritium / He-3 signature relative to the claimed heat falsifies a standard D–D interpretation at that rate.

## Hard falsifiers for neutron-capture transmutation

For a proposed chain such as

```text
139La(n,g)140La -> 140Ce
```

we require:

1. the appropriate precursor isotope;
2. a physically accounted neutron/event source;
3. daughter growth consistent with the intermediate half-life;
4. mass/charge bookkeeping closure;
5. exclusion of initial contamination at the required level.

## Liquidus test interpretation

Standard nuclear data do not predict that simply changing a metal from solid to liquid changes the identity of allowed nuclear reaction paths.

Therefore a sharp transmutation-rate discontinuity near liquidus, if observed after controls, would have to enter through a changed source term, spectrum, transport, local concentration, or another measurable state-dependent quantity.

A graph-code discontinuity by itself is not evidence.

## Data provenance rule

Every future result table should include enough metadata to reconstruct:

- data source/version;
- isotope composition;
- temperature/phase assumptions;
- spectrum/source assumption;
- random seed when applicable;
- code commit;
- solver tolerances.

## Repository philosophy

The goal is not to make the model produce a positive result. The goal is to create a model that can kill incorrect versions of the hypothesis quickly and leave only parameter regions that survive independent physics checks.
