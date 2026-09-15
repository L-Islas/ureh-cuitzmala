# Bayesian network model specification

This directory contains the probabilistic specification and provenance of the Bayesian network used for the UREH regionalization.

## Operational model files

The following CSV files are the model components read directly by the analysis scripts:

- `cpt_infiltration.csv`
  - Child node: infiltration
  - Parents: soil, land cover, slope and precipitation
  - Parent-state combinations: 180

- `cpt_runoff.csv`
  - Child node: runoff
  - Parents: infiltration and precipitation
  - Parent-state combinations: 9

- `cpt_hydrological_function.csv`
  - Child node: hydrological function
  - Parents: accumulated upstream contribution and runoff
  - Parent-state combinations: 12

- `cpt_ureh.csv`
  - Child node: UREH
  - Parents: hydrological function and land cover
  - Parent-state combinations: 15

- `model_states.csv`
  - Defines the ordered states of all nine Bayesian-network nodes.

All operational CPT rows were checked to ensure that:

1. all required parent-state combinations are present;
2. parent-state combinations are unique;
3. no probabilities are missing;
4. all probabilities fall within [0, 1];
5. probabilities sum to 1 for every parent-state combination;
6. all parent and child states match the network state definitions.

## Source model specification

The directory `source/` contains:

`ureh_bayesian_network_model_specification.xlsx`

This workbook is the original model-specification document used to construct the conditional probability tables. It includes the complete CPTs together with the decision rules, support rules, probability scale, methodological notes and model summary.

The workbook was historically named:

`CPT_completas_red_bayesiana_UREH_revision2.xlsx`

It was renamed in this repository to provide a stable and descriptive public filename.

## Model provenance

The directory `provenance/` contains CSV exports of the information used to construct and document the CPTs:

- `decision_rules.csv`
- `support_rules.csv`
- `probability_scale.csv`
- `methodological_notes.csv`
- `cpt_summary.csv`
- `cpt_infiltration_full.csv`
- `cpt_runoff_full.csv`
- `cpt_hydrological_function_full.csv`
- `cpt_ureh_full.csv`

The `*_full.csv` files preserve not only the probabilities but also fields such as dominant state, functional tendency, causal support and applied rule.

These provenance files document how the operational probability tables were constructed; they are not directly required during Bayesian inference.

## State-name harmonization

In the original hydrological-function CPT, the lowest upstream-contribution state was named:

`sin_agua_disponible`

In the operational model this state is named:

`aporte_muy_bajo`

This is a state-label harmonization only. The corresponding probabilities were not modified.

## Network structure

The Bayesian network contains nine nodes:

- soil
- land cover
- slope
- precipitation
- upstream contribution
- infiltration
- runoff
- hydrological function
- UREH

The directed structure is:

- soil → infiltration
- land cover → infiltration
- slope → infiltration
- precipitation → infiltration
- infiltration → runoff
- precipitation → runoff
- upstream contribution → hydrological function
- runoff → hydrological function
- hydrological function → UREH
- land cover → UREH

Root-node prior probabilities are calculated empirically from the analysis-ready microcatchment dataset during model construction and are therefore not stored as fixed CPT files in this directory.