# $\nu_\mu$ Inclusive Analysis

## Introduction

The purpose of this analysis is to extract 

## Data

This analysis uses the SBND Gen I production.

## Analysis

1. Create detector systematics using `detsys2.ipynb`
2. Create reweightable systematics (including statistical uncertainties) using `systematics.ipynb`
3. Merge the systematics using `systematics.ipynb` as well

* `interaction_plots.ipynb` creates interaction-level plots broken into signal and background. It can utilize systematics from `systematics.ipynb` to create plots, and overlays data and MC.
* `particle_fm.ipynb` creates particle-level plots and explores flash matching cuts. It also utilizes calorimetric variations for the cuts on track-level variables.
* `unfolding3.ipynb` extracts the cross section and does fake data studies. It depends on the systematics from `systematics.ipynb` and detector systematics from `detsys2.ipynb`.