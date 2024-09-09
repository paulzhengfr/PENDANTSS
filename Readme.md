# PENDANTSS: PEnalized Norm-ratios Disentangling Additive Noise, Trend and Sparse Spikes  


This repository is the official implementation in Matlab for the work on *PENDANTSS: PEnalized Norm-ratios Disentangling Additive Noise, Trend and Sparse Spikes* (Submitted). [Paper Link](https://hal.archives-ouvertes.fr/hal-03924136), [Paper Link IEEE](https://ieeexplore.ieee.org/document/10057984), [Poster](./poster.pdf).

## Abstract

Denoising, detrending, deconvolution: usual restoration tasks, traditionally decoupled. Coupled formulations entail complex ill-posed inverse problems. We propose PENDANTSS for joint trend removal and blind deconvolution of sparse peaklike signals. It blends a parsimonious prior with the hypothesis that smooth trend and noise can somewhat be separated by lowpass filtering. We combine the generalized pseudo-norm ratio SOOT/SPOQ sparse penalties $\ell_p/\ell_q$  with the BEADS ternary assisted source separation algorithm. This results in a both convergent and efficient tool, with a novel Trust-Region block alternating variable metric forward-backward approach. It outperforms comparable methods, when applied to typically peaked analytical chemistry signals. Reproducible code is provided in this repo.

## Structure

This repo contains:

- `main_reproduce_1results.m` : Matlab code to reproduce the result for one noise realizations and one case.
- `main_reproduce_all_results.m` : Matlab code to reproduce all the comparison results.

- *data* folder: the data used for confirming the performance.
- *param* folder: hyperparameters used in the codes.
- *functions* folder: All functions used.
  - *Core functions* folder: all the core functions used by PENDANTSS. The core algorithm is `SPOQ_BD_quadraFidel.m`. 
  - *baseline related* folder:  functions related to baseline methods.
  - *param search* folder: functions for hyperparameter search.
-  *results_reproduction_files* folder: The codes to run for either parameter search for different cases or for the comparison simulations of different methods in different cases with different noise sampling numbers to average. It also contains some utility code for treating the raw simulation results, and for plotting the simulation results. The files in this folder are listed as below:
  - `main_run.m` Code to run to reproduce all the comparison results, which uses `main_simulation_BD_*` and `main_treat_simulation_result.m`. 
  - `main_run_simplex_all.m` Code to run to look for hyperparameters for all cases, which uses `main_simplex_BD_*`.
  - Some grid search functions for some parameters `main_run_grid_*` 
  - `main_run_plot_simulations.m` to plot the simulation results.
  - *plot* folder: functions to plot the results.



## Results

The results will be stored in `\path\result\` and after treatment by `\path\results_reproduction_files\main_treat_simulation_result.m`  the statistic of results will be in the folder `\path\result\stat` . 

## Parameters effects
Some more extensive results regarding different hyperparameters of the algorithm is shown in [here](./extensive_results.md).


## To-do

- [ ] Have an end-to-end function to apply the method easily on an other problem.
- [ ] Add comments to all files.
- [ ] Add results with other signals.


### Citing PENDANTSS:
---
```
@ARTICLE{pendantss,
  author={Zheng, Paul and Chouzenoux, Emilie and Duval, Laurent},
  journal={IEEE Signal Processing Letters}, 
  title={{PENDANTSS: PEnalized Norm-Ratios Disentangling Additive Noise, Trend and Sparse Spikes}}, 
  year={2023},
  volume={30},
  number={},
  pages={215-219},
  doi={10.1109/LSP.2023.3251891}}
```
