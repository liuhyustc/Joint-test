# Joint-test

This repository contains R/Shell code for the methods in the paper:

- Liu, H., & Zhang, H. (2025). Powerful Rare-Variant Association Analysis of Secondary Phenotypes. Genetic epidemiology, 49(1), e22589. https://doi.org/10.1002/gepi.22589

It includes both **simulation code** and **real data analysis code**.

## Contents

- `joint_test.R`: core implementation of the proposed joint testing procedures (prospective / retrospective / EB / two-step).
- `pchisqsum2.R`: p-value calculation for mixtures of chi-square distributions (saddlepoint / integration / Liu).
- `Ind_test.R`: independence test between the environmental factor and genotypes (uses `aSPU`).
- `Generate_sample_joint.R`: generate population data under null/alternative for simulations.
- `run_generate_joint.sh`: example script to run `Generate_sample_joint.R` across parameter grids.
- `simu_once_result.R`: run simulation experiments and summarize type I error / power.
- `Joint_estimation.R`: parameter estimation experiments for the joint model.
- `realdata.R`: example workflow for real data analysis (requires external datasets; see below).

## Requirements

- R (>= 3.6 recommended)
- R packages: `Matrix`, `MASS`, `data.table`, `foreach`, `doParallel`, `dplyr`, `CompQuadForm`, `aSPU`

## Quick start

### 1) Simulate data

Generate one population dataset (saved under `joint_test/<null|alter>/`):

```bash
Rscript Generate_sample_joint.R \
  2000000 5 30 1.5 0 1 null
```

Or run a grid of settings:

```bash
bash run_generate_joint.sh
```

### 2) Run simulation evaluation

After generating the required `*.RData` files, run:

```bash
Rscript simu_once_result.R
```

### 3) Real data analysis

`realdata.R` shows how the method was applied to preterm birth data, but it **depends on external files** that are not included in this repository, e.g.

- `PretermData/PD_QC_MCpairs_phenotype.txt`
- `~/PretermData/MGdata.RData`

Please place the corresponding datasets in the expected paths (or edit the paths in `realdata.R`).

## Notes

- The scripts are research code released to support reproducibility of the paper.
- If you use this code, please cite the paper above.
