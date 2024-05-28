# quickmd-nf

![](https://img.shields.io/badge/current_version-0.0.14-blue)
![](https://github.com/stracquadaniolab/quickmd-nf/workflows/build/badge.svg)
## Overview
A workflow to anperform molecular dynamics (MD) simulations and analysis of proteins in solution

## Configuration

- param1: this is the parameter description (default: "hello")
- param2: this is the parameter description (default: "world")
- ...
- paramN: this is the parameter description (default: "flow")

## Running the workflow

### Install or update the workflow

```bash
nextflow pull stracquadaniolab/quickmd-nf
```

### Run the analysis

```bash
nextflow run stracquadaniolab/quickmd-nf
```

## Results

- `results/analysis.txt`: file with the analysis results.
- `results/tuning.txt`: file with the parameter tuning.
- ...

## Authors

- Josh David Littlefair
