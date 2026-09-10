# Changelog

This file records major computational and scientific milestones of the
Masterarbeit project.

## Unreleased

### Changed

- Porting the validated workflow from Stromboli to SuperMUC-NG.
- Adapting compiler, MPI, HDF5, filesystem paths, and Slurm configuration
  for the SuperMUC-NG environment.
- Preparing the qX ensemble for the next production stage.

## stromboli-baseline-2026-09-10

Validated computational baseline immediately before migration from
Stromboli to SuperMUC-NG.

### Environment

- Stromboli computing cluster
- GCC 10.2.1 toolchain
- OpenMPI 4.1.0-no_ucx
- Intel MKL
- Parallel HDF5 1.14.6 built locally
- Private QCD library
- BDIO

### Validated workflow

- Laplacian-eigenmode generation
- Laplace-static correlator measurements
- static-potential extraction and Cornell fits
- gluonic clover insertion
- tiny open-boundary validation
- temporal-source averaging and production diagnostics
- weighted-profile and GEVP validation

The tag `stromboli-baseline-2026-09-10` identifies the exact tracked
repository state used as the starting point for the SuperMUC-NG port.
