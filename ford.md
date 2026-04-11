---
project: Incremental Driver
summary: A Modern Fortran driver for calling and testing geotechnical constitutive
         models (Abaqus-style UMATs). Applies prescribed stress/strain paths to a
         UMAT and records the response — decoupled from any FEA solver so models
         can be calibrated and validated independently.
author: Jonathan Moore
author_email: jonathanm@vt.edu
project_url: https://github.com/CriticalSoilModels/Incremental_Driver
project_github: https://github.com/CriticalSoilModels/Incremental_Driver
src_dir: ./src
output_dir: ./build/doc
license: bsd
year: 2024
source: true
graph: true
coloured_edges: true
proc_internals: true
extra_vartypes: dp
---

A Modern Fortran incremental driver for calling and testing geotechnical
constitutive models (UMATs). Originally written by Andrzej Niemunis
([soilmodels.com/idriver](https://soilmodels.com/idriver/)). This repository
is an ongoing modernisation effort: refactoring legacy Fortran into modular,
testable, modern Fortran (2008+) while preserving mathematical correctness.

## Key modules

- **indr_constants** — named constants (string lengths, `voigt_len`, etc.)
- **indr_types** — core derived types: `step_config_t`, `material_state_t`
- **indr_loads** — parse each load-step type from `test.inp`
- **indr_file_io** — read parameters/initial-conditions files, write output
- **indr_run_model** — top-level driver loop (steps, increments, output)
- **indr_matrices** — isomorphic transform matrices (Roscoe, Rendulic, Cartesian)
- **indr_linalg** — linear algebra utilities
- **indr_solver** — mixed stress/strain equilibrium iterator (Newton loop)
