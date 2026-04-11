# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Modern Fortran incremental driver for calling and testing geotechnical constitutive models (UMATs). The driver applies prescribed stress/strain paths to a UMAT subroutine and records the response — decoupled from any FEA solver so models can be calibrated and validated independently.

Originally written by Andrzej Niemunis (soilmodels.com/idriver/). This repo is an ongoing modernization effort: refactoring legacy Fortran into modular, testable, modern Fortran (2008+) while preserving the mathematical correctness of the original.

## Build and Test Commands

```bash
# Build the project
fpm build

# Run all tests
fpm test

# Run a specific test suite
fpm test -- <suite_name>

# Run the main driver application
fpm run

# Generate documentation (requires ford)
ford ford.md

# Clean build artifacts (preserves dependencies — stdlib is NOT re-downloaded)
fpm clean

# Clean everything including dependencies (stdlib will be re-cloned and rebuilt)
fpm clean --all
```

### Environment Setup
```bash
conda env create -f environment.yml
conda activate fpm
```
Requires: gfortran, fpm (Fortran Package Manager), fortls, ford. `stdlib` is fetched
from GitHub on first build and cached in `build/dependencies/`.
Use `fpm clean` (not `fpm clean --all`) to avoid re-cloning it.

## Modernization Goals

See `MODERNIZATION_PLAN.md` for the full plan. In priority order:

1. **Tests first** — before refactoring any module, write regression tests that capture current behaviour.
2. **Modularise** — split concerns into focused modules (largely done; see Architecture below).
3. **Modern style** — replace `real(8)`, implicit typing, magic numbers, etc.
4. **Math verification** — add inline citations and validate calculations.

## Architecture

```
app/
  incrementalDriver.f90   — main program: imports UMAT and calls run_model

src/
  run_model.f90           — top-level driver loop (steps, increments, output)
  loads.f90               — parse each load-step type from test.inp
  file_operations.f90     — read parameters/initial-conditions files, write output
  step_params.f90         — descriptionOfStep type + helpers
  types.f90               — StressAlignment type
  constants.f90           — named constants (string lengths, voigt_len, etc.)
  matrices.f90            — isomorphic transform matrices (Roscoe, Rendulic, Cartesian)
  maps.f90                — debug write utilities (write66, write6)
  alignment.f90           — stress alignment from file
  command_line.f90        — command-line argument parsing
  write_formatting.f90    — mod_value_checks: set_zero_with_tol, check_stress_inc_size
  incrementalDriver_funcs.f90  — mod_inc_driver_funcs: get_increment, USOLVER, EXITNOW,
                                  splitaLine, PARSER, ReadStepCommons + math utilities
                                  (inv33, spectral decomposition, Abaqus imitations)
  elastic.f90             — mod_UMAT: example linear elastic UMAT

test/
  check.f90               — fpm test entry point
```

Key concepts:
- **Load step** — a block in the input file describing a stress/strain path (`*LinearLoad`, `*CirculatingLoad`, etc.)
- **Increment** — one sub-step within a load step; `get_increment` converts step commands into `ddstress`/`dstran` for that increment
- **UMAT interface** — `mod_UMAT` exposes the standard Abaqus UMAT signature; the main program imports a UMAT and passes it to `run_model`
- **USOLVER** — mixed stress/strain equilibrium iterator (Newton loop)

### Planned UMAT

The intended production UMAT is `UMAT_MCSS` from the `critical-soil-models` repo (Mohr-Coulomb strain softening). Currently blocked because NorSand in that repo requires `aba_param.inc` (Abaqus). Using local `elastic.f90` in the interim.

### Module prefix

All modules in this project use the `indr_` prefix (incremental driver):

```fortran
module indr_loads
module indr_step_params
module indr_file_io
```

The existing `mod_*` names on this branch are legacy and will be renamed to `indr_*` as part of the style sweep. New modules being extracted should use `indr_` from the start.

---

# Modern Fortran Style Guide

A practical style guide for modern Fortran (2008+) applied to this project.

## Naming Conventions

### Files

- Module files: `module_name.f90` (lowercase with underscores)
- Preprocessed files: `module_name.F90` (uppercase extension for files needing preprocessing)
- Test files: `test_module_name.f90`
- One module per file (submodules may be separate)

### Modules

Use the `indr_` prefix to namespace modules and avoid collisions on the fpm registry:

```fortran
! Good — prefixed with project abbreviation
module indr_loads
module indr_step_params

! Bad
module mod_loads            ! Too generic — collides on registry
module ShallowWaterSolver   ! CamelCase
module loads                ! No prefix, collision risk
```

### Derived Types

All derived types use the `_t` suffix with snake_case:

```fortran
type :: step_config_t
type :: umat_state_t
type :: load_path_t
```

### Working Precision

Use `dp` from `stdlib_kinds` consistently:

```fortran
use stdlib_kinds, only: dp

real(dp) :: stress(6)

! Never use literal kind numbers
real(8) :: x           ! Bad — non-portable, forbidden in new code
real(dp) :: x          ! Good
```

### Variable Naming for Mathematical Quantities

Key rules:
- Increments: `d"var"` — `dsig`, `deps`, `dlambda`
- Partial derivatives: `d"num"_by_d"denom"` — `dF_by_dsig`, `dG_by_dsig`
- Stress vector: `sig(6)`, strain increment: `deps(6)`, elastic stiffness: `D_e(6,6)`

### Variables and Procedures

Snake_case with descriptive names:

```fortran
integer :: num_steps
real(dp) :: total_time
subroutine parse_step(line, config)
```

**Exception:** Single-letter variables are acceptable for loop indices (`i`, `j`, `k`) and mathematical formulas matching published notation.

### Function and Subroutine Naming

Verb + noun patterns:

```fortran
subroutine parse_input(line, tokens)
subroutine apply_increment(state, dstran, ddstress)
function get_timestep(step_config) result(dt)

! Logical returns use is_/has_/can_ prefixes
function is_converged(residual, tol) result(converged)
```

### Units in Comments

Document physical units in variable declarations:

```fortran
type :: umat_state_t
   real(dp) :: stress(6)  !! Cauchy stress [kPa], Voigt notation
   real(dp) :: temp       !! Temperature [°C]
end type

real(dp), parameter :: DEFAULT_TOL = 1.0e-6_dp  ! [-]
```

### Constants

UPPERCASE with underscores:

```fortran
real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
real(dp), parameter :: CONV_TOL = 1.0e-8_dp
integer, parameter :: MAX_ITER = 1000
```

---

## Required Practices

### Use Statements with `only` Clause

```fortran
! Good
use stdlib_kinds, only: dp
use indr_step_params, only: descriptionOfStep

! Bad — pollutes namespace, hides dependencies
use mod_step_params
```

### Implicit None

Always in modules and programs:

```fortran
module mod_loads
   implicit none
   private
end module
```

### Intent Declarations

Always declare intent for all procedure arguments:

```fortran
subroutine apply_increment(state, dstran, ddstress)
   type(umat_state_t), intent(inout) :: state
   real(dp),           intent(in)    :: dstran(6)
   real(dp),           intent(in)    :: ddstress(6)
```

### Functions Should Have No Side Effects

```fortran
! Good — pure computation
function voigt_norm(v) result(n)
   real(dp), intent(in) :: v(6)
   real(dp) :: n
   n = sqrt(sum(v**2))
end function

! Bad — function mutates state (use a subroutine instead)
function update_and_return_norm(state) result(n)
   type(umat_state_t), intent(inout) :: state   ! Side effect hidden in function
```

### Private by Default

```fortran
module indr_loads
   implicit none
   private

   public :: read_linear_load
   public :: read_circulating_load
```

### Limit Procedure Arguments

Public procedures should have **6 or fewer arguments**. Group related arguments into derived types:

```fortran
! Bad
subroutine run_step(nx, stress, strain, dt, ninc, tol, max_iter, ...)

! Good
subroutine run_step(state, config)
```

---

## Forbidden Practices

- **No `goto`** — use structured control flow
- **No arithmetic IF** — use `if-then-else` or `select case`
- **No COMMON blocks** — use module variables or derived types
- **No EQUIVALENCE**
- **No fixed-form source** — all files must be `.f90` / `.F90`
- **No assumed-size arrays** (`arr(*)`) — use assumed-shape (`arr(:)`)
- **No `external` statements** — use modules for explicit interfaces
- **No `real(8)` or `real(4)`** — use `real(dp)` with `stdlib_kinds`
- **No implicit save** — avoid module-level variables with initialisation; if needed use explicit `save`

---

## Recommended Practices

### Prefer Allocatable Over Pointer

```fortran
! Good
real(dp), allocatable :: values(:)

! Bad — manual cleanup required
real(dp), pointer :: values(:) => null()
```

### Pure and Elemental Procedures

```fortran
pure function voigt_trace(v) result(tr)
   real(dp), intent(in) :: v(6)
   real(dp) :: tr
   tr = v(1) + v(2) + v(3)
end function
```

### Associate for Readability

```fortran
associate(sig  => state%stress,  &
          deps => state%dstran,  &
          dt   => config%dtime)
   call umat(sig, deps, dt, ...)
end associate
```

### No Magic Numbers

```fortran
real(dp), parameter :: CONV_TOL = 1.0e-8_dp

! Bad
if (res < 1.0e-8) return

! Good
if (res < CONV_TOL) return
```

This applies especially to physical quantities and domain values — any literal whose
meaning is not obvious from context must be given a named parameter with a unit comment:

```fortran
! Bad — what does 0.01 mean? Strain? Stress? Tolerance?
config%delta_load(1) = 0.01_dp
expected = (lam + 2.0_dp*mu) * 0.01_dp

! Good — intent and units are clear
real(dp), parameter :: TOTAL_EPS1 = 0.01_dp  !! total axial strain applied over the step [-]
config%delta_load(1) = TOTAL_EPS1
expected = (lam + 2.0_dp*mu) * TOTAL_EPS1
```

### Avoid Deep Nesting

Maximum 3-4 levels. Use early `cycle` and `return`.

### Documentation (FORD compatible)

```fortran
subroutine read_linear_load(file, config)
   !! Read a *LinearLoad step block from the open test file.
   !!
   !! Populates config with ninc, deltaTime, ifstress, and deltaLoad.
   integer,                   intent(in)  :: file
   type(descriptionOfStep),   intent(out) :: config
```

---

## Testing Strategy

Write tests **before** refactoring. Each test should:

1. Set up the minimal state needed
2. Call the procedure under test
3. Assert the result against a known analytic value

Use `test-drive` (available as a dev-dependency) for new test suites. For regression tests against the legacy code use standalone programs that exit with code 1 on failure.

---

## Common AI/LLM Mistakes in Fortran

### Pi is Not a Built-in Constant

```fortran
use stdlib_constants, only: pi_dp
! or:
real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
```

### `random_number` is a Subroutine

```fortran
call random_number(x)   ! correct
x = random_number()     ! wrong
```

### No I/O in `pure` Procedures

```fortran
pure function compute(x) result(y)
   real(dp), intent(in) :: x
   real(dp) :: y
   y = x**2
   ! print *, x   ! compiler error — I/O is a side effect
end function
```

### Declarations Before Executable Code

```fortran
subroutine foo()
   real(dp) :: x   ! declarations first
   x = 1.0_dp      ! then executable statements
end subroutine
```

---

## Summary Table

| Category | Do | Don't |
|---|---|---|
| Types | `step_config_t`, `umat_state_t` | `StepConfig`, `TModel` |
| Variables | `num_steps`, `conv_tol` | `nSteps`, `tol` |
| Constants | `MAX_ITER`, `CONV_TOL` | `maxIter`, `tol` |
| Imports | `use mod, only: x` | `use mod` |
| Arrays | `arr(:)` | `arr(*)` |
| Precision | `real(dp)` | `real(8)` |
| Memory | `allocatable` | `pointer` (unless needed) |
| Control | `if/select/do` | `goto` |
