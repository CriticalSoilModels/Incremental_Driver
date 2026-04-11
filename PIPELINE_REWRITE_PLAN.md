# Plan: Full Input Pipeline Rewrite + Model Abstraction

## Context

Two goals:
1. **File I/O is optional.** The entire incremental driver pipeline must work without
   any files. File I/O is retained for backwards compatibility but is not required.
2. **Support multiple model interfaces.** The driver must work with both Abaqus-style
   UMATs and modern object-oriented models (e.g. `euler_substep` from critical-soil-models).

The core driver responsibility — solving mixed stress/strain boundary conditions via
a Newton loop (USOLVER) — is independent of both file I/O and model type. The model
only needs to return updated stress/state and a stiffness matrix (`ddsdde`). In
practice this is always the elastic stiffness (possibly stress-path dependent), not
the full consistent elastoplastic tangent. USOLVER compensates via iteration.

---

## Architecture Overview

```
┌──────────────────────────────────────────────────────┐
│  File I/O layer  (optional, backwards-compatible)    │
│  parse_test_file → run_model → write_step_output     │
└───────────────────────────┬──────────────────────────┘
                            │ step_config_t, material_state_t
┌───────────────────────────▼──────────────────────────┐
│  Core computation layer  (no files)                  │
│  integrate_step(config, state, runner, results)      │
│  ├── get_increment   load → deps / dsig              │
│  └── USOLVER         Newton loop, mixed BCs          │
└───────────────────────────┬──────────────────────────┘
                            │ class(model_runner_t)
┌───────────────────────────▼──────────────────────────┐
│  Model abstraction                                   │
│  model_runner_t  (abstract)                          │
│  ├── umat_runner_t   wraps Abaqus UMAT               │
│  └── mcss_runner_t   wraps euler_substep             │
└──────────────────────────────────────────────────────┘
```

---

## New Types

### `material_state_t` — persistent model state across increments
Add to `src/step_params.f90`. Serves double duty as the live computation state
and as the per-increment result snapshot (combined from the former `increment_result_t`):

```fortran
type material_state_t
   real(dp)              :: sig(6)        !! Cauchy stress [kPa], Voigt
   real(dp)              :: eps(6)        !! total strain [-]
   real(dp), allocatable :: statev(:)     !! model state variables
   real(dp)              :: time(2)       !! [step time, total time] [s]
   real(dp)              :: dt            !! time increment for this snapshot [s]
   real(dp)              :: temp          !! temperature [°C]
   real(dp)              :: F_start(3,3)  !! deformation gradient, start of increment [-]
   real(dp)              :: F_end(3,3)    !! deformation gradient, end of increment [-]
end type material_state_t
```

`drot` is NOT a field — it is a local computed each increment by `get_increment`.
UMAT scratch outputs (`sse`, `spd`, `ddsig_by_ddeps`, etc.) remain locals inside `integrate_step`.

`integrate_step` returns `results(:) type(material_state_t)` — one snapshot per
increment. The caller (file writer, plotter, test) decides what to do with it.

### `step_record_t` — step config + output frequency
Defined in `src/test_parser.f90` (new):

```fortran
type step_record_t
   type(step_config_t) :: config
   integer             :: write_freq
end type step_record_t
```

`write_freq` is an output concern; it does not belong in `step_config_t`.

---

## Model Abstraction (`src/model_runner.f90`, new)

```fortran
module indr_model_runner
   use stdlib_kinds, only: dp
   use indr_step_params, only: material_state_t
   implicit none(type, external)
   private
   public :: model_runner_t

   type, abstract :: model_runner_t
   contains
      procedure(model_run_i), deferred :: run
   end type model_runner_t

   abstract interface
      subroutine model_run_i(this, state, deps, ddsig_by_ddeps, ndi, nshr, ntens)
         !! Integrate the model one increment.
         !! Returns updated state%sig, state%statev, and ddsig_by_ddeps (∂σ/∂ε).
         !! In practice this is the elastic stiffness (possibly stress-updated)
         !! since few soil models return the full consistent elastoplastic tangent.
         import model_runner_t, material_state_t, dp
         class(model_runner_t),  intent(inout) :: this
         type(material_state_t), intent(inout) :: state
         real(dp), intent(in)  :: deps(ntens)
         real(dp), intent(out) :: ddsig_by_ddeps(ntens, ntens)
         integer,  intent(in)  :: ndi, nshr, ntens
      end subroutine model_run_i
   end interface

end module indr_model_runner
```

### `umat_runner_t` — Abaqus UMAT wrapper
New file `src/umat_runner.f90`:

```fortran
type, extends(model_runner_t) :: umat_runner_t
   procedure(umat_interface), pointer, nopass :: proc => null()
   real(dp), allocatable :: props(:)
   character(80) :: cmname
   integer :: nstatv
contains
   procedure :: run => run_umat
end type umat_runner_t
```

`run_umat` unpacks `state` into the UMAT argument list, calls `this%proc(...)`,
and packs the outputs back into `state`. Fixed UMAT arguments (`noel=1`, `npt=1`,
`layer=1`, `kspt=1`, `coords=[0,0,0]`, `celent=1`, `predef`, `dpred`) are
constants set inside `run_umat`.

### `mcss_runner_t` — critical-soil-models wrapper
New file `src/mcss_runner.f90` (or in the `critical-soil-models` package):

```fortran
type, extends(model_runner_t) :: mcss_runner_t
   type(mcss_model_t)        :: model
   type(integrator_params_t) :: params
contains
   procedure :: run => run_mcss
end type mcss_runner_t
```

```fortran
subroutine run_mcss(this, state, deps, ddsig_by_ddeps, ndi, nshr, ntens)
   call euler_substep(this%model, state%sig, deps, this%params)
   ddsig_by_ddeps = this%model%ddsig_by_ddeps
end subroutine run_mcss
```

Note: `mcss_runner_t` may live in the `critical-soil-models` package since it
depends on types from that library. The abstract `model_runner_t` is what this
package exports; concrete implementations can live anywhere.

---

## Core Computation Layer (`src/step_runner.f90`, new)

```fortran
subroutine integrate_step(config, state, runner, results, align)
   type(step_config_t),                   intent(in)    :: config
   type(material_state_t),                intent(inout) :: state
   class(model_runner_t),                 intent(inout) :: runner
   type(material_state_t), allocatable,   intent(out)   :: results(:)
   type(StressAlignment), optional,       intent(in)    :: align
```

Five arguments (within the 6-arg guideline). Contains the extracted `do_kinc` loop.
No file I/O. No output writing. Returns results array for caller to handle.

**Local variable mapping (from current `run_model`):**

| Old local | New location |
|---|---|
| `stress`, `stran`, `statev`, `time`, `temp` | `state%*` |
| `dfgrd0`, `dfgrd1` | `state%*` |
| `dtime`, `drot`, `Qb33` | locals in `integrate_step` |
| `maxiter` | derived from `config%max_iter` + `config%ifstress` |
| `M`, `MmT` | locals, derived from `config%coord_sys` |
| `r_stress`, `r_statev` | `sig_saved`, `statev_saved` — local rollback buffers |
| UMAT scratch (`sse`, `spd`, `ddsig_by_ddeps`, ...) | locals |
| `import_file_id` | local; open/close within `integrate_step` |

The `umat_interface` abstract interface moves from `indr_run_model` to
`indr_step_params` (prerequisite — Phase 1 below).

---

## File I/O Layer

### `parse_test_file` — `src/test_parser.f90`, module `indr_test_parser`

```fortran
subroutine parse_test_file(file_id, heading, steps, n_steps, align)
   integer,                          intent(in)  :: file_id
   character(len=*),                 intent(out) :: heading
   type(step_record_t), allocatable, intent(out) :: steps(:)
   integer,                          intent(out) :: n_steps
   type(StressAlignment),            intent(out) :: align
```

`*Repetition nSteps nRepetitions` is **expanded upfront**: configs read once,
stored `nRepetitions` times. Buffer with a fixed local array (1000 entries),
then `allocate(steps(n_steps)); steps = buf(1:n_steps)`.

### `write_step_output` — add to `src/file_operations.f90`

```fortran
subroutine write_step_output(file_id, results, write_freq)
   integer,                 intent(in) :: file_id
   type(material_state_t),  intent(in) :: results(:)
   integer,                 intent(in) :: write_freq
```

### `run_model` — rewritten as thin orchestrator (~50 lines)

```fortran
call parse_test_file(test_file_id, heading, steps, n_steps, align)
! ... setup state from init conditions file ...
do kstep = 1, n_steps
   call integrate_step(steps(kstep)%config, state, runner, results, align)
   call write_step_output(output_file_id, results, steps(kstep)%write_freq)
end do
```

---

## Variable Naming Reference

Follows the critical-soil-models naming convention. Apply when writing new code
and when extracting locals into `integrate_step`.

| Old name | New name | Notes |
|---|---|---|
| `stress` | `sig` | Cauchy stress, Voigt |
| `stran` | `eps` | Total strain, Voigt |
| `dstran` / `DSTRAN` | `deps` | Strain increment |
| `ddstress` | `dsig` | Prescribed stress increment |
| `ddsdde` | `ddsig_by_ddeps` | ∂σ/∂ε — tangent from model |
| `ddsdde_bar` | `ddsig_by_ddeps_bar` | Transformed tangent (Roscoe / with constraints) |
| `dfgrd0` | `F_start` | Deformation gradient, start of increment |
| `dfgrd1` | `F_end` | Deformation gradient, end of increment |
| `delta_time` / `dtime` | `dt` | Time increment |
| `delta_temp` / `dtemp` / `dTemp` | `dtemp` | Temperature increment |
| `c_dstran` | `deps_corr` | Strain correction from USOLVER |
| `dstran_Cart` | `deps_Cart` | Strain increment in Cartesian coords |
| `a_dstress` | `dsig_approx` | Running Newton approximation of stress increment |
| `u_dstress` | `dsig_target` | Residual stress USOLVER needs to correct |
| `r_stress` | `sig_saved` | Rollback buffer (stress) |
| `r_statev` | `statev_saved` | Rollback buffer (state variables) |
| `stress_Rosc` | `sig_Rosc` | Stress in Roscoe coordinates |
| `r_stress_Rosc` | `sig_Rosc_saved` | Rollback buffer (Roscoe stress) |
| `Qb33` | `R_polar` | Polar decomposition rotation (DeformationGradient) |
| `drot` | `drot` | Keep — Abaqus convention, well known |
| `M` | `M` | Coordinate transform matrix (keep short) |
| `MmT` | `M_inv_T` | M inverse-transpose |

---

## Phase Sequence

| # | Work | New/changed files | Tests after |
|---|---|---|---|
| 0 | Add `material_state_t`, `increment_result_t` to `step_params.f90`; write `test_umat_state.f90` | `src/step_params.f90`, `test/test_umat_state.f90` | 13 pass |
| 1 | Move `umat_interface` from `run_model` to `step_params` | `src/step_params.f90`, `src/run_model.f90` | 13 pass |
| 2 | Create `src/model_runner.f90` with `model_runner_t` abstract type | `src/model_runner.f90` (new) | 13 pass |
| 3 | Create `src/umat_runner.f90` with `umat_runner_t`; write `test_umat_runner.f90` | `src/umat_runner.f90` (new), `test/test_umat_runner.f90` | 14 pass |
| 4a | Create stub `integrate_step` in `src/step_runner.f90` | `src/step_runner.f90` (new) | 14 pass |
| 4b | Implement `integrate_step` with extracted `do_kinc` logic; write `test_integrate_step.f90` (no file I/O) | `src/step_runner.f90`, `test/test_integrate_step.f90` | 15 pass |
| 5 | Wire `integrate_step` into `run_model` (replace inline `do_kinc`) | `src/run_model.f90` | 15 pass |
| 6 | Create `src/test_parser.f90` with `parse_test_file` (no `*Repetition` yet); write `test_parse_test_file.f90` | `src/test_parser.f90` (new), `test/test_parse_test_file.f90` | 16 pass |
| 7 | Extend `parse_test_file` to expand `*Repetition`; add test case | `src/test_parser.f90`, `test/test_parse_test_file.f90` | 16 pass |
| 8 | Add `write_step_output` to `file_operations.f90` | `src/file_operations.f90` | 16 pass |
| 9 | Rewrite `run_model` as thin orchestrator | `src/run_model.f90` | 16 pass |

---

## Key Test: `test_integrate_step.f90` (Phase 4b)

The proof that file I/O is not required for core computation:

```fortran
! No files — construct everything in memory
config%load_type     = '*LinearLoad'
config%coord_sys     = '*Cartesian'
config%n_inc         = 10
config%ifstress      = 0
config%delta_load(1) = 0.01_dp   ! total deps(1)
config%dt            = 1.0_dp

state%sig     = 0.0_dp;  state%eps = 0.0_dp
allocate(state%statev(1), source=0.0_dp)
state%time    = 0.0_dp;  state%temp = 0.0_dp
state%F_start = identity33;  state%F_end = identity33

runner%proc  = elastic_umat
runner%props = [E, nu]

call integrate_step(config, state, runner, results)

! Check final stress and result count — no file I/O anywhere
call check_real(state%sig(1), expected_sig11, tol, 'sig(1)', nfail)
call check_int(size(results), 10, 'result count', nfail)
```

---

## Tricky Parts

1. **Zero-dstran model call** before `do_kinc`: calls model with `dstran=0, dtime=0`
   for initial stiffness, then rolls back stress/statev. Must be preserved in `integrate_step`.

2. **`*ImportFile` per-increment reading**: open file at top of `integrate_step`
   (guarded by `config%load_type`), skip non-numeric headers, read one row per
   increment, always close after `do_kinc`.

3. **`*DeformationGradient` stress rotation**: after `do_kiter` per increment,
   apply rigid rotation using `Qb33` from `get_increment`.

4. **`*Perturbations*` undo**: after appending to `results`, roll back stran/statev/stress
   so perturbations do not accumulate.

5. **`mcss_runner_t` location**: depends on types from `critical-soil-models`. It
   likely lives in that package, not here. The abstract `model_runner_t` is what
   this package exports; concrete implementations can live anywhere.

---

## Deferred (post-pipeline)

- **csv-fortran** for `*ImportFile` reading and output writing — add after pipeline
  stabilizes; one-line `fpm.toml` dependency.
- **`mcss_runner_t`** implementation — depends on `critical-soil-models` resolving
  the `aba_param.inc` blocker and adding `model%ddsig_by_ddeps`.
- Style sweep (FORD comments, `real(dp)`, named constants) on new files — one
  commit per file after each phase lands.

---

## Verification

After each phase:
```bash
conda run -n fpm fpm build
conda run -n fpm fpm test
```

All existing tests must pass at every commit.
