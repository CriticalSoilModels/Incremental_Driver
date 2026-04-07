# Module Summary

## Dependency Flow

Arrows mean "imports from". Leaf nodes at the bottom have no dependencies on project code.

```
app/incrementalDriver
        │
        ├── mod_UMAT (elastic.f90)
        │
        └── mod_run_model
                │
                ├── mod_inc_driver_funcs ──────────────────────┐
                │   (splitaLine, ReadStepCommons, PARSER,      │
                │    get_increment, USOLVER, EXITNOW,          │
                │    ROTSIG, SINV, SPRINC, SPRIND, XIT,       │
                │    inv33, spectral_decomp, LU solver)        │
                │                                              │
                ├── mod_loads ──────────────── mod_alignment ──┤
                │   (13 load readers)              │           │
                │                                 └───────────►mod_inc_driver_funcs
                │                                              │
                ├── mod_file_io                                │
                │   (read params/IC, write output)             │
                │                                              │
                ├── mod_command_line                           │
                │   (CLI args → filenames)                     │
                │                                              │
                ├── mod_step_params ◄──── mod_loads            │
                │   (descriptionOfStep,                        │
                │    get/set_repetition_params)                │
                │                                              │
                ├── mod_types ◄──────────── mod_loads          │
                │   (StressAlignment)    ◄── mod_alignment     │
                │                                              │
                ├── mod_matrices                               │
                │   (MRosc, MCart, MRendul, …)                 │
                │                                              │
                ├── mod_maps                                   │
                │   (map2T/D/stress/stran, write66/6)          │
                │                                              │
                ├── mod_value_checks                           │
                │   (set_zero_with_tol,                        │
                │    check_stress_inc_size)                    │
                │                                              │
                └── mod_constants ◄────────────────────────────┘
                    (voight_len, max_fname_len, …)
                    (used by most modules above)
```

### Leaf modules (no project imports)
- `mod_constants` — only imports `stdlib_kinds`
- `mod_matrices` — only imports `stdlib_kinds`
- `mod_maps` — only imports `stdlib_kinds`
- `mod_UMAT` — only imports `stdlib_kinds`
- `mod_command_line` — no project imports
- `mod_value_checks` — only imports `stdlib_kinds`, `stdlib_optval`

### The coupling problem

`mod_inc_driver_funcs` is imported by three different modules (`mod_run_model`, `mod_loads`, `mod_alignment`). It's acting as a shared utilities layer, which is why it hasn't been split yet — everything depends on it. Breaking it apart (`indr_linalg`, `indr_parser`, `indr_solver`) will untangle this.

---

A description of every module in the project — what it contains, what it owns, and what's missing or worth discussing.

---

## `app/incrementalDriver.f90` — main program

**Role:** UMAT selector + launch point.

```
use mod_UMAT, only: UMAT        ← currently elastic; intended to be UMAT_MCSS
use mod_run_model, only: run_model
call run_model(umat)
```

**Status:** Good. Intentionally thin. The only thing that belongs here is the choice of UMAT.

---

## `mod_run_model` (`run_model.f90`)

**Role:** Top-level driver loop. Owns the entire simulation: reads inputs, iterates over steps and increments, calls UMAT, writes output.

**Public procedures:**
- `run_model(UMAT)` — the full simulation; takes the UMAT as a procedure argument via `umat_interface`

**What's inside `run_model`:**

| Concern | Lines (approx) |
|---------|---------------|
| File setup (filenames, open) | ~15 |
| Outer step loop (`do_keyword`) | whole thing |
| `*Repetition` handling — store/recall `ofStep(:)` | ~20 |
| Load reading dispatch (`if/else if` over 14 keywords) | ~80 |
| Zero-load UMAT call for initial stiffness | ~10 |
| Coordinate system selection (`M`, `MmT`) | ~10 |
| `*ImportFile` data reading (row by row) | ~25 |
| Inner increment loop (`do_kinc`) | ~100 |
| Equilibrium iteration loop (`do_kiter`) | ~50 |
| `*DeformationGradient` stress/strain rotation | ~10 |
| Output writing + time advance | ~15 |

**Local variables:** 50+ flat scalars and arrays — all of the UMAT arguments, all step parameters, all iteration workspace, all file IDs.

**Issues:**
- This subroutine is still a monolith (~430 lines). The flat variable list is the core problem.
- `every`/`ievery` (output frequency counter) is declared but used inconsistently — the write logic uses `ievery==1` but the actual write is a direct `write()` call without using `write_line_output_data`.
- `import_file_id` is closed unconditionally at the end even if `*ImportFile` was never used.
- Still uses bare `use mod_loads` and `use mod_maps` (no `only:` clause).

---

## `mod_loads` (`loads.f90`)

**Role:** One subroutine per load type. Reads from `test.inp` and populates the flat load parameters.

**Public procedures (13 load readers):**

| Subroutine | Load keyword | What it reads |
|-----------|-------------|---------------|
| `read_linear_load` | `*LinearLoad` | ninc, maxiter, dt, dT, freq, keyword, ifstress(6), deltaLoad(6) |
| `read_circulating_load` | `*CirculatingLoad` | + deltaLoadCirc(6), phase(6) |
| `read_deformation_gradient_load` | `*DeformationGradient` | deltaLoad(9) (full F tensor) |
| `read_file_load` | `*ImportFile` | columnsInFile(7), importFactor(7), filename |
| `read_oedometric_load` | `*OedometricE1` | deltaLoad(1) (strain-ctrl oedometer) |
| `read_oedometric_S1_load` | `*OedometricS1` | deltaLoad(1), ifstress(1) (stress-ctrl) |
| `read_triaxial_e1_load` | `*TriaxialE1` | deltaLoad(1) (strain-ctrl triaxial) |
| `read_triaxial_s1_load` | `*TriaxialS1` | deltaLoad(1), ifstress(1:3) |
| `read_triaxial_ueq_load` | `*TriaxialUEq` | deltaLoad(2) (Roscoe q, undrained) |
| `read_triaxial_uq_load` | `*TriaxialUq` | deltaLoad(2), ifstress(2) |
| `read_pure_relaxation_load` | `*PureRelaxation` | (no load, just timing) |
| `read_pure_creep_load` | `*PureCreep` | ifstress = all stress-ctrl |
| `read_undrained_creep` | `*UndrainedCreep` | ifstress(2:6)=1, ifstress(1)=0 |
| `read_obey_restrictions_load` | `*ObeyRestrictions` | 6 constraint lines → cMt, cMe, mb |
| `read_perturbations_S_load` | `*PerturbationsS` | deltaLoad(1), keyword3 |
| `read_perturbations_E_load` | `*PerturbationsE` | deltaLoad(1), keyword3 |
| `read_random_walk_load` | `*RandomWalk` | ifstress(6), deltaLoad(6) |

**Issues:**
- Every subroutine takes the same 5 scalar arguments: `ninc, maxiter, deltaTime, deltaTemp, write_freq`. These are begging to be a type.
- `private`/`public` are commented out — everything is accidentally public.
- `read_linear_load` hardcodes `read(1, ...)` instead of using `test_file_id`.
- `read_random_walk_load` also hardcodes `read(1, ...)`.

---

## `mod_inc_driver_funcs` (`incrementalDriver_funcs.f90`)

**Role:** Legacy grab-bag. Contains four distinct concerns that haven't been extracted yet.

**Public procedures:**

| Procedure | Concern | Notes |
|-----------|---------|-------|
| `get_increment` | Convert step description → per-increment ddstress/dstran | 14 arguments, 8 outputs |
| `USOLVER` | Newton equilibrium iterator (mixed stress/strain) | Contains LU decomposition internally |
| `splitaLine` | String split at separator character | Pure utility |
| `ReadStepCommons` | Read ninc/maxiter/deltaTime/deltaTemp from step header | Reads from file unit |
| `PARSER` | Parse `*ObeyRestrictions` constraint expressions | Reads from character array |
| `EXITNOW` | Evaluate exit condition string against current state | Pure logic |
| `ROTSIG` | Rotate a stress or strain tensor (Abaqus imitation) | |
| `SINV` | Compute mean stress p and deviatoric stress q (Abaqus) | |
| `SPRINC` | Principal values (Abaqus) | |
| `SPRIND` | Principal values + directions (Abaqus) | |
| `XIT` | Stop (Abaqus) | |

**Private (contained in `get_increment` or `USOLVER`):**
- `inv33` — 3×3 matrix inverse (contained in `get_increment`)
- `spectral_decomposition_of_symmetric` — Jacobi iteration eigenvalue solver
- `app_jacobian_similarity` — Givens rotation step
- `get_jacobian_rot` — find optimal rotation parameters
- `ludcmp`, `lubksb`, `xLittleUnsymmetricSolver` — LU decomposition (contained in `USOLVER`)

**Issues:**
- Four distinct concerns in one file: math/linalg, Abaqus API imitations, text parsing, increment/solver logic.
- `get_increment` has 14 arguments (8 outputs) — the outputs are a natural type.
- `inv33` buried in a `contains` block makes it untestable and unreusable.
- Module name `mod_inc_driver_funcs` is a placeholder; it should be split and deleted.

---

## `mod_file_io` (`file_operations.f90`)

**Role:** All file reading/writing for the three input files and the output file.

**Public procedures:**

| Procedure | What it does |
|-----------|-------------|
| `read_parameter_file` | Reads cmname, nprops, props(:) from `parameters.inp` |
| `read_init_conditions_file` | Reads initial stress, nstatv, statev(:), temperature from `initialconditions.inp` |
| `set_output_name_from_test_file` | Reads line 1 of `test.inp` for output filename + optional heading |
| `write_output_file_header` | Writes column header row (time, stran, stress, statev) |
| `write_line_output_data` | Writes one row of time/strain/stress/statev to output |

**Issues:**
- `read_init_conditions_file` returns both `state_vars` and `init_state_vars` — the caller only needs one; the copy is made internally and the distinction leaks out.
- `write_output_file_header` builds the header format internally — tightly coupled to the column layout.

---

## `mod_command_line` (`command_line.f90`)

**Role:** Parse CLI arguments to set input filenames and verbosity.

**Public procedures:**
- `set_inputs(parametersfilename, initialconditionsfilename, testfilename, outputfilename, verbose)` — sets defaults then overrides from command line

**Issues:**
- The 5 output arguments are all file-name strings plus one logical — a natural config type.

---

## `mod_step_params` (`step_params.f90`)

**Role:** The `descriptionOfStep` type + pack/unpack helpers for `*Repetition` support.

**Types:**
- `descriptionOfStep` — stores everything needed to replay a step on the next repetition:
  `ninc`, `maxiter`, `ifstress(6)`, `columnsInFile(7)`, `mImport`,
  `deltaLoadCirc(6)`, `phase0(6)`, `deltaLoad(9)`, `dfgrd0(3,3)`, `dfgrd1(3,3)`,
  `deltaTime`, `importFactor(7)`, `deltaTemp`,
  `keyword2`, `keyword3`, `exitCond`, `ImportFileName`,
  `cMt(6,6)`, `cMe(6,6)`, `mbinc(6)`, `existCond`

**Public procedures:**
- `get_repetition_params` — unpacks a `descriptionOfStep` into ~20 individual variables
- `set_repetition_params` — packs ~20 individual variables into a `descriptionOfStep`

**Issues:**
- The type name `descriptionOfStep` is non-standard style; should be `step_config_t`.
- The pack/unpack functions are a symptom of the flat variable architecture in `run_model` — if `run_model` held a `step_config_t` directly, these helpers would be unnecessary.
- The type bundles step-reading state (deltaLoad, deltaTime) with runtime state (dfgrd0, dfgrd1) — these are two different things.

---

## `mod_types` (`types.f90`)

**Role:** The `StressAlignment` type used by `*ImportFile` loads.

**Types:**
- `StressAlignment` — controls stress alignment from a reference file:
  `active`, `ImportFileName`, `kblank`, `nrec`, `kReversal`, `ncol`,
  `Reversal(100)`, `isig(6)`, `sigFac(6)`

**Issues:**
- Type name `StressAlignment` should be `stress_align_t`.
- `private`/`public` commented out.
- Module contains an empty `contains` block.

---

## `mod_constants` (`constants.f90`)

**Role:** Named integer parameters shared across modules.

**Contents:**
```
max_fname_len  = 40   ! max filename length
max_mater_len  = 80   ! max material name length
max_lname_len  = 40   ! max load keyword length
voight_len     = 6    ! Voigt vector length  ← misspelled (should be voigt_len)
max_head_len   = 260  ! max heading line length
iter_lower_limit = 5  ! minimum iterations for stress-controlled steps
```

**Issues:**
- `voight_len` → `voigt_len` (typo).
- No precision constant (`wp` or `dp`) — each module imports from `stdlib_kinds` directly, which is fine but means changing precision requires touching every file.

---

## `mod_matrices` (`matrices.f90`)

**Role:** The six 6×6 isomorphic transformation matrices as compile-time `parameter` arrays.

**Contents:**

| Constant | Represents |
|----------|-----------|
| `MRoscI`, `MRoscImt` | Isomorphic Roscoe (P, Q, Z variables) |
| `MRendul`, `MRendulmT` | Rendulic |
| `MRosc`, `MRoscmT` | Non-isomorphic Roscoe (p, q, z) |
| `MCart`, `MCartmT` | Cartesian (identity-like) |

Each matrix transforms a Voigt stress/strain vector between coordinate systems. `MmT` denotes the transpose of the inverse (needed because stress and strain transform differently).

**Status:** Clean and well-contained.

---

## `mod_maps` (`maps.f90`)

**Role:** Two unrelated things in one file — Voigt↔tensor conversions and Mathematica debug writers.

**Voigt↔tensor conversions:**

| Function | What it does |
|----------|-------------|
| `map2T(a, ntens)` | Voigt stress vector → 3×3 symmetric tensor |
| `map2stress(a, ntens)` | 3×3 tensor → Voigt stress vector |
| `map2D(a, ntens)` | Voigt strain vector → 3×3 strain tensor (halves shear components) |
| `map2stran(a, ntens)` | 3×3 strain tensor → Voigt strain vector (doubles shear components) |

**Debug writers:**
- `write66(a)` — writes a 6×6 matrix to `nic.m` in Mathematica format (append mode)
- `write6(a)` — writes a 6×1 vector to `nic.m`

**Issues:**
- The Voigt conversions and the debug writers are unrelated and should be in separate modules.
- The `ntens` argument on every conversion function is unnecessary — the functions always use all 6 components internally, and `ntens` just slices the output. Could use assumed-shape instead.
- Debug writers hardcode the file name `nic.m` and unit 12.

---

## `mod_alignment` (`alignment.f90`)

**Role:** Stress alignment against a reference file, used with `*ImportFile` loads.

**Public procedures:**
- `readAlignment(align, ImportFileName)` — reads a `.rev` file describing which stress components to align and by what factors
- `tryAlignStress(align, kinc, newState, mImport, stress, ntens)` — applies alignment to stress at a given increment

**Issues:**
- Depends on `splitaLine` from `mod_inc_driver_funcs` — will need updating when that's extracted.

---

## `mod_value_checks` (`write_formatting.f90`)

**Role:** Utility functions called in the output path of `run_model`.

**Public procedures:**
- `set_zero_with_tol(arr, ext_tol)` — zeros values below a tolerance to prevent underflow during formatted write (e.g. `1.3E-391`)
- `check_stress_inc_size(a_dstress, u_dstress)` — warns if the stress increment is diverging

**Issues:**
- File is misnamed (`write_formatting.f90` → should be `value_checks.f90` to match module name).

---

## `mod_UMAT` (`elastic.f90`)

**Role:** Example linear elastic UMAT satisfying the Abaqus subroutine interface.

**Contents:**
- `UMAT(...)` — computes DDSDDE from E and nu, then updates stress via `sig += DDSDDE * dstran`

**Issues:**
- Still uses `real(8)` throughout instead of `real(dp)`.
- `CHARACTER*80` (old-style).
- Missing `implicit none` at module level.

---

# Types Assessment

## Existing types

| Type | Module | Style name | Role |
|------|--------|-----------|------|
| `descriptionOfStep` | `mod_step_params` | → `step_config_t` | Stores step params for *Repetition |
| `StressAlignment` | `mod_types` | → `stress_align_t` | Controls *ImportFile stress alignment |
| `umat_interface` | `mod_run_model` | (abstract interface) | UMAT procedure signature |

## Types that are missing and would help

### `umat_state_t`
Groups all state that evolves per increment. Currently ~10 separate variables in `run_model`:
```fortran
type :: umat_state_t
   real(dp) :: stress(6)    !! Cauchy stress [kPa], Voigt
   real(dp) :: stran(6)     !! Total strain [-], Voigt
   real(dp), allocatable :: statev(:)  !! Internal state variables
   real(dp) :: time(2)      !! (step time, total time) [s]
   real(dp) :: dtime        !! Increment time [s]
   real(dp) :: temp         !! Temperature [°C]
   real(dp) :: dtemp        !! Temperature increment [°C]
end type
```
**Impact:** Cleans up `run_model` significantly. UMAT calls become `call umat(state%stress, state%statev, ...)` and the save/restore pattern (`r_stress`, `r_statev`) becomes `r_state = state`.

### `umat_config_t`
Groups the dimensional/material constants that don't change during a run:
```fortran
type :: umat_config_t
   integer  :: ntens, ndi, nshr, nstatv, nprops
   integer  :: noel, npt, layer, kspt
   character(len=80) :: cmname
   real(dp), allocatable :: props(:)
end type
```
**Impact:** The boilerplate `ndi, nshr, ntens, nstatv, nprops, noel, npt, layer, kspt, cmname, props` block that appears at every UMAT call collapses to passing `config`.

### `load_step_t`
The 5 common scalar arguments repeated across every load-reading subroutine in `mod_loads`:
```fortran
type :: load_step_t
   integer  :: ninc          !! Number of sub-increments
   integer  :: maxiter       !! Max equilibrium iterations
   real(dp) :: delta_time    !! Total step duration [s]
   real(dp) :: delta_temp    !! Temperature change over step [°C]
   integer  :: write_freq    !! Output every N increments
end type
```
**Impact:** Every subroutine in `mod_loads` drops 5 arguments. The common header-reading call (`ReadStepCommons`) returns a `load_step_t` directly.

### `increment_result_t`
What `get_increment` currently returns as 8 separate `intent(out)` arguments:
```fortran
type :: increment_result_t
   real(dp) :: dtime         !! Increment time step [s]
   real(dp) :: ddstress(6)   !! Prescribed stress increment [kPa], Voigt
   real(dp) :: dstran(6)     !! Prescribed strain increment [-], Voigt
   real(dp) :: dtemp         !! Temperature increment [°C]
   real(dp) :: Qb33(3,3)     !! Rotation matrix (Hughes-Winget)
   real(dp) :: dfgrd0(3,3)   !! Deformation gradient at step start
   real(dp) :: dfgrd1(3,3)   !! Deformation gradient at step end
   real(dp) :: drot(3,3)     !! Incremental rotation tensor
end type
```
**Impact:** `get_increment` becomes a function returning `increment_result_t`. Callers are much cleaner, and the result can be tested as a unit.

### `solver_workspace_t` (optional)
The iteration workspace inside `run_model`'s inner loop:
```fortran
type :: solver_workspace_t
   real(dp) :: r_stress(6)       !! Remembered stress at start of increment
   real(dp), allocatable :: r_statev(:)  !! Remembered statev at start of increment
   real(dp) :: a_dstress(6)      !! Approximated stress increment (Roscoe)
   real(dp) :: u_dstress(6)      !! Undesired stress increment
   real(dp) :: c_dstran(6)       !! Correction strain increment
   real(dp) :: ddsdde_bar(6,6)   !! Transformed stiffness (Roscoe-Roscoe)
end type
```
**Impact:** Isolates the Newton iteration state. Lower priority than the others — mainly useful if the equilibrium loop is extracted into its own subroutine.

---

# Summary of Issues

| Priority | Issue | Where |
|----------|-------|-------|
| High | `run_model` is a 430-line monolith with 50+ flat variables | `run_model.f90` |
| High | All 17 load subroutines repeat the same 5 scalar args | `loads.f90` |
| High | `get_increment` has 8 output arguments | `incrementalDriver_funcs.f90` |
| High | `mod_inc_driver_funcs` mixes 4 unrelated concerns | `incrementalDriver_funcs.f90` |
| Medium | `descriptionOfStep` / `StressAlignment` need renaming | `step_params.f90`, `types.f90` |
| Medium | `real(8)` still used in `elastic.f90` and the `umat_interface` | `run_model.f90`, `elastic.f90` |
| Medium | `voight_len` typo | `constants.f90` |
| Medium | `mod_maps` mixes Voigt conversions with debug writers | `maps.f90` |
| Medium | `read_linear_load` and `read_random_walk_load` hardcode unit `1` | `loads.f90` |
| Low | `private`/`public` commented out in `mod_loads`, `mod_types` | `loads.f90`, `types.f90` |
| Low | `write_formatting.f90` filename doesn't match module name | `write_formatting.f90` |
| Low | Debug writers hardcode filename `nic.m` | `maps.f90` |
