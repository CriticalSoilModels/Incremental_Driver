# Modernization Plan — Incremental Driver

## Current State (as of architecture-docs branch)

Significant modularization has already been done. The monolithic original has been split into focused modules. What remains is style cleanup, extracting the remaining legacy code from `mod_inc_driver_funcs`, and building a proper test suite.

### Module inventory

| Module | File | Status |
|--------|------|--------|
| `mod_run_model` | `run_model.f90` | Done — main driver loop |
| `mod_loads` | `loads.f90` | Done — load-step parsing |
| `mod_file_io` | `file_operations.f90` | Done — file read/write |
| `mod_step_params` | `step_params.f90` | Done — `descriptionOfStep` type |
| `mod_types` | `types.f90` | Done — `StressAlignment` type |
| `mod_constants` | `constants.f90` | Done — string lengths, `voigt_len`, etc. |
| `mod_matrices` | `matrices.f90` | Done — isomorphic transform matrices |
| `mod_maps` | `maps.f90` | Done — debug write utilities |
| `mod_alignment` | `alignment.f90` | Done — stress alignment |
| `mod_command_line` | `command_line.f90` | Done — CLI argument parsing |
| `mod_value_checks` | `write_formatting.f90` | Done — value checks/formatting |
| `mod_UMAT` | `elastic.f90` | Done — example elastic UMAT |
| `mod_inc_driver_funcs` | `incrementalDriver_funcs.f90` | **Partially done** — still a grab-bag; see below |

### What remains in `mod_inc_driver_funcs`

| Procedure | Concern | Target module |
|-----------|---------|---------------|
| `get_increment` | Increment computation (load path logic) | `mod_loads` or `mod_step` |
| `USOLVER` | Newton equilibrium iterator | `mod_solver` |
| `EXITNOW` | Exit condition evaluator | `mod_parser` |
| `PARSER` | `*ObeyRestrictions` parser | `mod_parser` |
| `splitaLine` | String split utility | `mod_parser` |
| `ReadStepCommons` | Step header reader | `mod_parser` |
| `ROTSIG` | Abaqus tensor rotation utility | `mod_abaqus_utils` |
| `SINV` | Abaqus stress invariants | `mod_abaqus_utils` |
| `SPRINC` | Abaqus principal values | `mod_abaqus_utils` |
| `SPRIND` | Abaqus principal values + directions | `mod_abaqus_utils` |
| `XIT` | Abaqus stop utility | `mod_abaqus_utils` |
| `inv33` | 3×3 matrix inverse (contained in `get_increment`) | `mod_linalg` |
| `spectral_decomposition_of_symmetric` | Jacobi eigenvalue solver | `mod_linalg` |
| `app_jacobian_similarity` | Givens rotation step | `mod_linalg` |
| `get_jacobian_rot` | Optimal rotation params | `mod_linalg` |

### Blocked dependency

`critical-soil-models` (intended UMAT: `UMAT_MCSS`) is blocked because NorSand in that repo requires `aba_param.inc` from Abaqus. Currently using local `elastic.f90`. Revisit once resolved upstream.

---

## Remaining Phases

### Phase 1 — Safety Net (tests)

**Rule: write tests before touching any module.**

- [x] `test_split_line` — `splitaLine` (4 cases)
- [x] `test_exit_now` — `EXITNOW` (5 cases)
- [x] `test_get_increment` — `*LinearLoad` branch (stress-ctrl, strain-ctrl, mixed)
- [x] `test_elastic_umat` — elastic UMAT (volumetric, uniaxial, shear)

Still needed before extracting from `mod_inc_driver_funcs`:
- [ ] `test_linalg` — `inv33` (identity, diagonal), `spectral_decomp` (diagonal, 2×2 block)
- [ ] `test_abaqus_utils` — `SINV` (hydrostatic, deviatoric), `ROTSIG` (identity rotation)
- [ ] `test_usolver` — pure stress control, pure strain control, mixed

---

### Phase 2 — Extract remaining concerns from `mod_inc_driver_funcs`

Each step: extract → update callers → `fpm test` → commit.

#### 2a. `mod_linalg`

Extract:
- `spectral_decomposition_of_symmetric`
- `app_jacobian_similarity`
- `get_jacobian_rot`
- `inv33` (move out of `contains` block in `get_increment`)

#### 2b. `mod_abaqus_utils`

Extract:
- `ROTSIG`, `SINV`, `SPRINC`, `SPRIND`, `XIT`

Keep original Abaqus names — UMATs call these directly.
Depends on `mod_linalg` for spectral decomposition.

#### 2c. `mod_parser`

Extract:
- `splitaLine` → `split_line`
- `ReadStepCommons` → `read_step_header`
- `PARSER`
- `EXITNOW` → `exit_now`

#### 2d. `mod_solver`

Extract:
- `USOLVER`

#### 2e. Merge `get_increment` into `mod_loads`

`get_increment` belongs with the other load logic already in `mod_loads`.

After 2a–2e, `mod_inc_driver_funcs` should be empty and removable.

---

### Phase 3 — Style Sweep

Apply to each module after it is extracted or revisited (do not batch all at once):

- [ ] `real(8)` → `real(dp)` (use `stdlib_kinds`)
- [ ] `real*8` → `real(dp)`
- [ ] `CHARACTER*80` → `character(len=80)`
- [ ] `integer(4)` → `integer`
- [ ] `enddo` → `end do`, `endif` → `end if`
- [ ] `parameter(ntens=6, ...)` block syntax → `integer, parameter :: ntens = 6`
- [ ] All arguments have explicit `intent`
- [ ] No magic numbers — introduce named `parameter` constants from `mod_constants`
- [ ] Add FORD `!!` doc comments to all public procedures and types
- [ ] `goto` → structured control flow (start with the numbered `do` label in the main program)

---

### Phase 4 — Math Verification

- [ ] `get_increment` `*DeformationGradient` branch: cite Hughes & Winget (1980)
- [ ] `SINV`: cite Cambridge p–q definition
- [ ] `spectral_decomp_sym`: cite Kielbasinski (Jacobi algorithm)
- [ ] `USOLVER`: document Newton iteration scheme and convergence criterion

---

### Phase 5 — critical-soil-models Integration

Once `aba_param.inc` issue is resolved upstream:
- [ ] Re-enable `critical-soil-models` in `fpm.toml`
- [ ] Switch `app/incrementalDriver.f90` back to `UMAT => UMAT_MCSS`
- [ ] Add integration test: run a known MCSS path and verify stress response

---

## Modernization Rules (quick reference)

1. **Tests before refactoring** — a failing test after extraction means you broke something
2. **One module per commit** — small, reviewable diffs
3. **Never change math and style in the same commit**
4. **`fpm test` must pass after every commit**
5. **Rename in a separate commit** from behavioural changes
