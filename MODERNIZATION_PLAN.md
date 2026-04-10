# Modernization Plan — Incremental Driver

## Current State (2026-04-09)

The codebase has been substantially modernized. All modules use the `indr_` prefix,
the old monolithic `mod_inc_driver_funcs` has been fully split out, and the main
driver loop (`run_model.f90`) has been refactored to use `step_config_t` throughout.
12 test suites pass.

### Module inventory

| Module | File | Status |
|--------|------|--------|
| `indr_run_model` | `run_model.f90` | Done — main driver loop; uses `step_config_t` |
| `indr_loads` | `loads.f90` | Done — all read_* routines populate `step_config_t` directly |
| `indr_file_io` | `file_operations.f90` | Done — file read/write |
| `indr_step_params` | `step_params.f90` | Done — `step_config_t` type (renamed from `descriptionOfStep`) |
| `indr_types` | `types.f90` | Done — `StressAlignment` type |
| `indr_constants` | `constants.f90` | Done — `voigt_len`, `iter_lower_limit`, string lengths |
| `indr_matrices` | `matrices.f90` | Done — isomorphic transform matrices |
| `indr_maps` | `maps.f90` | Done — debug write utilities |
| `indr_alignment` | `alignment.f90` | Done — stress alignment |
| `indr_command_line` | `command_line.f90` | Done — CLI argument parsing |
| `indr_value_checks` | `write_formatting.f90` | Done — `set_zero_with_tol`, `check_stress_inc_size` |
| `indr_parser` | `parser.f90` | Done — `splitaLine`, `ReadStepCommons`, `PARSER`, `EXITNOW` |
| `indr_solver` | `solver.f90` | Done — `USOLVER` Newton equilibrium iterator |
| `indr_linalg` | `linalg.f90` | Done — `inv33`, spectral decomposition, Jacobi helpers |
| `indr_abaqus_utils` | `abaqus_utils.f90` | Done — `ROTSIG`, `SINV`, `SPRINC`, `SPRIND`, `XIT` |
| `indr_umat` | `elastic.f90` | Done — example linear elastic UMAT |

The old `mod_inc_driver_funcs` (`incrementalDriver_funcs.f90`) has been fully extracted
and removed.

### Test suite (all passing)

| Suite | File | Coverage |
|-------|------|----------|
| `test_elastic_umat` | `test/test_elastic_umat.f90` | Volumetric, uniaxial, shear loading |
| `test_get_increment` | `test/test_get_increment.f90` | `*LinearLoad` stress/strain/mixed ctrl |
| `test_step_params` | `test/test_step_params.f90` | `step_config_t` field assignment and copy |
| `test_maps` | `test/test_maps.f90` | Debug write utilities |
| `test_linalg` | `test/test_linalg.f90` | `inv33`, spectral decomposition |
| `test_exit_now` | `test/test_exit_now.f90` | `EXITNOW` (5 cases) |
| `check` | `test/check.f90` | Baseline build confirmed |
| `test_split_line` | `test/test_split_line.f90` | `splitaLine` (4 cases) |
| `test_value_checks` | `test/test_value_checks.f90` | `set_zero_with_tol`, `check_stress_inc_size` |
| `test_usolver` | `test/test_usolver.f90` | Pure stress/strain/mixed control |
| `test_read_loads` | `test/test_read_loads.f90` | All 7 `read_*` routines via `step_config_t` |

### Blocked dependency

`critical-soil-models` (intended UMAT: `UMAT_MCSS`) is blocked because NorSand in that
repo requires `aba_param.inc` from Abaqus. Currently using local `elastic.f90` in the
interim.

---

## Remaining Work

### Near term

- [ ] **Commit the `step_config_t` read-routine refactor** — `loads.f90`,
  `run_model.f90`, and `test_read_loads.f90` are updated and all 12 tests pass,
  but the changes have not yet been committed. Run `fpm test` then commit.

- [ ] **`get_increment` argument count** — `get_increment` currently accepts a
  `step_config_t` but still takes several individual args (`time`, `dtime`,
  `ddstress`, `dstran`, `dTemp`, `Qb33`, `dfgrd0`, `dfgrd1`, `drot`).
  Consider whether these belong in a `umat_state_t` to further reduce the
  argument list. Write tests before changing the signature.

- [ ] **`*ImportFile` branch in `run_model.f90`** — still mutates `config%delta_time`
  inside the kinc loop (line ~322). This is a side-effect on stored config;
  consider whether that is the correct semantic or whether a local `dtime`
  override is cleaner.

### Style sweep (apply per module, not all at once)

- [ ] `real(8)` → `real(dp)` throughout remaining legacy-style code
- [ ] `CHARACTER*80` → `character(len=80)`
- [ ] `enddo` / `endif` → `end do` / `end if`
- [ ] Named `parameter` constants instead of magic numbers
- [ ] FORD `!!` doc comments on all public procedures and types
- [ ] Remove remaining `goto` statements (currently in `parser.f90` — the `555`
  error label in `EXITNOW` and `splitaLine`)

### Phase — Math Verification

- [ ] `get_increment` `*DeformationGradient` branch: cite Hughes & Winget (1980)
- [ ] `SINV`: cite Cambridge p–q definition
- [ ] `spectral_decomp_sym`: cite Kielbasinski (Jacobi algorithm)
- [ ] `USOLVER`: document Newton iteration scheme and convergence criterion

### Phase — stdlib Linear Algebra Integration

Replace hand-rolled linear algebra in `indr_linalg` with stdlib equivalents once
all tests are stable.

Candidates:
- `inv33` → `stdlib_linalg` `inv`
- `spectral_decomposition_of_symmetric` + helpers → `stdlib_linalg` `eigh`
- `xLittleUnsymmetricSolver` (inside `USOLVER`) → `stdlib_linalg` `solve`

**Rule:** one stdlib swap per commit; tests must pass before and after.

### Phase — critical-soil-models Integration

Once `aba_param.inc` issue is resolved upstream:
- [ ] Re-enable `critical-soil-models` in `fpm.toml`
- [ ] Switch `app/incrementalDriver.f90` back to `UMAT => UMAT_MCSS`
- [ ] Add integration test: run a known MCSS path and verify stress response

---

## Modernization Rules (quick reference)

1. **Tests before refactoring** — a failing test after extraction means you broke something
2. **One concern per commit** — small, reviewable diffs
3. **Never change math and style in the same commit**
4. **`fpm test` must pass after every commit**
5. **Rename in a separate commit** from behavioural changes
